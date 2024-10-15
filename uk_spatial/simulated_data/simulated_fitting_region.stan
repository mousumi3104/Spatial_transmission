functions {
  array[,] real geometric_random_walk(int M_regions, int N, array[] real init_R, array[,] real rw_noise, real rw_sd) {
    array[N, M_regions] real Rt; 
    for ( m in 1 :M_regions) {
      Rt[1, m] = init_R[m]; // Initialize with the given initial value for each region
      for (i in 2: N) {
        Rt[i, m] = Rt[i - 1, m] + rw_noise[i - 1, m] * rw_sd; // Add noise multiplied by rw_sd
      }
    }
    return exp(Rt); // Return the exponential of the resulting array
  }
}

data {
  int<lower=1> M_regions;
  matrix[M_regions,M_regions] C;
  int<lower=1> final_time;    //number of days to simulate
  int<lower =1> N;
  int<lower=1> initial_seeding_day;
  int<lower=1> data_length;
  // array[data_length,M_regions] int<lower =0> data_deaths;       //number of infected individual at any time
  array[data_length,M_regions] int<lower =0> data_inf;       //number of infected individual at any time
  // array[data_length,M_regions] int<lower =0> data_cases;
  // matrix[final_time, M_regions+1] cases;
  vector[N] SI;
  // vector[N] f1;
  // vector[N] f2;
  vector[M_regions] pop;
  int fitting_start;
  int prediction_horizon;
}

transformed data {
  vector[N] SI_rev; // SI in reverse order
  for(t in 1:N){
    SI_rev[t] = SI[N-t+1];
  }
}

parameters {
  real<lower=0> phi1;
  array[M_regions] real init_R;
  array[N-1,M_regions] real rw_noise; // random walk noise
  real<lower = 0> rw_sd; // random walk standard deviation
}

transformed parameters{
  matrix[N, M_regions] infection = rep_matrix(0, N, M_regions);    // daily initialization
  matrix[N , M_regions] SI_regions = rep_matrix(SI_rev,M_regions);  // serial interval for every region
  array[N, M_regions] real<lower=0> Rt = geometric_random_walk(M_regions,N,init_R, rw_noise, rw_sd);
  
   
    for (m in 1:M_regions){                  // for initial seeding
      infection[1:initial_seeding_day,m] = rep_vector(30,initial_seeding_day);
    }
  
    
  for (t in (initial_seeding_day+1):N){ //for loop over time
   
    row_vector[M_regions] convolution_inf = columns_dot_product(infection[1:(t-1),:] , SI_regions[(N-t+2):N,:]);  //infections at each region "k"
    vector[M_regions] total_inf = C * convolution_inf';     //total infection at region "j"
    vector[M_regions] eff_pop = C * pop;
    
    for (m in 1:M_regions){      // for loop over region "i" (final infection at region "i")   
      real sus = pop[m] - sum(infection[1:(t-1),m]);
      
      row_vector[M_regions] RR ;
      for (mm in 1:M_regions){
        RR[mm] = Rt[t,mm];
      }
      // print(RR);
      infection[t,m] = dot_product(C[:,m]', (((rep_vector(sus , M_regions) ./ eff_pop)' .* RR).* total_inf'));  
     
    }
  }
}
  

model {
 phi1 ~ normal(0,5);
 init_R ~ normal(-0.1, 0.5); // Approximately Normal(1, 0.5)
 for (m in 1:M_regions){
   for (i in 1:(N-1)){
     rw_noise[i,m] ~ normal(0,1);
   }
 }

 rw_sd ~ normal(0, 0.05) T[0,];
 // for (m in 1:M_regions){
 //   for (i in 1:initial_seeding_day){
 //      Rt[i,m] ~ normal(2,1) ;
 //   }
 // for (i in initial_seeding_day:final_time){
 //   Rt[i,m] ~ normal(Rt[i-1,m],0.3) ;
 //   }
 // }
   for (t in fitting_start:final_time){
     for (m in 1:M_regions){
       // data_inf[t, m] ~ poisson(infection[t, m]);
     target +=  neg_binomial_2_lpmf(data_inf[t, m]| infection[t, m], phi1);
    }
  }
}

generated quantities{
  array[prediction_horizon, M_regions] real<lower=0> infection_forecast;     // daily initialization
  if (prediction_horizon >0){
    for (tt in 1:prediction_horizon){
      for (mm in 1:M_regions){
        infection_forecast[tt,mm] = poisson_rng(infection[final_time+tt, mm]);
      }
    }
  }
}
