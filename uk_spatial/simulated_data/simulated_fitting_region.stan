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
  array[data_length,M_regions] int<lower =0> data_inf;       //number of infected individual at any time
  vector[N] SI;
  vector[M_regions] pop;
  int fitting_start;
  // int prediction_horizon;
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
  // matrix[final_time, M_regions] infection_in_own = rep_matrix(0, N, M_regions);
  // matrix[final_time, M_regions] infection_in_mob = rep_matrix(0, N, M_regions);
  // matrix[final_time, M_regions] infection_out_mob = rep_matrix(0, N, M_regions);
  
   
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
      // infection_in_own[t,m] = (sus * C[m,m] / eff_pop[m]) * Rt[t,m] * (C[m,m] * convolution_inf[m]) ;
      // infection_in_mob[t,m] = (sus * C[m,m] / eff_pop[m]) * Rt[t,m] * (total_inf[m] - (C[m,m] * convolution_inf[m])) ;
      // array[M_regions-1] int index; 
      // int ind = 1;
      // for (n  in 1:M_regions){
      //   if (n != m){
      //     index[ind] = n;
      //     ind += 1;
      //   }
      // }
      // infection_out_mob[t,m] = dot_product(C[index,m]', (((rep_vector(sus , M_regions-1) ./ eff_pop[index])' .* Rt[t,index]).* total_inf[index]'));
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
   // for (t in fitting_start:final_time){
     for (m in 1:M_regions){
     data_inf[fitting_start:final_time, m] ~ neg_binomial_2( infection[fitting_start:final_time, m], phi1);
    }
  // }
}

// generated quantities{
//   array[prediction_horizon, M_regions] real<lower=0> infection_forecast;     // daily initialization
//   if (prediction_horizon >0){
//     for (tt in 1:prediction_horizon){
//       for (mm in 1:M_regions){
//         infection_forecast[tt,mm] = poisson_rng(infection[final_time+tt, mm]);
//       }
//     }
//   }
// }
generated quantities {
  array[final_time,M_regions] real log_lik;
  for (n in 1:final_time) {
    for (m in 1:M_regions){
      log_lik[n,m] = neg_binomial_2_lpmf(data_inf[n,m] | infection[n, m], phi1);
    }
  }
}
