data {
  int<lower=1> M_regions;
  matrix[M_regions,M_regions] C;
  int<lower=1> final_time;    //number of days to simulate
  int<lower=1> initial_seeding_day;
  int<lower=1> data_length;
  // array[data_length,M_regions] int<lower =0> data_deaths;       //number of infected individual at any time
  array[data_length,M_regions] int<lower =0> data_inf;       //number of infected individual at any time
  // array[data_length,M_regions] int<lower =0> data_cases;
  // matrix[final_time, M_regions+1] cases;
  vector[final_time]  SI;
  vector[final_time] f1;
  vector[final_time] f2;
  vector[M_regions] pop;
  int fitting_start;
}

transformed data {
  vector[final_time] SI_rev; // SI in reverse order
  vector[final_time] f1_rev; // f in reversed order
  vector[final_time] f2_rev; 
  for(t in 1:final_time){
    SI_rev[t] = SI[final_time-t+1];
  }
  for(t in 1:final_time) {
     f1_rev[t] = f1[final_time-t+1];       // infection to death
     
    }
    for(t in 1:final_time) {
     f2_rev[t] = f2[final_time-t+1];       // infection to death
    }
}

parameters {
  // array[M_regions] real<lower=0> mu;      //parameters for Rt
  // array[M_regions] real<lower=0> initial_seeding;
  // real<lower=0> tau;
  real<lower=0> phi1;
  // array[M_regions] real<lower=0>ifr_noise;
  // array[M_regions] real<lower=0>iar_noise;
  // real<lower=0> weekly_var;
  array[final_time,M_regions] real<lower=0> Rt;
  
}

transformed parameters{
  matrix[final_time, M_regions] infection = rep_matrix(0, final_time, M_regions);    // daily initialization
  matrix[final_time, M_regions] daily_deaths = rep_matrix(0,final_time, M_regions);      // daily deaths (infection becomes a case after latent period)
  // matrix[final_time, M_regions] daily_cases = rep_matrix(0,final_time, M_regions);      // daily deaths (infection becomes a case after latent period)

  // matrix[final_time/7, M_regions] weekly_deaths = rep_matrix(0, final_time/7, M_regions);
  // matrix[W,M_regions] weekly_effect = rep_matrix(0,W,M_regions);
  // matrix[final_time,M_regions] Rt = rep_matrix(0,final_time,M_regions);

  //-------------------------------------------------------------------------------------------------//
  {
    
  matrix[final_time, M_regions] SI_regions = rep_matrix(SI_rev,M_regions);  // serial interval for every region
  // matrix[final_time, M_regions] f_regions = rep_matrix(f1_rev,M_regions);    // infection to death distribution for every region
  
    for (m in 1:M_regions){                  // for initial seeding
      infection[1:initial_seeding_day,m] = rep_vector(30,initial_seeding_day);
      // daily_deaths[1:initial_seeding_day,m] = 1e-15 * infection[1:initial_seeding_day,m];      // daily deaths (infection becomes a case after latent period)

      
      // weekly_effect[:, m] = weekly_var * cumulative_sum(weekly_effect_d[:, m]) ;    // weekly effect
      // Rt[1:initial_seeding_day,m] =  rep_vector(mu[m] .* 2 * inv_logit(- weekly_effect[1,m]),initial_seeding_day);   // why this is weekly_effect[m]??
    }
    
  for (t in (initial_seeding_day+1):final_time){ //for loop over time
  
    row_vector[M_regions] convolution_inf = columns_dot_product(infection[1:(t-1),:] , SI_regions[(final_time-t+2):final_time,:]);  //infections at each region "k"
    vector[M_regions] total_inf = C * convolution_inf';     //total infection at region "j"
    vector[M_regions] eff_pop = C * pop;
    
    for (m in 1:M_regions){      // for loop over region "i" (final infection at region "i")   

      real sus = pop[m] - sum(infection[1:(t-1),m]);
      
      // Rt[t,m] = mu[m] * 2 * inv_logit(- weekly_effect[day_week_index[t],m]);   // why this is weekly_effect[m]??
      row_vector[M_regions] RR ;
      for (mm in 1:M_regions){
        RR[mm] = Rt[t,mm];
      }
      infection[t,m] = dot_product(C[:,m]', (((rep_vector(sus , M_regions) ./ eff_pop)' .* RR).* total_inf'));  
      // daily_cases[t,m] = iar_noise[m]*dot_product(sub_col(infection, 1 , m, t-1), tail(f2_rev, t-1));
      // daily_deaths[t,m] =  dot_product(sub_col(infection, 1 , m, t-1), tail(f1_rev, t-1));
    }
  }
  
  // for (t in 2:final_time){
  //   for (m in 1:M_regions){
  //     daily_deaths[t,m] =  dot_product(sub_col(infection, 1 , m, t-1), tail(f_rev, t-1));
  //   }
  // }
  // for (m in 1:M_regions){
  //   for (t in 1:(final_time %/% 7)) { //weekly_deaths
  //       weekly_deaths[t,m] = sum(daily_deaths[(7*(t-1)+1):7*t,m]);
  //   }
  // }
  }
}

model {
 phi1 ~ normal(0,5);
 // tau ~ exponential(0.03);
 // initial_seeding ~ exponential(1/tau);

 // mu ~ normal(3.28,1);
 // ifr_noise ~ normal(1,0.1);
 // iar_noise ~ normal(1,0.1);
 // weekly_var ~ normal(0,0.2);
 // 
 // weekly_effect_d[1,M_regions] ~ normal(0,1);
 
 for (m in 1:M_regions){
   for (i in 1:initial_seeding_day){
      Rt[i,m] ~ normal(2,1) ;
   }
 for (i in initial_seeding_day:final_time){
   Rt[i,m] ~ normal(Rt[i-1,m],0.3) ;
   }
 }

 for (t in fitting_start:final_time){
   for (m in 1:M_regions){
    target +=  neg_binomial_2_lpmf(data_inf[t, m]| infection[t, m], phi1);
   }
 }
}
