data {
  int<lower=1> M_regions;
  int<lower=1> final_time;    //number of days to simulate
  int<lower=1> W;
  matrix[final_time, M_regions] gmobility;
  int<lower=1> initial_seeding_day;
  int<lower=1> death_data_length;
  array[death_data_length, M_regions] int<lower=0> death;       //number of infected individual at any time
  vector[final_time]  SI;
  vector[final_time] f;
  vector[M_regions] pop;
  matrix[M_regions,M_regions] C_base;
  matrix[M_regions,M_regions] C_lockdown;
  array[final_time] int day_week_index;
  matrix[final_time, 4] I;
  int first_lockdown_end;
  int second_lockdown_end;
  int infection_gen_time;
}

transformed data {
  vector[final_time] SI_rev; // SI in reverse order
  vector[final_time] f_rev; // f in reversed order
  
  for(t in 1:final_time){
    SI_rev[t] = SI[final_time-t+1];
  }
  for(t in 1:final_time) {
     f_rev[t] = f[final_time-t+1];       // infection to death
    }
}

parameters {
  array[M_regions] real<lower=0> mu;      //parameters for Rt
  // real<lower = 0> kappa;
  array[M_regions] real<lower=1> initial_seeding;
  real<lower=0> tau;
  real<lower=0> weekly_var;       //weekly variance
  matrix[W, M_regions] weekly_effect_d;     //parameters for Rt (why W+1).     ?????
  real<lower=0> phi1;
  array[M_regions] real<lower=0> ifr_noise;
  real<lower = 0> gamma;
  vector[M_regions] x1;
  vector[M_regions] y1;
  vector[M_regions] z1;
}

transformed parameters{
   
   matrix[final_time, M_regions] infection;    // daily initialization
   matrix[final_time, M_regions] daily_deaths;      // daily deaths (infection becomes a case after latent period)
   matrix[W,M_regions] weekly_effect;   // check why this +1 is needed
   matrix[ final_time, M_regions] Rt;   //reproduction number
    
   //-------------------------------------------------------------------------------------------------//
   
   matrix[final_time, M_regions] SI_regions = rep_matrix(SI_rev,M_regions);  // serial interval for every region
   matrix[final_time, M_regions] f_regions = rep_matrix(f_rev,M_regions);    // infection to death distribution for every region
   
   for (m in 1:M_regions){                  // for initial seeding
      infection[1:initial_seeding_day,m] = rep_vector(initial_seeding[m], initial_seeding_day);  //initial_seeding[m]
      weekly_effect[:,m] = weekly_var * cumulative_sum(weekly_effect_d[:,m]);
      Rt[1:initial_seeding_day, m] = rep_vector(mu[m] * 2 * inv_logit(- weekly_effect[1, m]) , initial_seeding_day);
    }
   for (m in 1:M_regions){
      for (t in (initial_seeding_day+1):final_time){
          Rt[t,m] = mu[m] * 2 * inv_logit(- weekly_effect[day_week_index[t], m]
                                     - (x1[m] * gmobility[t,m] * I[t,2])
                                     - (y1[m] * gmobility[t,m] * I[t,3])
                                     - (z1[m] * gmobility[t,m] * I[t,4]));
      }
     }
   
   for (t in (initial_seeding_day + 1) : final_time){ //for loop over time
      vector[M_regions] convolution_inf = (columns_dot_product(infection[1:(t-1),:] , SI_regions[(final_time-t+2):final_time,:]))';  //infections at each region "k"
      matrix[M_regions,M_regions] C;
      
      if (I[t,1] == 0){
        C = C_base ;
      }else{
        C = C_lockdown ;
      }
      
      vector[M_regions] total_inf = C * convolution_inf;     //total infection at region "j"
      vector[M_regions] eff_pop = C * pop;
      
      for (m in 1:M_regions){
        real sus = pop[m] - sum(infection[1:(t-1),m]);
        infection[t,m] = dot_product(C[:,m]', (((rep_vector(sus , M_regions) ./ eff_pop)'.* Rt[t,:]).* total_inf'));
      }
    }
    daily_deaths[1,:] = 1e-15 * infection[1,:];
    for (t in 2:final_time){
      daily_deaths[t,:] = (to_vector(ifr_noise))' .* columns_dot_product(infection[1:(t-1),:] , f_regions[(final_time-t+2):final_time,:]);
    }
 }

model {
 phi1 ~ normal(0,5);
 tau ~ exponential(0.01);
 weekly_var ~ normal(0,.2);
 // kappa ~ normal(0,0.5);
 mu ~ normal(3.28,0.5);
 initial_seeding ~ exponential(1/tau);
 to_vector(ifr_noise) ~ normal(1,0.1);
 gamma ~ normal(0,0.5);
 
 x1 ~ normal(0, gamma);
 y1 ~ normal(0, gamma);
 z1 ~ normal(0, gamma);
 to_vector(weekly_effect_d) ~ normal(0, 1);
 
  for (t in 1:death_data_length){
    for (m in 1:M_regions){
      target += neg_binomial_2_lpmf(death[t,m]| daily_deaths[t + infection_gen_time ,m], phi1);
    }
  }
}
