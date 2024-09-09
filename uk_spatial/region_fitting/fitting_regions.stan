data {
  int<lower=1> M_regions;
  int<lower=1> final_time;    //number of days to simulate
  int W;      // number weeks to simulate
  array[final_time, 5, M_regions] real gmobility;
  int<lower=1> initial_seeding_day;
  int<lower=1> death_data_length;
  array[death_data_length, M_regions] int<lower=0> death;       //number of infected individual at any time
  vector[final_time]  SI;
  vector[final_time] f;
  vector[M_regions] pop;
  matrix[M_regions,M_regions] C_base;
  matrix[M_regions,M_regions] C_lockdown;
  array[final_time] int<lower = 0> day_week_index;
  array[final_time, 4] int I;
  int first_lockdown_end;
  int second_lockdown_end;
  int fitting_death_start;
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
  array[M_regions] real x1;
  array[M_regions] real x2;
  array[M_regions] real x3;
  array[M_regions] real x4;
  array[M_regions] real x5;
  array[M_regions] real y1;
  array[M_regions] real y2;
  array[M_regions] real y3;
  array[M_regions] real y4;
  array[M_regions] real y5;
  array[M_regions] real z1;
  array[M_regions] real z2;
  array[M_regions] real z3;
  array[M_regions] real z4;
  array[M_regions] real z5;
}

transformed parameters{
   matrix[final_time, M_regions] infection = rep_matrix(0, final_time, M_regions);    // daily initialization
   matrix[final_time, M_regions] daily_deaths = rep_matrix(0,final_time, M_regions);      // daily deaths (infection becomes a case after latent period)
   matrix[death_data_length, M_regions] weekly_deaths = rep_matrix(0, death_data_length, M_regions);
   matrix[W,M_regions] weekly_effect = rep_matrix(0,W,M_regions);   // check why this +1 is needed
   matrix[ final_time, M_regions] Rt = rep_matrix(0, final_time, M_regions);   //reproduction number
    
   //-------------------------------------------------------------------------------------------------//
   matrix[final_time, M_regions] SI_regions = rep_matrix(SI_rev,M_regions);  // serial interval for every region
   matrix[final_time, M_regions] f_regions = rep_matrix(f_rev,M_regions);    // infection to death distribution for every region
    
   for (m in 1:M_regions){                  // for initial seeding
      infection[1:initial_seeding_day,m] = rep_vector(initial_seeding[m], initial_seeding_day);  //initial_seeding[m]
      weekly_effect[:,m] = weekly_var * cumulative_sum(weekly_effect_d[:,m]);
      Rt[:,m] = mu[m] * 2 * inv_logit(- weekly_effect[day_week_index, m] 
                                      - x1[m] * (to_vector(gmobility[:,1,m]) .* to_vector(I[:,2])) - x2[m] * (to_vector(gmobility[:,2,m]) .* to_vector(I[:,2])) - x3[m] * (to_vector(gmobility[:,3,m]) .* to_vector(I[:,2])) - x4[m] * (to_vector(gmobility[:,4,m]) .* to_vector(I[:,2])) - x5[m] * (to_vector(gmobility[:,5,m]) .* to_vector(I[:,2]))
                                      - y1[m] * (to_vector(gmobility[:,1,m]) .* to_vector(I[:,3])) - y2[m] * (to_vector(gmobility[:,2,m]) .* to_vector(I[:,3])) - y3[m] * (to_vector(gmobility[:,3,m]) .* to_vector(I[:,3])) - y4[m] * (to_vector(gmobility[:,4,m]) .* to_vector(I[:,3])) - y5[m] * (to_vector(gmobility[:,5,m]) .* to_vector(I[:,3]))
                                      - z1[m] * (to_vector(gmobility[:,1,m]) .* to_vector(I[:,4])) - z2[m] * (to_vector(gmobility[:,2,m]) .* to_vector(I[:,4])) - z3[m] * (to_vector(gmobility[:,3,m]) .* to_vector(I[:,4])) - z4[m] * (to_vector(gmobility[:,4,m]) .* to_vector(I[:,4])) - z5[m] * (to_vector(gmobility[:,5,m]) .* to_vector(I[:,4])));
    }
    
   for (t in (initial_seeding_day + 1) : final_time){ //for loop over time
      row_vector[M_regions] convolution_inf = columns_dot_product(infection[1:(t-1),:] , SI_regions[(final_time-t+2):final_time,:]);  //infections at each region "k"
      matrix[M_regions,M_regions] C;
      if (I[t,1] == 0){
        C = C_base ;
      }else{
        C = C_lockdown ;
      }

      // C = C_base;
      
      vector[M_regions] total_inf = C * convolution_inf';     //total infection at region "j"
      vector[M_regions] eff_pop = C * pop;
      
      for (m in 1:M_regions){
        real sus = pop[m] - sum(infection[1:(t-1),m]);
        infection[t,m] = dot_product(C[:,m]', (((rep_vector(sus , M_regions) ./ eff_pop)'.* Rt[t,:]).* total_inf'));
      }
    }
   for (m in 1:M_regions){
      daily_deaths[1,m] = 1e-15 * infection[1,m];
      for (t in 2:final_time){
        daily_deaths[t,m] = ifr_noise[m] * dot_product(sub_col(infection, 1 , m, t-1), tail(f_rev, t-1));  // ifr_noise[m] * 
      }
    }
  
   for (m in 1:M_regions){
     for (t in 1:(final_time %/% 7)) { //weekly_deaths
        weekly_deaths[t,m] = sum(daily_deaths[(7*(t-1)+1):7*t,m]);
    }
    weekly_deaths[(final_time %/% 7)+1,m] = sum(daily_deaths[(7*(final_time %/% 7)+1) : final_time,m]);
  }
}

model {
 phi1 ~ normal(0,5);
 tau ~ exponential(0.03);
 weekly_var ~ normal(0,.2);
 // kappa ~ normal(0,0.5);
 mu ~ normal(3.28,0.5);
 initial_seeding ~ exponential(1/tau);
 ifr_noise ~ normal(1,0.1);
 gamma ~ normal(0,0.5);
 x1 ~ normal(0,gamma);
 x2 ~ normal(0,gamma);
 x3 ~ normal(0,gamma);
 x4 ~ normal(0,gamma);
 x5 ~ normal(0,gamma);
 y1 ~ normal(0,gamma);
 y2 ~ normal(0,gamma);
 y3 ~ normal(0,gamma);
 y4 ~ normal(0,gamma);
 y5 ~ normal(0,gamma);
 z1 ~ normal(0,gamma);
 z2 ~ normal(0,gamma);
 z3 ~ normal(0,gamma);
 z4 ~ normal(0,gamma);
 z5 ~ normal(0,gamma);
  
  for (m in 1:M_regions){
    weekly_effect_d[1,m] ~ normal(0, 1);
    for (week in 2:W){
      weekly_effect_d[week,m] ~ normal(0, 1);
    }
  }
  
 for (week in fitting_death_start:death_data_length){
    target += neg_binomial_2_lpmf(death[week, :]| weekly_deaths[week, :], phi1);
 }
}
