data {
  int<lower=1> M_regions;
  int<lower=1> final_time;    //number of days to simulate
  int W;      // number weeks to simulate
  int<lower=1> initial_seeding_day;
  //array[M_regions] int<lower=1> initial_seeding;
  array[W , M_regions] int<lower=0> death;       //number of infected individual at any time
  vector[final_time]  SI;
  vector[final_time] f;
  vector[M_regions] pop;
  matrix[M_regions,M_regions] C;
  array[final_time] int<lower=0> week_index;
  real<lower=0> ifr;
}

transformed data {
  vector[final_time] SI_rev; // SI in reverse order
  vector[final_time] f_rev; // f in reversed order
  for(i in 1:final_time){
    SI_rev[i] = SI[final_time-i+1];
  }
  for(i in 1:final_time) {
     f_rev[i] = f[final_time-i+1];
    }
}

parameters {
  array[M_regions] real<lower=0> mu;      //parameters for Rt
  array[M_regions] real <lower=0> initial_seeding;
  real<lower=0> tau;
  real<lower=0> kappa;
  real weekly_var;       //weekly variance
  matrix[W, M_regions] weekly_effect_d;     //parameters for Rt (why W+1).     ?????
  real<lower=0> phi2;
  //real<lower=0,upper =1> ifr;  
}

transformed parameters{
  matrix[final_time, M_regions] infection = rep_matrix(0, final_time, M_regions);    // daily initialization
  matrix[final_time, M_regions] daily_deaths = rep_matrix(0,final_time, M_regions);      // daily deaths (infection becomes a case after latent period)
  matrix[W, M_regions] weekly_deaths = rep_matrix(0, W, M_regions);
  matrix[W,M_regions] weekly_effect = rep_matrix(0,W,M_regions);   // check why this +1 is needed 
  matrix[ final_time, M_regions] Rt = rep_matrix(0, final_time, M_regions);   //reproduction number
  
  //////////////////////////////////////////
  {
  matrix[final_time, M_regions] SI_regions = rep_matrix(SI_rev,M_regions);  // serial interval for every region
  matrix[final_time, M_regions] f_regions = rep_matrix(f_rev,M_regions);    // infection to death distribution for every region
  
  
    for (m in 1:M_regions){                  // for initial seeding
      infection[1:initial_seeding_day,m] = rep_vector(initial_seeding[m],initial_seeding_day);      // learn the number of cases in the first initial seeding days
      weekly_effect[:, m] = weekly_var * weekly_effect_d[:, m] ;    // weekly effect
      Rt[:,m] = mu[m] * 2 * inv_logit(- weekly_effect[week_index,m]);   // why this is weekly_effect[m]??
    }
  
  for (t in (initial_seeding_day+1):final_time){ //for loop over time
  
    row_vector[M_regions] convolution_inf = columns_dot_product(infection[1:(t-1),:] , SI_regions[(final_time-t+2):final_time,:]);  //infections at each region "k"
    vector[M_regions] total_inf = C * convolution_inf';     //total infection at region "j"
    vector[M_regions] eff_pop = C * pop;
    
    for (m in 1:M_regions){      // for loop over region "i" (final infection at region "i")   
    
      real sus = pop[m]-sum(infection[1:(t-1),m]);
      infection[t,m] = dot_product(C[:,m]', (((rep_vector(sus , M_regions) ./ eff_pop)'.* Rt[t,:]).* total_inf'));
    }
  }
  
  daily_deaths[1,:] = 1e-15 * infection[1,:];
  for (t in 2:final_time){
    daily_deaths[t,:] = ifr * columns_dot_product(infection[1:(t-1),:], f_regions[(final_time-t+2):final_time,:]); 
  }
  for (t in 1:(final_time/7)) { //weekly_deaths
    for (m in 1:M_regions){
        weekly_deaths[t,m] = sum(daily_deaths[(7*(t-1)+1):7*t,m]);
    }
  }
  }
}

model {
 phi2 ~ normal(0,5);
 tau ~ exponential(0.03);
 
 for (m in 1:M_regions){
    initial_seeding[m] ~ exponential(1/tau);
 }
 //ifr ~ normal(0.0103,0.1);

 mu ~ normal(3.28,kappa);
 kappa ~ normal(0,0.5);
 
 weekly_var ~ normal(0,.2);
 to_vector(weekly_effect_d[1:W, ]) ~ normal(0, 1);      

 for (t in 1:W){
    target += neg_binomial_2_lpmf(death[t, ]| weekly_deaths[t, ], phi2);
 }
}
