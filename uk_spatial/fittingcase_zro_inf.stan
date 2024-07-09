data {
  int<lower=1> M_regions;
  int<lower=1> final_time;    //number of days to simulate
  int<lower=1> initial_seeding_day;
  int<lower=1> initial_seeding[M_regions];
  int<lower=0> death[final_time/7,M_regions];       //number of infected individual at any time
  real SI[final_time];
  real f[final_time];
  vector[M_regions] pop;
  matrix[M_regions,M_regions] C;
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
  matrix<lower=0.1>[final_time,M_regions] rt;    //time varying reproduction number
  real<lower=0> phi2;
  real<lower=0> ifr; 
}

transformed parameters{
  
  matrix[final_time, M_regions] infection = rep_matrix(0, final_time, M_regions);    // daily initialization
  matrix[final_time, M_regions] daily_deaths = rep_matrix(0,final_time, M_regions);      // daily deaths (infection becomes a case after latent period)
  matrix[final_time/7, M_regions] weekly_deaths = rep_matrix(0, final_time/7, M_regions);
  // matrix[final_time, M_regions] final_infection ;     // total infection till date
  
  for (i in 1:M_regions){                  // for initial seeding
    infection[1:initial_seeding_day,i] = rep_vector(initial_seeding[i],initial_seeding_day);      // learn the number of cases in the first initial seeding days
    daily_deaths[1:initial_seeding_day,i] = 1e-15 * infection[1:initial_seeding_day,i];
  }
  matrix[final_time,M_regions] SI_regions = rep_matrix(SI_rev,M_regions); 
  
  for (t in (initial_seeding_day+1):final_time){ //for loop over time
    row_vector[M_regions] convolution_inf = columns_dot_product(infection[1:(t-1),:] , SI_regions[(final_time-t+2):final_time,:]);  //infections at each region "k"
    vector[M_regions] total_inf = C * convolution_inf';     //total infection at region "j"
    vector[M_regions] eff_pop = C * pop;
      
    for (i in 1:M_regions){      // for loop over region "i" (final infection at region "i")   
      real sus = pop[i]-sum(infection[1:(t-1),i]);
      infection[t,i] = dot_product(C[:,i]', (((rep_vector(sus , M_regions) ./ eff_pop)' .* rt[t,:]).* total_inf'));
      daily_deaths[t,i] = ifr * dot_product(infection[1:(t-1),i], tail(f_rev,t-1));       
      
    }
  }
  for (t in 1:(final_time/7)) { //weekly_deaths
    for (i in 1:M_regions){
        weekly_deaths[t,i] = sum(daily_deaths[(7*(t-1)+1):7*t,i]);
    }
  }
}

model {
 phi2 ~ normal(0,5);
 ifr ~ normal(0.5,0.1);
 for (t in 1:initial_seeding_day){
     rt[t,1] ~ normal(1,0.01);
     rt[t,2] ~ normal(1,0.01);
 }
 for (t in (initial_seeding_day + 1):final_time){
     rt[t,1] ~ normal(rt[t-1,1],0.02);
     rt[t,2] ~ normal(rt[t-1,2],0.02);

     if((t % 7) == 0){
       target += neg_binomial_2_lpmf(death[t/7,:]|weekly_deaths[t/7,:],phi2);
   }
  }
}




