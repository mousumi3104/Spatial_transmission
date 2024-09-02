data {
  int<lower=1> final_time;    //number of days to simulate
  int<lower=1> initial_seeding_day;
  int data_length;
  array[final_time] int<lower=0> data_fit;       //number of deaths
  vector[final_time]  SI;
  vector[final_time] f1;
  vector[final_time] f2;
  int pop;
  int fitting_start;
}

transformed data {
  vector[final_time] SI_rev; // SI in reverse order
  vector[final_time] f1_rev; // f in reversed order
  vector[final_time] f2_rev;
  
  for(i in 1:final_time){
    SI_rev[i] = SI[final_time-i+1];
  }
  for(i in 1:final_time) {
     f1_rev[i] = f1[final_time-i+1];
    }
    for(i in 1:final_time) {
     f2_rev[i] = f2[final_time-i+1];
    }
}

parameters {
  // real<lower=0> mu;      //parameters for Rt
  // real <lower=10> initial_seeding;
  // real<lower=0> tau;
  // real<lower=0> kappa;
  // real<lower=0> weekly_var;       //weekly variance
  // vector[W] weekly_effect_d;     //parameters for Rt (why W+1).     ?????
  real<lower=0> phi2;
  array[final_time] real<lower=0> Rt;
// real ifr_noise;
}

transformed parameters{
  vector[final_time] infection = rep_vector(0, final_time);    // daily initialization
  // vector[final_time] cases = rep_vector(0, final_time); 
  // vector[final_time] deaths = rep_vector(0,final_time);
  
  infection[1:initial_seeding_day] = rep_vector(30, initial_seeding_day);      // learn the number of cases in the first initial seeding days
  // deaths[1:initial_seeding_day] = 1e-15 * infection[1:initial_seeding_day];
  // cases[1:initial_seeding_day] = 1e-15 * infection[1:initial_seeding_day];
  
  for (t in (initial_seeding_day+1):final_time){ //for loop over time
      infection[t] = ((pop - sum(infection[1:(t-1)]))/pop) * Rt[t] * dot_product(infection[1: (t-1)], tail(SI_rev, t-1));
    }
    // deaths[1] = 1e-15 * infection[1];
    // cases[1] = 1e-15 * infection[1];
    // for (t in (initial_seeding_day+1):final_time){
    //   deaths[t] = dot_product(infection[1: (t-1)], tail(f1_rev, t-1));
    //   cases[t] =  dot_product(infection[1: (t-1)], tail(f2_rev, t-1));
    // }
}

model {
 phi2 ~ normal(0,5);
 // tau ~exponential(0.03);
 // initial_seeding ~ exponential(1/tau);
 // ifr_noise ~ normal(1,0.1);

 Rt[1] ~ normal(3.28,1);
 for (i in 2:final_time){
   Rt[i] ~ normal(Rt[i-1],0.5) ;
 }
 
 for (t in fitting_start:final_time){
  //  target += poisson_lpmf(data_fit[t]| deaths[t]);
    target += neg_binomial_2_lpmf(data_fit[t]| infection[t],phi2);
 }
}
