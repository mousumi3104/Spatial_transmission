functions {
  array[] real geometric_random_walk( int N, real init_R, array[] real rw_noise, real rw_sd) {
  array[N] real Rt; 
  Rt[1] = init_R; // Initialize with the given initial value for each region
  for (i in 2: N) {
    Rt[i] = Rt[i - 1] + rw_noise[i - 1] * rw_sd; // Add noise multiplied by rw_sd
  }
  return exp(Rt); // Return the exponential of the resulting array
  }
}

data {
  int<lower=1> final_time;    //number of days to simulate
  int N;
  int<lower=1> initial_seeding_day;
  int data_length;
  array[N] int<lower=0> data_fit;       //number of deaths
  vector[N]  SI;
  int pop;
  int fitting_start;
  int prediction_horizon; 
}

transformed data {
  vector[N] SI_rev; // SI in reverse order

  for(i in 1:N){
    SI_rev[i] = SI[N-i+1];
  }
}

parameters {
  real<lower=0> phi2;
  real init_R;
  array[N-1] real rw_noise; // random walk noise
  real<lower = 0> rw_sd; // random walk standard deviation
}

transformed parameters{
  vector[N] infection = rep_vector(0, N);    // daily initialization
  
  infection[1:initial_seeding_day] = rep_vector(30, initial_seeding_day);      // learn the number of cases in the first initial seeding days
  array[N] real<lower=0> Rt = geometric_random_walk(N,init_R, rw_noise, rw_sd);
  
  for (t in (initial_seeding_day+1):N){ //for loop over time
      infection[t] = ((pop - sum(infection[1:(t-1)]))/pop) * Rt[t] * dot_product(infection[1: (t-1)], tail(SI_rev, t-1));
    }
}

model {
 phi2 ~ normal(0,5);
 init_R ~ normal(-0.1, 0.5); // Approximately Normal(1, 0.5)
 
   for (i in 1:(N-1)){
     rw_noise[i] ~ normal(0,1);
   }
 

 rw_sd ~ normal(0, 0.05) T[0,];
 // Rt[1] ~ normal(3.28,1);
 // for (i in 2:final_time){
 //   Rt[i] ~ normal(Rt[i-1],0.5) ;
 // }
 
 for (t in fitting_start:final_time){
  //  target += poisson_lpmf(data_fit[t]| deaths[t]);
    target += neg_binomial_2_lpmf(data_fit[t]| infection[t],phi2);
 }
}

generated quantities{
  vector[prediction_horizon] infection_forecast = rep_vector(0, prediction_horizon);    // daily initialization
  if (prediction_horizon >0){
    for (tt in 1:prediction_horizon){
      infection_forecast[tt] = poisson_rng(infection[final_time+tt]);
    }
  }
}
