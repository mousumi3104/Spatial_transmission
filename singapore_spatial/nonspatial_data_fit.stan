 data {
  int<lower=0> final_time;    //number of days to simulate
  int<lower=0> initial_seeding_day;
  vector[initial_seeding_day] initial_seeding;
  int incidence[final_time];       //number of infected individual at any time
  real SI[final_time];
  int pop;
}

transformed data {
  vector[final_time] SI_rev; // SI in reverse order
  for(i in 1:final_time)
    SI_rev[i] = SI[final_time-i+1];
}

parameters {
  real<lower=0> Rt[final_time];    //time varying reproduction number
  real<lower=0> phi2;     // the 
 // real<lower=0> y;
 // real<lower=0> tau;
}

transformed parameters{
  vector[final_time] prediction = rep_vector(0,final_time);
  vector[final_time] Rt_adj = rep_vector(0,final_time);
  vector[final_time] cumm_sum =rep_vector(0, final_time);
  
  prediction[1:initial_seeding_day] = initial_seeding;
  for (t in 2:initial_seeding_day){
    cumm_sum[t] = cumm_sum[t-1]+prediction[t-1];
    Rt_adj[t] = Rt[t];    // have to check the dimension
  }
  
  for (t in (initial_seeding_day+1):final_time){
    real convolution = dot_product(prediction[1:t-1] , tail(SI_rev,(t-1)));
    cumm_sum[t] = cumm_sum[t-1] + prediction[t-1];
    Rt_adj[t] = ((pop-cumm_sum[t])/pop)*Rt[t];
    prediction[t] =  Rt_adj[t]* convolution;   
  }
}

model {
  //tau ~ exponential(0.03);
  //y ~ exponential(1/tau);
  
  phi2 ~ normal(0,5);
  //iar ~ normal(1,0.1);
  Rt[1] ~ normal(2,0.1);
  for (t in 2:final_time){
    Rt[t] ~ normal(Rt[t-1],0.2);
    incidence[t] ~ neg_binomial_2(prediction[t] , phi2);
  }
}
