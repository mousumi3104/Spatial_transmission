data {
  int<lower=1> final_time;    //number of days to simulate
  int W;      // number weeks to simulate
  int<lower=1> initial_seeding_day;
  int data_length;
  array[data_length] int<lower=0> death;       //number of deaths
  vector[final_time]  SI;
  vector[final_time] f;
  int pop;
  array[final_time] int<lower=0> day_week_index;
  int fitting_start;
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
  real<lower=0> mu;      //parameters for Rt
  real <lower=0> initial_seeding;
  real<lower=0> tau;
  // real<lower=0> kappa;
  real<lower=0> weekly_var;       //weekly variance
  vector[W] weekly_effect_d;     //parameters for Rt (why W+1).     ?????
  real<lower=0> phi2;
}

transformed parameters{
  vector[final_time] infection = rep_vector(0, final_time);    // daily initialization
  vector[final_time] daily_deaths = rep_vector(0,final_time);      // daily deaths (infection becomes a case after latent period)
  vector[data_length] weekly_deaths = rep_vector(0, data_length);
  vector[W] weekly_effect = rep_vector(0,W);   // check why this +1 is needed
  vector[final_time] Rt = rep_vector(0, final_time);   //reproduction number

  //////////////////////////////////////////
  {

 // for initial seeding
      infection[1:initial_seeding_day] = rep_vector(initial_seeding, initial_seeding_day);      // learn the number of cases in the first initial seeding days
      weekly_effect = weekly_var * cumulative_sum(weekly_effect_d) ;    // weekly effect
      Rt[1:initial_seeding_day] = rep_vector(mu * 2 * inv_logit(- weekly_effect[1]),initial_seeding_day);   // why this is weekly_effect[m]??

  for (t in (initial_seeding_day+1):final_time){ //for loop over time
      // weekly_effect = weekly_var * cumulative_sum(weekly_effect_d) ;    // weekly effect
      Rt[t] = mu * 2 * inv_logit(- weekly_effect[day_week_index[t]]);   // why this is weekly_effect[m]??

      real convolution = dot_product(infection[1: (t-1)], tail(SI_rev, t-1));

      real sus = pop - sum(infection[1:(t-1)]);
      infection[t] = (sus/pop) * Rt[t] * convolution ;
    }

  daily_deaths[1] = 1e-15 * infection[1];
  for (t in 2:final_time){
    daily_deaths[t] =  dot_product(infection[1:(t-1)], tail(f_rev,t-1));
    }

  for (t in 1:(final_time/7)) {        //weekly_deaths
    weekly_deaths[t] = sum(daily_deaths[(7*(t-1)+1):7*t]);
    }
    weekly_deaths[(final_time/7)+1] = sum(daily_deaths[(7*(final_time/7)+1) : final_time]);
  }
}

model {
 phi2 ~ normal(0,5);
 tau ~ exponential(0.03);
 initial_seeding ~ exponential(1/tau);

 // kappa ~ normal(0,0.5);
 mu ~ normal(3.28,0.5);

 weekly_var ~ normal(0,.2);

 weekly_effect_d[1] ~ normal(0,1);

 for (i in 2:W){
   weekly_effect_d[i] ~ normal(0,1);
 }

 for (week in fitting_start:data_length){
    target += neg_binomial_2_lpmf(death[week]| weekly_deaths[week], phi2);
 }
}
