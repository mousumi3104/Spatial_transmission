// functions {
//   // Partial sum function for reduce_sum
//   real convolution(int slice_week, int start, int end, int M_regions, matrix death, matrix weekly_deaths, real phi1) {
//     real target_acc = 0;
//     for (week in start:end) {
//       target_acc += neg_binomial_2_lpmf(death[week,:] | weekly_deaths[week,:], phi1);
//     }
//     return target_acc;
//   }
// }

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
  matrix[final_time, M_regions] gm = rep_matrix(0, final_time, M_regions);
  for (m in 1:M_regions){
    gm[:,m] = (to_vector(gmobility[:,1,m]) + to_vector(gmobility[:,2,m]) + to_vector(gmobility[:,3,m]) + to_vector(gmobility[:,4,m]) + to_vector(gmobility[:,5,m]))/5;
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
  array[M_regions] real x;
  array[M_regions] real y;
  array[M_regions] real z;
}

transformed parameters{
   matrix[final_time, M_regions] infection = rep_matrix(0, final_time, M_regions);    // daily initialization
   matrix[final_time, M_regions] daily_deaths = rep_matrix(0,final_time, M_regions);      // daily deaths (infection becomes a case after latent period)
   matrix[death_data_length, M_regions] weekly_deaths = rep_matrix(0, death_data_length, M_regions);
   matrix[W,M_regions] weekly_effect = rep_matrix(0,W,M_regions);   // check why this +1 is needed
   matrix[final_time, M_regions] Rt = rep_matrix(0, final_time, M_regions);   //reproduction number
   
    for (m in 1:M_regions){
      infection[1:initial_seeding_day,m] = rep_vector(initial_seeding[m],initial_seeding_day);
      weekly_effect[:,m] = weekly_var * cumulative_sum(weekly_effect_d[:,m]) ;
      Rt[1:initial_seeding_day,m] = rep_vector(mu[m] * 2 * inv_logit(weekly_effect[1,m]) , initial_seeding_day) ; 
    }
   
   for (t in (initial_seeding_day+1):final_time){
     
     matrix[M_regions,M_regions] C ;
      if (I[t,1] == 0){
        C = C_base ;
      }else{
        C = C_lockdown ;
      }
      
     vector[M_regions] conv = rep_vector(0,M_regions);
      for (m in 1:M_regions){
       Rt[t,m] = mu[m] * 2 * inv_logit(weekly_effect[day_week_index[t],m] - x[m] * gm[t,m] * I[t,2] - y[m] * gm[t,m] * I[t,3] - z[m] * gm[t,m] * I[t,4]);
       conv[m] = dot_product(infection[1:(t-1),m], tail(SI_rev, (t-1)));
       daily_deaths[t,m] = ifr_noise[m] * dot_product(infection[1:(t-1),m], tail(f_rev, (t-1)));
      }
      for (m in 1: M_regions){
        for (j in 1:M_regions){
          real sus_jm = C[j,m] * (1 - (sum(infection[1:(t-1),m])/pop[m]));
          real total_inf_j = 0;
          for (k in 1:M_regions){
            total_inf_j = total_inf_j + C[j,k] * conv[k];
          }
          real total_pop_j = 0;
          for (l in 1:M_regions){
            total_pop_j = total_pop_j + C[j,l] * pop[l];
          }
          infection[t,m] = infection[t,m] + (sus_jm / total_pop_j) * Rt[t,j] * total_inf_j ;
        }
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
 tau ~ exponential(0.01);
 weekly_var ~ normal(0,.2);
 mu ~ normal(3.28,0.5);
 initial_seeding ~ exponential(1/tau);
 ifr_noise ~ normal(1,0.1);
 gamma ~ normal(0,0.5);

  for (m in 1:M_regions){
    x[m] ~ normal(0,gamma);
    y[m] ~ normal(x[m],gamma);
    z[m] ~ normal(y[m],gamma);
    weekly_effect_d[1,m] ~ normal(0, 1);
    for (week in 2:W){
      weekly_effect_d[week,m] ~ normal(0, 1);
    }
  }
  
 for (week in fitting_death_start:death_data_length){
    target += neg_binomial_2_lpmf(death[week, :]| weekly_deaths[week, :], phi1);
 }
}









