data {
  int<lower=1> M_regions;
  int<lower=1> final_time;    //number of days to simulate
  int W;      // number weeks to simulate
  // array[final_time, 5, M_regions] real gmobility;
  int<lower=1> initial_seeding_day;
  int<lower=1> death_data_length;
  array[death_data_length,M_regions] int<lower =0> death;       //number of infected individual at any time
  vector[final_time]  SI;
  vector[final_time] f;
  vector[M_regions] pop;
  matrix[M_regions,M_regions] C_base;
  matrix[M_regions,M_regions] C_lockdown;
  array[final_time]int<lower = 0> day_week_index;
  array[final_time, 4] int<lower = 1> I;
  int fitting_death_start;
}

transformed data {
  vector[final_time] SI_rev; // SI in reverse order
  vector[final_time] f_rev; // f in reversed order
  vector[final_time] f_rev1;
  for(t in 1:final_time){
    SI_rev[t] = SI[final_time-t+1];
  }
  for(t in 1:final_time) {
     f_rev[t] = f[final_time-t+1];       // infection to death
    }
}

parameters {
  row_vector<lower=0>[M_regions] mu;      //parameters for Rt
  row_vector<lower=0>[M_regions] initial_seeding;
  real<lower=0> tau;
  real<lower=0> weekly_var;       //weekly variance
  matrix[W, M_regions] weekly_effect_d;     //parameters for Rt (why W+1).     ?????
  real<lower=0> phi1;
  row_vector<lower=0>[M_regions] ifr_noise;
  // matrix[final_time, M_regions] X11;
  // matrix[final_time, M_regions] X12;
  // matrix[final_time, M_regions] X13;
  // matrix[final_time, M_regions] X21;
  // matrix[final_time, M_regions] X22;
  // matrix[final_time, M_regions] X23;
  // matrix[final_time, M_regions] X31;
  // matrix[final_time, M_regions] X32;
  // matrix[final_time, M_regions] X33;
  // matrix[final_time, M_regions] X41;
  // matrix[final_time, M_regions] X42;
  // matrix[final_time, M_regions] X43;
  // matrix[final_time, M_regions] X51;
  // matrix[final_time, M_regions] X52;
  // matrix[final_time, M_regions] X53;
}

transformed parameters{
  matrix[final_time, M_regions] infection = rep_matrix(0, final_time, M_regions);    // daily initialization
  matrix[final_time, M_regions] daily_deaths = rep_matrix(0,final_time, M_regions);      // daily deaths (infection becomes a case after latent period)
  // matrix[final_time, M_regions] daily_cases = rep_matrix(0,final_time,M_regions);
  matrix[death_data_length, M_regions] weekly_deaths = rep_matrix(0, death_data_length, M_regions);
  matrix[W,M_regions] weekly_effect = rep_matrix(0,W,M_regions);   // check why this +1 is needed
  matrix[ final_time, M_regions] Rt = rep_matrix(0, final_time, M_regions);   //reproduction number

  //-------------------------------------------------------------------------------------------------//
  // {
  // 
  // matrix[final_time, M_regions] SI_regions = rep_matrix(SI_rev,M_regions);  // serial interval for every region
  // matrix[final_time, M_regions] f_regions = rep_matrix(f_rev,M_regions);    // infection to death distribution for every region
  // 
  //   for (m in 1:M_regions){                  // for initial seeding
  //     infection[1:initial_seeding_day,m] = rep_vector(initial_seeding[m],initial_seeding_day);
  // 
  //     weekly_effect[:, m] = weekly_var * cumulative_sum(weekly_effect_d[:, m]) ;    // weekly effect
  //     Rt[1:initial_seeding_day,m] = rep_vector(mu[m] * 2 * inv_logit(- weekly_effect[1,m]),initial_seeding_day);   // why this is weekly_effect[m]??
  //   }
  // 
  // for (t in (initial_seeding_day+1):final_time){ //for loop over time
  // 
  //   row_vector[M_regions] convolution_inf = columns_dot_product(infection[1:(t-1),:] , SI_regions[(final_time-t+2):final_time,:]);  //infections at each region "k"
  // 
  //   // matrix[M_regions,M_regions] C;
  // 
  //   // if (I[t,2] == 0){
  //   //    C = C_base;
  //   // }else{
  //   //    C = C_lockdown;
  //   // }
  // 
  //   vector[M_regions] total_inf = C_base * convolution_inf';     //total infection at region "j"
  //   vector[M_regions] eff_pop = C_base * pop;
  // 
  //   for (m in 1:M_regions){ 
  //     Rt[t,m] = mu[m] * 2;// for loop over region "i" (final infection at region "i")
  //     // Rt[t,m] = mu[m] * 2 * inv_logit(- X11[t,m] * gmobility[t,1,m] * I[t,2] - X12[t,m] * gmobility[t,1,m] * I[t,3] - X13[t,m] * gmobility[t,1,m] * I[t,4]
  //     //                                 - X21[t,m] * gmobility[t,2,m] * I[t,2] - X22[t,m] * gmobility[t,2,m] * I[t,3] - X23[t,m] * gmobility[t,2,m] * I[t,4]
  //     //                                 - X31[t,m] * gmobility[t,3,m] * I[t,2] - X32[t,m] * gmobility[t,3,m] * I[t,3] - X33[t,m] * gmobility[t,3,m] * I[t,4]
  //     //                                 - X41[t,m] * gmobility[t,4,m] * I[t,2] - X42[t,m] * gmobility[t,4,m] * I[t,3] - X43[t,m] * gmobility[t,4,m] * I[t,4]
  //     //                                 - X51[t,m] * gmobility[t,5,m] * I[t,2] - X52[t,m] * gmobility[t,5,m] * I[t,3] - X53[t,m] * gmobility[t,5,m] * I[t,4]
  //     //                                 - weekly_effect[day_week_index[t],m]);
  //     real sus = pop[m] - sum(infection[1:(t-1),m]);
  //     infection[t,m] = dot_product(C[:,m]', (((rep_vector(sus , M_regions) ./ eff_pop)'.* Rt[t,:]).* total_inf'));
  //    }
  // }
  // 
  // daily_deaths[1,:] = 1e-15 * infection[1,:];
  // for (t in 2:final_time){
  //   for (m in 1:M_regions){
  //     daily_deaths[t,m] = ifr_noise[m] * dot_product(sub_col(infection, 1 , m, t-1), tail(f_rev, t-1));
  //   }
  // }
  // for (m in 1:M_regions){
  //   for (t in 1:(final_time %/% 7)) { //weekly_deaths
  //       weekly_deaths[t,m] = sum(daily_deaths[(7*(t-1)+1):7*t,m]);
  //   }
  //   weekly_deaths[(final_time %/% 7)+1,m] = sum(daily_deaths[(7*(final_time %/% 7)+1) : final_time,m]);
  // }
  // }
}

model {
 phi1 ~ normal(0,5);
 tau ~ exponential(0.03);
 weekly_var ~ normal(0,.2);

    for (m in 1:M_regions){
      // for (t in 1:final_time){
      //     X11[t,m] ~ normal(0,1);
      //     X12[t,m] ~ normal(0,1);
      //     X13[t,m] ~ normal(0,1);
      //     X21[t,m] ~ normal(0,1);
      //     X22[t,m] ~ normal(0,1);
      //     X23[t,m] ~ normal(0,1);
      //     X31[t,m] ~ normal(0,1);
      //     X32[t,m] ~ normal(0,1);
      //     X33[t,m] ~ normal(0,1);
      //     X41[t,m] ~ normal(0,1);
      //     X42[t,m] ~ normal(0,1);
      //     X43[t,m] ~ normal(0,1);
      //     X51[t,m] ~ normal(0,1);
      //     X52[t,m] ~ normal(0,1);
      //     X53[t,m] ~ normal(0,1);
      // }
      mu[m] ~ normal(3.28,0.2);
      initial_seeding[m] ~ exponential(1/tau);
      ifr_noise[m] ~ normal(1,0.1);
      weekly_effect_d[1,m] ~ normal(0, 1);
      for (week in 2:W){
        weekly_effect_d[week,m] ~ normal(0, 1);
      }
    }

 for (week in fitting_death_start:death_data_length){
   for (m in 1:M_regions){
    target += neg_binomial_2_lpmf(death[week, m]| weekly_deaths[week, m], phi1);
   }
 }
}
