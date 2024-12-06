// Reproduction number(Rt) flollows from the region where the population are having contact with the people
data {
  int<lower = 1> M;        //number of regions
  vector[M] pop;                 //population of three regions
  int<lower = 1> final_time;        // total time point to run the simulation
  int<lower = 0> initial_seeding_day;
  array[M] int init_seed;       //initial seeding time    
  vector[final_time] SI;      //serial interval
  matrix[final_time,M] Rt;
  matrix[M, M] C;                // mobility matrix
  //array[final_time] int<lower=0> day_week_index;
  // int week;
}

transformed data {
  vector[final_time] SI_rev = reverse(SI);      //reverse of the SI (required for simulation)
  // vector[final_time] f1_rev = reverse(f1);       //reverse of the infection-death distribution
  // vector[final_time] f2_rev = reverse(f2);       //reverse of the infection-onset distribution
}

parameters{

}

generated quantities{
  matrix[final_time, M] infection = rep_matrix(0, final_time, M);    // daily initialization
  matrix[final_time, M] infection_in_own = rep_matrix(0, final_time, M);
  matrix[final_time, M] infection_in_mob = rep_matrix(0, final_time, M);
  matrix[final_time, M] infection_out_mob = rep_matrix(0, final_time, M);
  matrix[final_time, M] cumm_sum = rep_matrix(0, final_time, M);     // total infection to date
  matrix[final_time, M] daily_cases = rep_matrix(0,final_time, M);
  matrix[final_time, M] daily_deaths = rep_matrix(0,final_time, M);      // daily cases (infection becomes a case after latent period)
  matrix[final_time %/% 7, M] weekly_deaths = rep_matrix(0, final_time %/% 7, M); 
  
  
  
  for (i in 1:M){                  // for initial seeding
    infection[1:initial_seeding_day,i] = rep_vector(init_seed[i],initial_seeding_day);      // learn the number of cases in the first initial seeding days
    cumm_sum[2:initial_seeding_day,i] = cumulative_sum(infection[2:initial_seeding_day,i]);    // cumulative infection
    daily_cases[1:initial_seeding_day,i] = 1e-15 * infection[1:initial_seeding_day,i];
    daily_deaths[1:initial_seeding_day,i] = 1e-15 * infection[1:initial_seeding_day,i];
  }
  matrix[final_time,M] SI_regions = rep_matrix(SI_rev,M); 
  
  ///////////////
  for (t in (initial_seeding_day+1):final_time){ //for loop over time
    vector[M] convolution_inf = (columns_dot_product(infection[1:(t-1),:] , SI_regions[(final_time-t+2):final_time,:]))';  //infections at each region "k"
    vector[M] total_inf = C * convolution_inf;     //total infection at region "j"
    vector[M] eff_pop = C * pop;
  
    for (i in 1:M){      // for loop over region "i" (final infection at region "i")  
      real iar_noise = 1;//normal_rng(1,0.2);
      real ifr_noise = 1;//normal_rng(1,0.1);
      real sus = pop[i]-sum(infection[1:(t-1),i]);
  
      infection[t,i] = poisson_rng(dot_product(C[:,i]', (((rep_vector(sus , M) ./ eff_pop)' .* Rt[t,:]).* total_inf')));
      
      infection_in_own[t,i] = (sus * C[i,i] / eff_pop[i]) * Rt[t,i] * (C[i,i] * convolution_inf[i]) ;
      infection_in_mob[t,i] = (sus * C[i,i] / eff_pop[i]) * Rt[t,i] * (total_inf[i] - (C[i,i] * convolution_inf[i])) ;
      array[M-1] int index; 
      int ind = 1;
      for (n  in 1:M){
        if (n != i){
          index[ind] = n;
          ind += 1;
        }
      }
      infection_out_mob[t,i] = dot_product(C[index,i]', (((rep_vector(sus , M-1) ./ eff_pop[index])' .* Rt[t,index]).* total_inf[index]'));
      
      // daily_cases[t,i] = iar_noise * dot_product(infection[1:(t-1),i], tail(f2_rev,t-1));
      // daily_deaths[t,i] =  ifr_noise * dot_product(infection[1:(t-1),i], tail(f1_rev,t-1)); 
      // cumm_sum[t,i] = cumm_sum[t-1,i] + infection[t-1,i];
    }
  }
  // for (t in 1:(final_time %/% 7)) {
  //   for (i in 1:M){
  //     weekly_deaths[t,i] = neg_binomial_2_rng(sum(daily_deaths[(7*(t-1)+1):(7*t),i]),1);
  //   }
  // }
}

