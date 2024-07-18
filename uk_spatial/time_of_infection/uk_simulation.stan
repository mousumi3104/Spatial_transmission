// Reproduction number(Rt) flollows from the region where the population are having contact with the people
data {
  int<lower = 1> M;        //number of regions
  vector[M] pop;                 //population of three regions
  int<lower = 1> final_time;        // total time point to run the simulation
  int<lower = 0> seed_time;       //initial seeding time
  real init_seed[M];              //initial seeding
  vector[final_time] Rt_complete;
  vector[final_time] SI;      //serial interval
  vector[final_time] f;      //infection ot observation distribution
  matrix[M, M] C;                // mobility matrix
  real<lower=0> iar;         // infection ascertainment ratio (reporting rate)
}

transformed data {
  vector[final_time] SI_rev = reverse(SI);      //reverse of the SI (required for simulation)
  vector[final_time] f_rev = reverse(f);       //reverse of the infection-observation distribution
}

parameters{

}

generated quantities{
  matrix[final_time, M] infection = rep_matrix(0, final_time, M);    // daily initialization
  matrix[final_time, M] cumm_sum = rep_matrix(0, final_time, M);     // total infection to date
  matrix[final_time, M] daily_deaths = rep_matrix(0,final_time, M);      // daily cases (infection becomes a case after latent period)
  matrix[final_time/7, M] weekly_deaths = rep_matrix(0, final_time/7, M); 
  
  //for (s in 1:n_simulation){
    matrix[final_time, M] Rt_regions = rep_matrix(Rt_complete, M);
  
    for (i in 1:M){                  // for initial seeding
      infection[1:seed_time,i] = rep_vector(init_seed[i],seed_time);      // learn the number of cases in the first initial seeding days
      cumm_sum[2:seed_time,i] = cumulative_sum(infection[2:seed_time,i]);    // cumulative infection
      daily_deaths[1:seed_time,i] = 1e-15 * infection[1:seed_time,i];
    }
    
    matrix[final_time,M] SI_regions = rep_matrix(SI_rev,M); 
 
///////////////
    for (t in (seed_time+1):final_time){ //for loop over time
      row_vector[M] convolution_inf = columns_dot_product(infection[1:(t-1),:] , SI_regions[(final_time-t+2):final_time,:]);  //infections at each region "k"
    
      vector[M] total_inf = C * convolution_inf';     //total infection at region "j"
      vector[M] eff_pop = C * pop;
      
      for (i in 1:M){      // for loop over region "i" (final infection at region "i")   
        real sus = pop[i]-sum(infection[1:(t-1),i]);
        infection[t,i] = dot_product(C[:,i]', (((rep_vector(sus , (M)) ./ eff_pop)' .* Rt_regions[t,:]).* total_inf'));
        daily_deaths[t,i] = iar * dot_product(infection[1:(t-1),i], tail(f_rev,t-1)); 
        cumm_sum[t,i] = cumm_sum[t-1,i] + infection[t-1,i];
      }
    }
    for (t in 1:(final_time/7)) {
     for (i in 1:M){
        weekly_deaths[t,i] = sum(daily_deaths[(7*(t-1)+1):(7*t),i]);
      }
    }
}
