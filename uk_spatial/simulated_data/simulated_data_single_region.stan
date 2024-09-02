// Reproduction number(Rt) flollows from the region where the population are having contact with the people
data {
  int pop;                 //population of three regions
  int<lower = 1> final_time;        // total time point to run the simulation
  int<lower = 0> initial_seeding_day;
  int init_seed;       //initial seeding time    
  vector[final_time] SI;      //serial interval
  vector[final_time] f1;      //infection ot observation distribution
  vector[final_time] f2; 
  vector[final_time] Rt;
  real iar;
}

transformed data {
  vector[final_time] SI_rev = reverse(SI);      //reverse of the SI (required for simulation)
  vector[final_time] f1_rev = reverse(f1);       //reverse of the infection-death distribution
  vector[final_time] f2_rev = reverse(f2);       //reverse of the infection-onset distribution
}

parameters{

}

generated quantities{
  vector[final_time] infection = rep_vector(0, final_time);    // daily initialization
  vector[final_time] cumm_sum = rep_vector(0, final_time);     // total infection to date
  vector[final_time] daily_cases = rep_vector(0,final_time);
  vector[final_time] daily_deaths = rep_vector(0,final_time);      // daily cases (infection becomes a case after latent period)
  vector[final_time/7] weekly_deaths = rep_vector(0, final_time/7); 
  
      infection[1:initial_seeding_day] = rep_vector(init_seed,initial_seeding_day);      // learn the number of cases in the first initial seeding days
      cumm_sum[2:initial_seeding_day] = cumulative_sum(infection[2:initial_seeding_day]);    // cumulative infection
      daily_deaths[1:initial_seeding_day] = 1e-15 * infection[1:initial_seeding_day];
 
///////////////
    for (t in (initial_seeding_day+1):final_time){ //for loop over time

        real iar_noise = normal_rng(1,0.2);
        real ifr_noise = normal_rng(1,0.1);
        real sus = (pop-sum(infection[1:(t-1)]))/pop;
        
        infection[t] = sus * Rt[t] * dot_product(infection[1:(t-1)], tail(SI_rev,t-1)); 
        daily_cases[t] = iar * dot_product(infection[1:(t-1)], tail(f2_rev,t-1));
        daily_deaths[t] =  dot_product(infection[1:(t-1)], tail(f1_rev,t-1)); 
        cumm_sum[t] = cumm_sum[t-1] + infection[t-1];
      }
    for (t in 1:(final_time/7)) {
        weekly_deaths[t] = sum(daily_deaths[(7*(t-1)+1):(7*t)]);
     }
}
