data {
  int<lower =1> M;  //number of regions
  vector[M] pop;             //population of three regions
  int<lower = 1> final_time;        // total time point to run the simulation
  int<lower = 0> seed_time;       //initial seeding time
  matrix[M,M] C;              //initialize the susceptible population
  matrix[final_time,M] Rt;       //reproduction number
  vector[final_time] SI;     //serial interval
  vector[final_time] f;      //infection ot observation distribution
  real<lower=1> iar;      // infection ascertainment ratio (reporting rate)
  real init_seed[M];      //initial seeding
}
transformed data {
  vector[final_time] SI_rev = reverse(SI);      //reverse of the SI (required for simulation)
  vector[final_time] f_rev = reverse(f);       //reverse of the infection-observation distribution
}
model{
  //this is not required
}
generated quantities{
  matrix[final_time, M] infection = rep_matrix(0,final_time,M);    //initialization
  //matrix[final_time, M] cases  = rep_matrix(0,final_time,M);
  //matrix[final_time,M] cumm_sum = rep_matrix(0,final_time,M);
  
  for (i in 1:M){           // for initial seeding
    infection[1:seed_time,i] = rep_vector(init_seed[i],seed_time);      // learn the number of cases in the first initial seeding days
    //cumm_sum[2:seed_time,i] = cumulative_sum(infection[2:seed_time,i]);    // cumulative infection
  }
  
  for (t in (seed_time+1):final_time){         //for loop over time
    
    matrix[final_time,M] SI_rep = rep_matrix(SI_rev,M);      // infections at region "k"
    row_vector[M] convolution = columns_dot_product(infection[1:t-1,:] , SI_rep[final_time-t+2:final_time,:]);  //infections at each region "k"
    vector[M] total_inf =  C * convolution';     //total infection at region "j" 
    vector[M] eff_pop = C * pop;
    
    for (i in 1:M){                   // for loop over region "i" (final infection at region "i")
      real sus = pop[i]-sum(infection[1:(t-1),i]);
      infection[t,i] = dot_product(C[1:M,i]', (((rep_vector(sus , M) ./ eff_pop)' .* Rt[t,1:M]) .* total_inf'));
      //cumm_sum[t,i] = cumm_sum[t-1,i] + infection[t-1,i];
    }
  }
}
