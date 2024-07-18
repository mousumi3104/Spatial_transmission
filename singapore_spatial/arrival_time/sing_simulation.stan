// Reproduction number(Rt) flollows from the region where the population are having contact with the people

data {
  int<lower = 1> M_res;        //number of regions
  int<lower = 1> M_poi;
  vector[M_res] pop;                 //population of three regions
  int<lower=1> final_time;        // total time point to run the simulation
  int<lower = 0> seed_time;       //initial seeding time
  array[M_res] real init_seed;    //initial seeding
  matrix[final_time,M_res] Rt_res;
  matrix[final_time,M_poi] Rt_poi;
  matrix[M_res + M_poi, M_res] C_day;                // mobility matrix
  matrix[M_res + M_poi, M_res] C_end;
  vector[final_time] SI;      //serial interval
  vector[final_time] f;      //infection ot observation distribution
  real<lower=0> iar;         // infection ascertainment ratio (reporting rate)
  int<lower=0> n_simulation;
}

transformed data {
  vector[final_time] SI_rev = reverse(SI);      //reverse of the SI (required for simulation)
  vector[final_time] f_rev = reverse(f);       //reverse of the infection-observation distribution
  matrix[M_res, M_res] C_day_res = C_day[1:M_res,:];
  matrix[M_poi, M_res] C_day_poi = C_day[(M_res+1):(M_res+M_poi),:];
  matrix[M_res, M_res] C_end_res = C_end[1:M_res,:];
  matrix[M_poi, M_res] C_end_poi = C_end[(M_res+1):(M_res+M_poi),:];
}



//------------------------------------
    // if (sum(C_day_poi) == 0){ 
    // 
    // // infection only at residential place (in case of no movement)
    // 
    //   for (t in (seed_time+1):final_time){          //for loop over time
    //     
    //     row_vector[M_res] convolution_inf = columns_dot_product(infection[1:(t-1),:] , SI_res[(final_time-t+2):final_time,:]);  //infections at each region "k"
    // 
    //     int remainder = t % 7;
    //     if (remainder == 0 || remainder == 6){       // for weekends
    //       vector[M_res ] total_inf = C_end_res * convolution_inf';     //total infection at residential region "j" 
    //       vector[M_res ] eff_pop = C_end_res * pop;
    // 
    //       for (i in 1:M_res){      // for loop over region "i" (final infection at region "i")   
    //         
    //         real sus = pop[i]-sum(infection[1:(t-1),i]);
    //         infection[t,i] = dot_product(C_end_res[:,i]', (((rep_vector(sus , M_res) ./ eff_pop)' .* Rt[t,1:M_res]).* total_inf'));
    //         cases[t,i] = iar * dot_product(infection[1:(t-1),i], tail(f_rev,(t-1)));       
    //         cumm_sum[t,i] = cumm_sum[t-1,i] + infection[t-1,i];
    //       }
    //     }
    //     else{            // for weekdays
    //       vector[M_res ] total_inf = C_day_res * convolution_inf';     //total infection at residential region "j" 
    //       vector[M_res ] eff_pop = C_day_res * pop;
    //       
    //       for (i in 1:M_res){      // for loop over region "i" (final infection at region "i")   
    //         
    //         real sus = pop[i]-sum(infection[1:(t-1),i]);
    //         infection[t,i] = dot_product(C_day_res[:,i]', (((rep_vector(sus , M_res) ./ eff_pop)' .* Rt[t,1:M_res]).* total_inf'));
    //         cases[t,i] = iar * dot_product(infection[1:(t-1),i], tail(f_rev,(t-1)));       
    //         cumm_sum[t,i] = cumm_sum[t-1,i] + infection[t-1,i];
    //       }
    //     }
    //   }
    // }
    //else{
      
      
generated quantities{

  matrix[final_time, M_res] infection = rep_matrix(0, final_time, M_res);    // daily initialization
  matrix[final_time, M_res] cumm_sum = rep_matrix(0, final_time, M_res);     // total infection to date
  matrix[final_time, M_res] cases = rep_matrix(0,final_time, M_res);      // daily cases (infection becomes a case after latent period)
  matrix[final_time, M_res] daily_deaths = rep_matrix(0,final_time, M_res);
  //matrix[final_time, M_res] final_infection[n_simulation] ;     // total infection to date
  
  matrix[final_time, M_res+M_poi] Rt = append_col(Rt_res,Rt_poi);
  
    for (m in 1:M_res){                  // for initial seeding
      infection[1:seed_time,m] = rep_vector(init_seed[m],seed_time);      // learn the number of cases in the first initial seeding days
      cumm_sum[2:seed_time,m] = cumulative_sum(infection[2:seed_time,m]);    // cumulative infection
      daily_deaths[1:seed_time,m] = 1e-15 *infection[1:seed_time,m];
    }
    
    matrix[final_time,M_res] SI_res = rep_matrix(SI_rev,M_res); 
 
//---------------------------------------------------------------------------------------------------------//    

  for (t in (seed_time+1):final_time){ //for loop over time
    row_vector[M_res] convolution_inf = columns_dot_product(infection[1:(t-1),:] , SI_res[(final_time-t+2):final_time,:]);  //infections at each region "k"
    
    int remainder = t % 7;
    if (remainder == 0 || remainder == 6){       // for weekends
      vector[M_res + M_poi] total_inf = C_end * convolution_inf';     //total infection at residential region "j" 
      vector[M_res + M_poi] eff_pop = C_end * pop;
    
      for (m in 1:M_res){      // for loop over region "i" (final infection at region "i")   
        real sus = pop[m] - sum(infection[1:(t-1),m]);
        infection[t,m] = dot_product(C_end[:,m]', (((rep_vector(sus , (M_res + M_poi)) ./ eff_pop)' .* Rt[t,:]).* total_inf'));
        cases[t,m] = iar * dot_product(infection[1:(t-1),m], tail(f_rev,(t-1)));       
        cumm_sum[t,m] = cumm_sum[t-1,m] + infection[t-1,m];
      }
    }
    else{
      vector[M_res + M_poi] total_inf = C_day * convolution_inf';     //total infection at residential region "j" 
      vector[M_res + M_poi] eff_pop = C_day * pop;
      
      for (m in 1:M_res){      // for loop over region "i" (final infection at region "i")   
        real sus = pop[m]-sum(infection[1:(t-1),m]);
        infection[t,m] = dot_product(C_day[:,m]', (((rep_vector(sus , (M_res + M_poi)) ./ eff_pop)' .* Rt[t,:]).* total_inf'));
        cases[t,m] = iar * dot_product(infection[1:(t-1),m], tail(f_rev,(t-1)));       
        cumm_sum[t,m] = cumm_sum[t-1,m] + infection[t-1,m];
      }
    }
  }
}
