data {
  int<lower=1> M_regions;
  int<lower=1> final_time;    //number of days to simulate
  int W;      // number weeks to simulate
  matrix[final_time, M_regions] gmobility;
  int<lower=1> initial_seeding_day;
  vector[M_regions] initial_seeding;
  vector[final_time]  SI;
  vector[final_time] f;
  vector[M_regions] pop;
  matrix[M_regions,M_regions] C_base;
  matrix[M_regions,M_regions] C_lockdown;
  matrix[ final_time, M_regions] Rt;   //reproduction number
  matrix[final_time, 4] I;
  vector[M_regions] ifr_noise;
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
}

generated quantities{
   
   matrix[final_time, M_regions] infection;    // daily initialization
   matrix[final_time, M_regions] infection_in_own;
   matrix[final_time, M_regions] infection_in_mob ;
   matrix[final_time, M_regions] infection_out_mob;
   matrix[final_time, M_regions] daily_deaths;      // daily deaths (infection becomes a case after latent period)
   matrix[W, M_regions] weekly_effect;
   matrix[(final_time %/% 7) +1, M_regions] weekly_deaths;
   
   //-------------------------------------------------------------------------------------------------//
   matrix[final_time, M_regions] SI_regions = rep_matrix(SI_rev,M_regions);  // serial interval for every region
   matrix[final_time, M_regions] f_regions = rep_matrix(f_rev,M_regions);    // infection to death distribution for every region
    
   for (m in 1:M_regions){                  // for initial seeding
      infection[1:initial_seeding_day,m] = rep_vector(initial_seeding[m], initial_seeding_day);  //initial_seeding[m]
      infection_in_own[1:initial_seeding_day,m] = rep_vector(initial_seeding[m],initial_seeding_day);
      infection_in_mob[1:initial_seeding_day,m] = rep_vector(0,initial_seeding_day);
      infection_out_mob[1:initial_seeding_day,m] = rep_vector(0,initial_seeding_day);
    }
   
   for (t in (initial_seeding_day + 1) : final_time){ //for loop over time
      vector[M_regions] convolution_inf = (columns_dot_product(infection[1:(t-1),:] , SI_regions[(final_time-t+2):final_time,:]))';  //infections at each region "k"
      matrix[M_regions,M_regions] C;
      
      if (I[t,1] == 0){
        C = C_base ;
      }else{
        C = C_lockdown ;
      }
      
      vector[M_regions] total_inf = C * convolution_inf;     //total infection at region "j"
      vector[M_regions] eff_pop = C * pop;
      
      for (m in 1:M_regions){
        
        real sus = pop[m] - sum(infection[1:(t-1),m]);
        infection[t,m] = dot_product(C[:,m]', (((rep_vector(sus , M_regions) ./ eff_pop)'.* Rt[t,:]).* total_inf'));
        infection_in_own[t,m] = ((sus * C[m,m]) / eff_pop[m]) * Rt[t,m] * (C[m,m] * convolution_inf[m]) ;
        infection_in_mob[t,m] = ((sus * C[m,m]) / eff_pop[m]) * Rt[t,m] * (total_inf[m] - (C[m,m] * convolution_inf[m])) ;
        array[M_regions - 1] int index; 
        int ind = 1;
        
        for (n  in 1:M_regions){
          if (n != m){
            
            index[ind] = n;
            ind += 1;
          }
        }
      infection_out_mob[t,m] = dot_product(C[index,m]', (((rep_vector(sus , M_regions-1) ./ eff_pop[index])' .* Rt[t,index]).* total_inf[index]'));
      }
    }
    daily_deaths[1,:] = 1e-15 * infection[1,:];
     for (t in 2:final_time){
        daily_deaths[t,:] = (ifr_noise)' .* columns_dot_product(infection[1:(t-1),:] , f_regions[(final_time-t+2):final_time,:]);
      }
  
   for (m in 1:M_regions){
     for (t in 1:(final_time %/% 7)) { //weekly_deaths
        weekly_deaths[t,m] = sum(daily_deaths[(7*(t-1)+1):7*t,m]);
    }
    weekly_deaths[(final_time %/% 7)+1,m] = sum(daily_deaths[(7*(final_time %/% 7)+1) : final_time,m]);
  }
}
