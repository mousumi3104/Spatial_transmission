data {
  int<lower=1> M_regions;
  int<lower=1> final_time;    //number of days to simulate
  int W;      // number weeks to simulate
  int<lower=1> initial_seeding_day;
  vector[M_regions] initial_seeding;
  vector[final_time]  SI;
  vector[final_time] f;
  vector[M_regions] pop;
  matrix[ final_time, M_regions] Rt;   //reproduction number
  matrix[final_time, 4] I;
  vector[M_regions] ifr_noise;
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

generated quantities{
   matrix[final_time, M_regions] infection;        // daily initialization
   matrix[W, M_regions] weekly_effect;

  
   matrix[final_time, M_regions] SI_regions = rep_matrix(SI_rev,M_regions);  // serial interval for every region
   matrix[final_time, M_regions] f_regions = rep_matrix(f_rev,M_regions);    // infection to death distribution for every region
    
  for (m in 1:M_regions){                  // for initial seeding
    infection[1:initial_seeding_day,m] = rep_vector(initial_seeding[m], initial_seeding_day);  //initial_seeding[m]
  }
   
  for (t in (initial_seeding_day + 1) : final_time){ //for loop over time
    vector[M_regions] convolution_inf = (columns_dot_product(infection[1:(t-1),:] , SI_regions[(final_time-t+2):final_time,:]))';  //infections at each region "k"
    for (m in 1:M_regions){
        
      real sus = pop[m] - sum(infection[1:(t-1),m]);
      infection[t,m] = (sus/pop[m]) * Rt[t,m] * convolution_inf[m];
    }
  }
}
  