data {
  int<lower = 1> M_res;        //number of regions
  int<lower = 1> M_poi;
  vector[M_res] pop;                 //population of three regions
  int<lower = 1> final_time;        // total time point to run the simulation
  int<lower = 0> seed_time;       //initial seeding time
  vector[M_res] Rt_res;
  vector[M_poi] Rt_poi;
  real init_seed[M_res];              //initial seeding
  matrix[M_res + M_poi, M_res] C;                // mobility matrix
  vector[final_time] SI;      //serial interval
  vector[final_time] f;      //infection ot observation distribution
  real<lower=0> iar;         // infection ascertainment ratio (reporting rate)
  int<lower=0> n_simulation;
}

transformed data {
  vector[final_time] SI_rev = reverse(SI);      //reverse of the SI (required for simulation)
  vector[final_time] f_rev = reverse(f);       //reverse of the infection-observation distribution
  matrix[M_res, M_res] C_res= C[1:M_res,:];
  matrix[M_poi, M_res] C_poi= C[(M_res+1):(M_res + M_poi),:];
}

parameters{
  //real<lower=0.7,upper=1.2> Rt_res_data;
  //real<lower=2.7,upper=3.2> Rt_poi_data;
}

generated quantities{

  matrix[final_time, M_res] infection = rep_matrix(0, final_time, M_res);    // daily initialization
  matrix[final_time, M_res] cumm_sum = rep_matrix(0, final_time, M_res);     // total infection to date
  matrix[final_time, M_res] cases = rep_matrix(0,final_time, M_res);      // daily cases (infection becomes a case after latent period)
  matrix[final_time, M_res] final_infection[n_simulation] ;     // total infection to date
  matrix[final_time, M_res+M_poi] final_Rt[n_simulation];      // daily cases (infection becomes a case after latent period)
  
  
  for (s in 1:n_simulation){
    row_vector[ M_res ] Rt_res_data1; 
    row_vector[ M_poi ] Rt_poi_data1;
    
  for (j in 1:M_res){
        Rt_res_data1[j] = normal_rng(Rt_res[j],0.2); 
      }
      for (j in 1:M_poi){
        Rt_poi_data1[j] = normal_rng(Rt_poi[j],0.2); 
      }
      matrix[final_time, M_res+M_poi] Rt = append_col(rep_matrix(Rt_res_data1,final_time),rep_matrix(Rt_poi_data1,final_time));  
    
    for (i in 1:M_res){                  // for initial seeding
      infection[1:seed_time,i] = rep_vector(init_seed[i],seed_time);      // learn the number of cases in the first initial seeding days
      cumm_sum[2:seed_time,i] = cumulative_sum(infection[2:seed_time,i]);    // cumulative infection
    }
    matrix[final_time,M_res] SI_res = rep_matrix(SI_rev,M_res); 
 
/////////////////// for inter area mobility (connection matrix is the identity matrix) ///////////////////////////////////////         
    if (sum(C[(M_res+1):(M_res+M_poi)]) == 0){ 
    
    // infection only at residential place (in case of no movement)
    
    for (t in (seed_time+1):final_time){ //for loop over time
    
      row_vector[M_res] convolution_res = columns_dot_product(infection[1:(t-1),:] , SI_res[(final_time-t+2):final_time,:]);  //infections at each region "k"
    
      vector[M_res + M_poi] total_inf = C * convolution_res';     //total infection at residential region "j" 
      vector[M_res + M_poi] eff_pop = C * pop;
  
      for (i in 1:M_res){      // for loop over region "i" (final infection at region "i")   
        real sus = pop[i]-sum(infection[1:(t-1),i]);
        infection[t,i] = dot_product(C[1:M_res,i]', (((rep_vector(sus , M_res) ./ eff_pop[1:M_res])' .* Rt[t,1:M_res]).* total_inf[1:M_res]'));
        
        cases[t,i] = iar * dot_product(infection[1:(t-1),i], tail(f_rev,(t-1)));       
        cumm_sum[t,i] = cumm_sum[t-1,i] + infection[t-1,i];
      }
    }
  }
//////////////////////////////////////////////////////////////////////////////////////      
  else{
    
  for (t in (seed_time+1):final_time){ //for loop over time
    row_vector[M_res] convolution_inf = columns_dot_product(infection[1:(t-1),:] , SI_res[(final_time-t+2):final_time,:]);  //infections at each region "k"
    
    vector[M_res + M_poi] total_inf = C * convolution_inf';     //total infection at residential region "j" 
    vector[M_res + M_poi] eff_pop = C * pop;
    
    
    for (i in 1:M_res){      // for loop over region "i" (final infection at region "i")   
      real sus = pop[i]-sum(infection[1:(t-1),i]);
      infection[t,i] = dot_product(C[:,i]', (((rep_vector(sus , (M_res + M_poi)) ./ eff_pop)' .* Rt[t,:]).* total_inf'));
      
      cases[t,i] = iar * dot_product(infection[1:(t-1),i], tail(f_rev,(t-1)));       
      cumm_sum[t,i] = cumm_sum[t-1,i] + infection[t-1,i];
    }
  }
}
final_infection[s,,]= infection;
final_Rt[s,,]=Rt;
}
}
