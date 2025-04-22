
data {
  int<lower=1> M_regions;
  int<lower=1> final_time;    //number of days to simulate
  int<lower=1> W;
  matrix[final_time, M_regions] gmobility;
  int<lower=1> initial_seeding_day;
  int<lower=1> death_data_length;
  array[death_data_length, M_regions] int<lower=0> death;       //number of infected individual at any time
  vector[final_time]  SI;
  vector[final_time] f;
  vector[M_regions] pop;
  array[final_time] int day_week_index;
  matrix[final_time, 4] I;
  int fitting_start;
}

transformed data {
  vector[final_time] SI_rev;  // SI in reverse order
  vector[final_time] f_rev;   // f in reversed order
  for(i in 1:final_time){
    SI_rev[i] = SI[final_time-i+1];
  }
  for(i in 1:final_time) {
     f_rev[i] = f[final_time-i+1];
    }
}

parameters {
  vector<lower=0>[M_regions] mu;      //parameters for Rt
  array[M_regions] real <lower=0> initial_seeding;
  real<lower=0> tau;
  real<lower=0> weekly_var;       //weekly variance
  matrix[W, M_regions] weekly_effect_d;     //parameters for Rt (why W+1).     ?????
  real<lower=0> phi1;
  vector<lower=0>[M_regions] ifr_noise;
  real<lower=0> gamma;
  vector[M_regions] x;
  vector[M_regions] y;
  vector[M_regions] z;
}

transformed parameters{
  matrix[final_time, M_regions] infection ;    // daily initialization
  matrix[final_time, M_regions] daily_deaths ;      // daily deaths (infection becomes a case after latent period)
  matrix[death_data_length, M_regions] weekly_deaths ;
  matrix[W, M_regions] weekly_effect ;   // check why this +1 is needed 
  matrix[final_time, M_regions] Rt ;   //reproduction number
  
  //////////////////////////////////////////
  matrix[final_time, M_regions] SI_regions = rep_matrix(SI_rev,M_regions);  // serial interval for every region
  matrix[final_time, M_regions] f_regions = rep_matrix(f_rev,M_regions);    // infection to death distribution for every region

 // for initial seeding
 
 for (m in 1:M_regions){
   
    infection[1:initial_seeding_day,m] = rep_vector(initial_seeding[m], initial_seeding_day);      // learn the number of cases in the first initial seeding days
    weekly_effect[:,m] = weekly_var * cumulative_sum(weekly_effect_d[:,m]) ;    // weekly effect
    Rt[1:initial_seeding_day,m] = rep_vector( 2 * mu[m]* inv_logit(- weekly_effect[1,m]) , initial_seeding_day);   // why this is weekly_effect[m]??
 }
  
  for (t in (initial_seeding_day+1):final_time){ 
    vector[M_regions] convolution_inf = (columns_dot_product(infection[1:(t-1),:] , SI_regions[(final_time-t+2):final_time,:]))';	
    
		Rt[t,:] = 2 * (mu' .* inv_logit(- weekly_effect[day_week_index[t],:]
                                 - ((x' .* gmobility[t,:]) * I[t,2])
				                         - ((y' .* gmobility[t,:]) * I[t,3])
				                         - ((z' .* gmobility[t,:]) * I[t,4])));
				                         
			
		vector[M_regions] sus = pop -  (rep_row_vector(1, t-1) * (infection[1:(t-1),:]))';	
		infection[t,:] = ((sus ./ pop) .* Rt[t,:]' .* convolution_inf)' ;
  }
  
  
  daily_deaths[1] = 1e-15 * infection[1];
  for (t in 2:final_time){
    daily_deaths[t,:] =  ifr_noise' .* columns_dot_product(infection[1:(t-1),:], f_regions[(final_time-t+2):final_time,:]); 
    }
  
  for (t in 1:(final_time %/%7)) {        //weekly_deaths
    weekly_deaths[t,:] = rep_row_vector(1,7) * daily_deaths[(7*(t-1)+1):7*t,:]; 
    }
  //   weekly_deaths[(final_time %/% 7)+1] = sum(daily_deaths[(7*(final_time %/% 7)+1) : final_time]);
}

model {
 phi1 ~ normal(0,5);
 tau ~ exponential(0.01);
 initial_seeding ~ exponential(1/tau);
 mu ~ normal(3.28,0.5)T[0,];
 gamma ~ normal(0,0.5);
 ifr_noise ~ normal(1,0.1);
 
 x ~ normal(0,gamma);
 y ~ normal(0,gamma);
 z ~ normal(0,gamma);

 weekly_var ~ normal(0,.2);
 to_vector(weekly_effect_d) ~ normal(0,1);
 
 for (week in fitting_start : death_data_length){
    target += neg_binomial_2_lpmf(death[week,]| weekly_deaths[week,], phi1);
 }
}
