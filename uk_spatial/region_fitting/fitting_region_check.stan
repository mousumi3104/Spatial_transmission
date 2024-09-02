data {
  int<lower=1> M_regions;
  int<lower=1> final_time;    //number of days to simulate
  int W;      // number weeks to simulate
  vector[final_time] mobility_change;
  int<lower=1> initial_seeding_day;
  int<lower=1> death_data_length;
  matrix[death_data_length , M_regions] death;       //number of infected individual at any time


  
}

transformed data {

}
parameters {
  row_vector<lower=0>[M_regions] mu;      //parameters for Rt
}

model{
  for (m in 1:M_regions){
      mu[m] ~ normal(3.28,0.2);
  }
}
