
Run_control = function(){
  pacman::p_load(cmdstanr, tidyverse,posterior, bayesplot, tidybayes, furrr, rstan, brms, faux,LRO.utilities,reticulate)
  
  source(here::here("control analysis", "Simulate_control_agents.R"))
  
  
  subjects = 30
  
  HR_intervention = seq(0,15,5)
  
  real_alpha_dif = seq(0,10,5)
  
  replicate = 1:100
  
  parameters = expand.grid(subjects = subjects,
                           HR_intervention = HR_intervention,
                           real_alpha_dif = real_alpha_dif,
                           replicate = replicate) %>% 
    mutate(id = 1:nrow(.))
  
  
  data_list <- split(parameters, parameters$id)
  
  
  cores = 10
  
  plan(multisession, workers = cores)
  
  # qq = sim_and_fit(data_list[[1]])
  
  possfit_model = possibly(.f = sim_and_fit, otherwise = "Error")
  
  results <- future_map(data_list, ~possfit_model(.x), .options = furrr_options(seed = TRUE))
  
  saveRDS(results,here::here("control analysis", "testing_real.RDS"))
  
}

Run_control()
