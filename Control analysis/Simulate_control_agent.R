

#Idea is that the intervention changes the average HR of the participant.

#Therefore a significant difference in the threshold might still mean that the participant has the same belief about their HR

#So because a decrease in HR from the drugs and then a significant increase in the threshold is observed this could be excatly what one would expect
#even without a "real" difference from in belief in their HR. 

# We therefore need to first simulate that this can happen i.e. a significant threshold change based upon just a mean HR change of intervention

# after we have done this we need to then investigate that we can then observe a real latent effect even when this is present by "controling" for mean HR.


transform_data_to_stan = function(data){
  
  data = data %>% group_by(X,participant_id, sessions,alpha,beta, avg_HR)%>%
    summarise(X=mean(X),
              npx=n(),
              resp=sum(resp))%>%
    ungroup()
  
  return(data)
}


sim_and_fit = function(parameters){
  
  
  subjects = parameters$subjects
  
  HR_intervention = parameters$HR_intervention
  
  real_alpha_dif = parameters$real_alpha_dif
  
  alphas = data.frame(cbind(realalphas_ses1 = rnorm(subjects,-10,5),
                            realalphas_ses2 = rnorm(subjects,-10-real_alpha_dif,5)))
  
  real_thresholds = data.frame(cbind(alphas_ses1 = rnorm(subjects,55,5),
                                     alphas_ses2 = rnorm(subjects,55-real_alpha_dif,5)))
  
  
  betas = data.frame(cbind(beta_ses1 = rlnorm(subjects,2,0.2),
                           beta_ses2 = rlnorm(subjects,2,0.2)))
  
  
  lapses = data.frame(cbind(lapse_ses1 = brms::inv_logit_scaled(rnorm(subjects,-5,0.5)),
                            lapse_ses2 = brms::inv_logit_scaled(rnorm(subjects,-5,0.5))))
  
  
  avg_HR = data.frame(cbind(avg_HR_ses1 = rnorm(subjects,65,5),
                            avg_HR_ses2 = rnorm(subjects,65 - HR_intervention,5)))
  
  
  obs_alpha = data.frame(alpha_ses1 = avg_HR$avg_HR_ses1-real_thresholds$alphas_ses1,
                         alpha_ses2 = avg_HR$avg_HR_ses2-real_thresholds$alphas_ses2)
  
  
  df = data.frame(alphas,betas,lapses,avg_HR,obs_alpha) %>% mutate(participant_id = 1:n(), trials = 60) %>% 
    pivot_longer(
      cols = -c(participant_id,trials),
      names_to = c(".value", "sessions"),
      names_sep = "_ses"
    ) %>% arrange(sessions,participant_id) %>% mutate(sessions =as.numeric(sessions))
  
  dd = get_psi_stim(df)
  
  # dd %>% ggplot(aes(x = X, y = resp, col = as.factor(sessions)))+facet_wrap(~participant_id )+geom_point()
  
  data = inner_join(df,dd %>% dplyr::select(resp,X,participant_id,sessions,Estimatedthreshold, Estimatedslope,q5_threshold,q95_threshold,q5_slope,q95_slope), by = c("participant_id","sessions"))
  
  
  
  mod_noncent = cmdstanr::cmdstan_model(here::here("Manuscript","Control analysis","control stan","standard.stan"),stanc_options = list("O1"))
  
  data = data %>% 
    mutate(session = ifelse(sessions == 1, 0 ,1)) %>%
    mutate(id = parameters$id,
           subjects = parameters$subjects,
           HR_intervention = parameters$HR_intervention,
           real_alpha_dif = parameters$real_alpha_dif)
  
  data_fit = transform_data_to_stan(data)
  
  
  data_stan = list(T = nrow(data_fit),
                   S = length(unique(data_fit$participant_id)),
                   S_id = as.numeric(data_fit$participant_id ),
                   X = data_fit %>% .$X,
                   X_lapse = as.matrix(data.frame(int = rep(1,nrow(data_fit)))),
                   X_alpha = as.matrix(data.frame(int = rep(1,nrow(data_fit)),
                                                  session = data_fit %>% .$sessions)),
                   X_beta = as.matrix(data.frame(int = rep(1,nrow(data_fit)),
                                                 session = data_fit %>% .$sessions)),
                   N_alpha = 2,
                   N_beta = 2,
                   N_lapse = 1,
                   Y = data_fit %>% .$resp,
                   npx = data_fit %>% .$npx
  )
  
  
  
  
  #fitting
  fit_norm <- mod_noncent$sample(
    data = data_stan,
    iter_sampling = 1000,
    iter_warmup = 1000,
    chains = 4,
    parallel_chains = 4,
    refresh = 500,
    adapt_delta = 0.80,
    max_treedepth = 12
  )
  
  group_variables = fit_norm$summary("gm")
  
  group_variables = group_variables %>% mutate(variable = c("beta","beta_dif","lapse","alpha","alpha_dif"))%>%
    mutate(id = parameters$id,
           subjects = parameters$subjects,
           HR_intervention = parameters$HR_intervention,
           real_alpha_dif = parameters$real_alpha_dif)
  
  non_confits = fit_norm$summary(c("alpha_p","beta_p","lapse_p"))
  
  non_confits = non_confits %>% separate(variable, into = c("name", "indices"), sep = "\\[") %>%
    separate(indices, into = c("index1", "index2"), sep = ",|\\]") %>%
    mutate(subject = as.numeric(index1),
           sessions = as.numeric(index2))%>% select(name,mean,sd,q5,q95,subject,sessions)%>%
    mutate(id = parameters$id,subjects = parameters$subjects,
           HR_intervention = parameters$HR_intervention,
           real_alpha_dif = parameters$real_alpha_dif)
  
  
  diag_non_control = data.frame(fit_norm$diagnostic_summary())%>%
    mutate(id = parameters$id,
           subjects = parameters$subjects,
           HR_intervention = parameters$HR_intervention,
           real_alpha_dif = parameters$real_alpha_dif)
  
  mod_control = cmdstanr::cmdstan_model(here::here("Manuscript","Control analysis","control stan","control.stan"),stanc_options = list("O1"))
  
  
  data_stan = list(T = nrow(data_fit),
                   S = length(unique(data_fit$participant_id)),
                   S_id = as.numeric(data_fit$participant_id ),
                   X = data_fit %>% .$X,
                   X_lapse = as.matrix(data.frame(int = rep(1,nrow(data_fit)))),
                   X_alpha = as.matrix(data.frame(int = rep(1,nrow(data_fit)),
                                                  session = data_fit %>% .$sessions,
                                                  avghr = scale(data_fit$avg_HR)[,1])),
                   X_beta = as.matrix(data.frame(int = rep(1,nrow(data_fit)),
                                                 session = data_fit %>% .$sessions,
                                                 avghr = scale(data_fit$avg_HR)[,1])),
                   N_alpha = 3,
                   N_beta = 3,
                   N_lapse = 1,
                   Y = data_fit %>% .$resp,
                   npx = data_fit %>% .$npx
  )
  
  
  
  
  #fitting
  fit_control <- mod_control$sample(
    data = data_stan,
    iter_sampling = 1000,
    iter_warmup = 1000,
    chains = 4,
    init = 0,
    parallel_chains = 4,
    refresh = 50,
    adapt_delta = 0.9,
    max_treedepth = 12
  )
  
  diag_control = data.frame(fit_control$diagnostic_summary())%>%
    mutate(id = parameters$id,
           subjects = parameters$subjects,
           HR_intervention = parameters$HR_intervention,
           real_alpha_dif = parameters$real_alpha_dif)
  
  
  group_variables_con = fit_control$summary("gm")
  
  group_variables_con = group_variables_con %>% mutate(variable = c("alpha","alpha_dif","alpha_con","beta","beta_dif","beta_con","lapse"))%>% 
    mutate(id = parameters$id,
           subjects = parameters$subjects,
           HR_intervention = parameters$HR_intervention,
           real_alpha_dif = parameters$real_alpha_dif)
  
  
  confits = fit_control$summary(c("alpha_p","beta_p","lapse_p"))
  
  confits = confits %>% separate(variable, into = c("name", "indices"), sep = "\\[") %>%
    separate(indices, into = c("index1", "index2"), sep = ",|\\]") %>%
    mutate(subject = as.numeric(index1),
           sessions = as.numeric(index2)) %>% select(name,mean,sd,q5,q95,subject,sessions)%>% 
    mutate(id = parameters$id,
           subjects = parameters$subjects,
           HR_intervention = parameters$HR_intervention,
           real_alpha_dif = parameters$real_alpha_dif)
  
  
  real_simulated = df %>% group_by(sessions) %>% summarize(mean_realalpha = mean(realalphas),se_realalpha = sd(realalphas)/sqrt(n()),
                                                           mean_beta = mean(beta),se_beta = sd(beta)/sqrt(n()),
                                                           mean_lapse = mean(lapse), se_lapse = sd(lapse)/sqrt(n()),
                                                           mean_avg_HR = mean(avg_HR),se_avg_HR = sd(avg_HR)/sqrt(n()),
                                                           mean_alpha = mean(alpha), se_alpha = sd(alpha)/sqrt(n()))%>% 
    mutate(id = parameters$id,
           subjects = parameters$subjects,
           HR_intervention = parameters$HR_intervention,
           real_alpha_dif = parameters$real_alpha_dif)
  
  return(list(sim = real_simulated,
              controlled_fits = confits, group_controlled_fits = group_variables_con,diag_control = diag_control,
              uncontrolled_fits = non_confits, group_uncontrolled_fits = group_variables,diag_non_control = diag_non_control,
              full_data = data))  
  
}


get_psi_stim = function(parameters){
  
  python_script <- here::here("Control analysis","PSI2.py")
  
  alpha = parameters$alpha
  
  beta = parameters$beta
  
  lapse = parameters$lapse
  
  trials = as.integer(parameters$trials)
  
  ids = as.integer(parameters$participant_id)
  
  subjects = max(as.integer(unique(parameters$participant_id)))
  
  sessions = max(as.integer(unique(parameters$sessions)))
  
  library(reticulate)
  
  # Use reticulate to run the Python script with arguments
  
  source_python(python_script, convert = FALSE)
  
  d = get_stim(lapse, alpha, beta, ids, trials, subjects, sessions)
  
  #Produces a warning about will be removed in future versions (but seems to be a thing with reticulate and not my code)
  dd = reticulate::py_to_r(d)
  dd = reticulate::py_to_r(d)
  
  d = data.frame(lapse = dd$lapse, alpha = dd$alpha, beta = dd$beta, participant_id = unlist(dd$participant_id),
                 trials = unlist(dd$trials), subs = unlist(dd$subs), X = unlist(dd$X), resp = unlist(dd$resp), sessions = unlist(dd$sessions),
                 Estimatedthreshold = unlist(dd$Estimatedthreshold), Estimatedslope = unlist(dd$Estimatedslope), q5_threshold =  unlist(dd$q5_threshold),
                 q95_threshold =  unlist(dd$q95_threshold), q5_slope =  unlist(dd$q5_slope), q95_slope =  unlist(dd$q95_slope))
  
  
  #plot for it to make sense:
  # d  %>%  group_by(participant_id,sessions) %>% mutate(trial = 1:n()) %>%
  #   ggplot(aes(x = trial, y = Estimatedslope, col = sessions))+geom_pointrange(aes(ymin = q5_slope, ymax = q95_slope))+
  #   facet_wrap(~participant_id)+geom_hline(aes(yintercept = beta, col = sessions))
  # 
  # 
  # d  %>%  group_by(participant_id,sessions) %>% mutate(trial = 1:n()) %>%
  #   ggplot(aes(x = trial, y = Estimatedthreshold, col = sessions))+geom_pointrange(aes(ymin = q5_threshold, ymax = q95_threshold))+
  #   facet_wrap(~participant_id)+geom_hline(aes(yintercept = alpha, col = sessions))
  # 
  
  return(d)
  
}





#function to get confidence intervals used in plot_intervals.
ci = function(x){
  list = list(which(cumsum(x)/sum(x) > 0.025)[1],
              last(which(cumsum(x)/sum(x) < 0.975)))
  return(list)
  
}
