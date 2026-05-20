
## PSI trajectories
plot_1 = function(){
  
  ## helper functions:
  # function to get confidence intervals used in plot_intervals.
  ci <- function(x) {
    list <- list(
      which(cumsum(x) / sum(x) > 0.025)[1],
      last(which(cumsum(x) / sum(x) < 0.975))
    )
    return(list)
  }
  
  # get the lines and the confidence intervals
  # function that uses ci()
  
  get_line_intervals <- function(data) {
    # lower and lower conf. interval
    upper <- array(NA, np$size(data[, 1, 1]))
    lower <- array(NA, np$size(data[, 1, 1]))
    
    # loop through the np array
    for (i in 1:np$size(data[, 1, 1])) {
      
      confidence <- ci(rowMeans(data, dims = 2)[i, ])
      rg <- seq(-50.5, 50.5, by = 1)
      upper[i] <- rg[confidence[[1]]]
      lower[i] <- rg[confidence[[2]]]
    }
    
    data <- data.frame(upper = upper, lower = lower, x = seq(0, nrow(upper), length.out = nrow(upper)))
    
    return(data)
  }
  
  
  # function that uses get_line_intervals(.)
  plot_interval <- function(df, interoPost = NA) {
    d <- data.frame()
    dd <- data.frame()
    
    
    if (!is.na(interoPost)[1] == TRUE) {
      dd <- get_line_intervals(interoPost)
    }
    d <- rbind(d, dd) %>% mutate(trials = 1:n())
    
    dataline <- df %>%
      filter(TrialType == "psi") %>%
      filter(Modality == "Intero") %>%
      mutate(x = seq(0, nrow(.), length.out = (nrow(.))))
    
    
    # making trials go from 0-60 in each modality
    df <- df %>%
      group_by(Modality) %>%
      mutate(trials = 1:n()) %>%
      ungroup()
    
    q = inner_join(d,df, by = "trials")
    return(q)
  }
  
  
  # function that uses plot_interval(.)
  intervals = function(df, interoPost = NA){
    
    df$ResponseCorrect <- ifelse(df$ResponseCorrect == "", NA, df$ResponseCorrect)
    df$Decision <- ifelse(df$Decision == "", NA, df$Decision)
    
    df <- df %>%
      mutate(
        Decision = as.character(df$Decision),
        ConfidenceRT = as.numeric(df$ConfidenceRT),
        DecisionRT = as.numeric(df$DecisionRT),
        Confidence = as.numeric(df$Confidence),
        Condition = as.character(df$Condition),
        listenBPM = as.numeric(df$listenBPM),
        responseBPM = as.numeric(df$responseBPM),
        ResponseCorrect = as.numeric(as.factor(df$ResponseCorrect)) - 1,
        EstimatedThreshold = as.numeric(df$EstimatedThreshold),
        EstimatedSlope = as.numeric(df$EstimatedSlope),
      )
    
    # check for NA's in the critical columns of the data:
    # specify different ways of coding missing values:
    missing_pres <- c(NA, "na", "N/A", "NaN", NaN, "n/a")
    
    trials_missing <- df %>%
      select(Decision, Confidence, ConfidenceRT, Confidence, Condition, listenBPM, responseBPM, nTrials) %>%
      filter_all(any_vars(. %in% missing_pres)) %>%
      .$nTrials
    
    
    
    if (length(trials_missing) != 0) {
      print(paste("Number of NA's = ", length(trials_missing), " detected in trials : "))
      print(as.character(trials_missing))
    }
    
    # remove the NA's
    df1 <- df %>% filter(!nTrials %in% trials_missing)
    
    interval <- plot_interval(df1, interoPost)
    
    return(interval)
    
  }
  
  
  #psychometric equation HRD.
  psycho_fit = function(x, alpha,lapse,beta){
    
    return(lapse + (1 - 2 * lapse) * (0.5 + 0.5 *pracma::erf((x-alpha)/(beta*sqrt(2)))))
  }
  
  #getting indvidual subject posterior fits HRD (for the choosen subject)
  get_indi = function(model,draw_id){
    
    individual_params = as_draws_df(model$draws(paste0("param[11,",1:11,"]")))
    
    alpha_beta_lapse = individual_params %>% 
      select(-contains(".")) %>% 
      mutate(draw = 1:n()) %>% 
      pivot_longer(-draw, names_to = "variable") %>% 
      filter(draw %in% draw_id) %>% 
      mutate(parameter = as.numeric(sub(".*,(\\d+)\\]", "\\1", variable))) %>% 
      mutate(subject = as.numeric(gsub("^param\\[(\\d+),\\d+\\]", "\\1", variable))) %>% 
      mutate(parameter = ifelse(parameter == 1, "alpha_placebo1",
                                ifelse(parameter == 2, "alpha_biso_dif",
                                       ifelse(parameter == 3, "alpha_prop_dif",
                                              ifelse(parameter == 4, "alpha_control",
                                                     ifelse(parameter == 5, "beta_placebo1",
                                                            ifelse(parameter == 6, "beta_biso_dif",
                                                                   ifelse(parameter == 7, "beta_prop_dif",
                                                                          ifelse(parameter == 8, "beta_control",
                                                                                 ifelse(parameter == 9, "lapse_placebo1",
                                                                                        ifelse(parameter == 10, "lapse_biso_dif",
                                                                                               ifelse(parameter == 11, "lapse_prop_dif",NA)))))))))))) %>%
      select(parameter,draw,value) %>% pivot_wider(names_from = "parameter", values_from = "value") %>% 
      mutate(alpha_placebo = alpha_placebo1,
             alpha_biso = alpha_placebo1+alpha_biso_dif,
             alpha_prop = alpha_placebo1+alpha_prop_dif,
             beta_placebo= exp(beta_placebo1),
             beta_biso = exp(beta_placebo1+beta_biso_dif),
             beta_prop = exp(beta_placebo1+beta_prop_dif),
             lapse_placebo = brms::inv_logit_scaled(lapse_placebo1)/2,
             lapse_biso = brms::inv_logit_scaled(lapse_placebo1+lapse_biso_dif)/2,
             lapse_prop = brms::inv_logit_scaled(lapse_placebo1+lapse_prop_dif)/2
      ) %>% select(-contains("dif")) %>% 
      dplyr::select(contains(c("placebo","biso","prop","draw"))) %>% select(-c("alpha_placebo1","beta_placebo1","lapse_placebo1")) %>% 
      pivot_longer(cols = starts_with(c("alpha", "beta", "lapse")), 
                   names_to = c(".value", "status"), 
                   names_sep = "_")
    
    
    
    individual_means = alpha_beta_lapse %>% 
      ungroup() %>% 
      rowwise() %>%
      mutate(prob = list(psycho_fit(seq(-35,55,by = 0.1),alpha, lapse,beta)),
             Alpha = list(seq(-35,55,by = 0.1))) %>% 
      unnest() %>% mutate(status = ifelse(status == "placebo","placebo",ifelse(status == "biso","bisoprolol","propranolol"))) %>% 
      mutate(status = factor(status, levels = c("bisoprolol","placebo", "propranolol")))
    
    
  }
  
  # psychometric equation (RRST)
  psycho_fit_rrst = function(x, alpha,lapse,beta){
    return(0.5 + (1 - 0.5 - lapse) * (1-exp(-10^(beta * (x-alpha)))))
  }
  
  #getting indvidual subject posterior fits RRST (for the choosen subject)
  get_indi_rrst = function(model,draw_id){
    
    individual_params = as_draws_df(model$draws(paste0("param[11,",1:9,"]")))
    
    
    individual_means = individual_params %>% 
      select(-contains(".")) %>% 
      mutate(draw = 1:n()) %>% 
      pivot_longer(-draw, names_to = "variable") %>% 
      filter(draw %in% draw_id) %>% 
      mutate(parameter = as.numeric(sub(".*,(\\d+)\\]", "\\1", variable))) %>% 
      mutate(subject = as.numeric(gsub("^param\\[(\\d+),\\d+\\]", "\\1", variable))) %>% 
      mutate(parameter = ifelse(parameter == 1, "beta_placebo1",
                                ifelse(parameter == 2, "beta_biso_dif",
                                       ifelse(parameter == 3, "beta_prop_dif",
                                              ifelse(parameter == 4, "lapse_placebo1",
                                                     ifelse(parameter == 5, "lapse_biso_dif",
                                                            ifelse(parameter == 6, "lapse_prop_dif",
                                                                   ifelse(parameter == 7, "alpha_placebo1",
                                                                          ifelse(parameter == 8, "alpha_biso_dif",
                                                                                 ifelse(parameter == 9, "alpha_prop_dif",NA)))))))))) %>%
      select(parameter,draw,value) %>% pivot_wider(names_from = "parameter", values_from = "value") %>% 
      mutate(beta_placebo= exp(beta_placebo1),
             beta_biso = exp(beta_placebo1+beta_biso_dif),
             beta_prop = exp(beta_placebo1+beta_prop_dif),
             lapse_placebo = brms::inv_logit_scaled(lapse_placebo1)/2,
             lapse_biso = brms::inv_logit_scaled(lapse_placebo1+lapse_biso_dif)/2,
             lapse_prop = brms::inv_logit_scaled(lapse_placebo1+lapse_prop_dif)/2,
             alpha_placebo = (alpha_placebo1),
             alpha_biso = (alpha_placebo1+alpha_biso_dif),
             alpha_prop = (alpha_placebo1+alpha_prop_dif))%>% 
      select(-contains("dif")) %>% 
      dplyr::select(contains(c("placebo","biso","prop","draw"))) %>% 
      select(-c("beta_placebo1","lapse_placebo1","alpha_placebo1")) %>% 
      pivot_longer(cols = starts_with(c("alpha", "beta", "lapse")), 
                   names_to = c(".value", "status"), 
                   names_sep = "_")
    
    
    
    
    individual_means = individual_means %>% 
      ungroup() %>% 
      rowwise() %>%
      mutate(prob = list(psycho_fit_rrst(seq(-2.5,22.5,by = 0.1),alpha,lapse,beta)),
             Alpha = list(seq(-2.5,22.5,by = 0.1))) %>% 
      unnest() %>% mutate(status = ifelse(status == "placebo","placebo",ifelse(status == "biso","bisoprolol","propranolol"))) %>% 
      mutate(status = factor(status, levels = c("bisoprolol","placebo", "propranolol")))
    
    return(individual_means)
    
  }
  
  
  ############# HRD plot1
  plot1_hrd = function(){
    
    # load subjet data for first visist
    df = read_delim(here::here("data","sub_3025","visit_0001",
                               "InteroceptionTasks","HRD","0251HRD","0251HRD_final.txt"))
    # load subjet numpy array (PSI trajectory) for first visist    
    interoPost = np$load(here::here("data","sub_3025","visit_0001",
                                    "InteroceptionTasks","HRD","0251HRD","0251Intero_posterior.npy"))
    # load subjet data for second visit
    df1 = read_delim(here::here("data","sub_3025","visit_0002",
                                "InteroceptionTasks","HRD","0252HRD","0252HRD_final.txt"))
    
    # load subjet numpy array (PSI trajectory) for second visist    
    interoPost1 = np$load(here::here("data","sub_3025","visit_0002",
                                     "InteroceptionTasks","HRD","0252HRD","0252Intero_posterior.npy"))
    # load subjet data for last visist
    df2 = read_delim(here::here("data","sub_3025","visit_0003",
                                "InteroceptionTasks","HRD","0253HRD","0253HRD_final.txt"))
    
    # load subjet numpy array (PSI trajectory) for last visist    
    interoPost2 = np$load(here::here("data","sub_3025","visit_0003",
                                     "InteroceptionTasks","HRD","0253HRD","0253Intero_posterior.npy"))
    
    #get the intervals for each of the vists (this subject had Bisoprolol, propanolol and lastly placebo).
    v1 = intervals(df,interoPost)%>% mutate(drug = "Biso")
    v2 = intervals(df1,interoPost1)%>% mutate(drug = "Prop")
    v3 = intervals(df2,interoPost2)%>% mutate(drug = "Placebo")
    
    
    
    #get the subjects fitted psychometric:
    model <- readRDS(here::here("STAN models","HRD.RDS"))
    
    draw_id = sample(1:8000,50)
    
    # get this subjects individual trajectories.
    individual_means = get_indi(model,draw_id)
    
    ##################  all three conditions plot:
    plot_combined = rbind(v1,v2,v3) %>% ggplot(aes()) +
      geom_point(aes(x = trials, y = Alpha, color = drug, shape = Decision), size = 2.5) +
      # shapes' shape
      scale_shape_manual(values = c(0, 16)) +
      # 0 line
      geom_hline(yintercept = 0, linetype = "dashed") +
      # line that is inside Confidence interval
      geom_line(aes(x = x, y = EstimatedThreshold, color = drug)) +
      # confidence interval
      geom_ribbon(aes(x = x, ymin = lower, ymax = upper, fill = drug), alpha = 0.3) +
      # themes and text
      guides(color = "none", alpha = "none") +
      scale_x_continuous(name = "#Trials", limits = c(0, nrow(df)), breaks = seq(0, nrow(df), by = 10)) +
      scale_y_continuous(name = expression(paste("Intensity  (", Delta, "BPM)")), limits = c(-35,55), breaks = seq(-35,55,10))+
      theme_minimal()+
      theme(legend.position = c(0.6,0.8),
            legend.box = "horizontal")
    
    # the second part:
    plot2_combined = individual_means %>% ggplot() +
      geom_line(aes(x = Alpha, y = prob, group = interaction(status,draw), col = status),
                linewidth = 0.5, alpha = 0.25)+coord_flip()+
      scale_x_continuous(name = expression(paste("Intensity  (", Delta, "BPM)")), limits = c(-35,55), breaks = seq(-35,55,10))+
      scale_y_continuous(name = "P(response = Faster)", limits = c(0,1), breaks = seq(0,1,0.25))+theme_minimal()+
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            legend.position = "none") +  # Remove the legend
      guides(col = "none")
    
    # combine it
    plot_draws_all = plot_combined+plot2_combined+ plot_layout(widths = c(2, 1))
    
    
    
    ############### Placebo only plot for Manuscript:
    
    fontsize = 28
    ## Only placebo first part
    plot1_placebo = rbind(v3) %>% 
      mutate(Response = ifelse(Decision == "More","Faster","Slower")) %>% 
      ggplot(aes()) +
      geom_point(aes(x = trials, y = Alpha, shape = Response), size = 3) +
      scale_shape_manual(values = c(15, 17)) +
      geom_line(aes(x = x, y = EstimatedThreshold), size = 1.1) +
      geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
      geom_ribbon(aes(x = x, ymin = lower, ymax = upper, fill = drug), alpha = 0.3) +
      guides(color = "none", alpha = "none", fill = "none",
             shape = guide_legend(direction = "horizontal")) +
      scale_x_continuous(name = "#Trials", limits = c(0, nrow(df)), breaks = seq(0, nrow(df), by = 20)) +
      scale_y_continuous(name = expression(paste("Intensity  (", Delta, "BPM)")), limits = c(-40,50), breaks = seq(-40,50,20))+
      theme_classic()+
      theme(legend.position = c(0.6,0.8),
            text = element_text(size = fontsize),         
            axis.title = element_text(size = fontsize),   
            axis.text = element_text(size = fontsize),    
            legend.text = element_text(size = fontsize),  
            legend.title = element_text(size = fontsize)  
      )
    
    dd = v3 %>% mutate(Decision_num = ifelse(Decision == "More",1,0),
                       Decision = ifelse(Decision == "More","Faster","Slower"))
    # second part only placebo
    plot2_placebo = individual_means %>% filter(status == "placebo") %>% ggplot() +
      geom_line(aes(x = Alpha, y = prob, group = interaction(status,draw), col = status),
                linewidth = 0.75, alpha = 0.25)+coord_flip()+
      scale_x_continuous(name = expression(paste("Intensity  (", Delta, "BPM)")), limits = c(-35,55), breaks = seq(-35,55,10))+
      scale_y_continuous(name = "P(Faster)", limits = c(0,1), breaks = seq(0,1,0.5))+
      theme_classic()+
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            legend.position = "none") +
      guides(col = "none")+
      theme(text = element_text(size = fontsize),         
            axis.title = element_text(size = fontsize),   
            axis.text = element_text(size = fontsize),    
            legend.text = element_text(size = fontsize), 
            legend.title = element_text(size = fontsize)  
      )
    #combine
    plot_draws_placebo = plot1_placebo+plot2_placebo+ plot_layout(widths = c(4, 1))
    #save
    ggsave(here::here("Figures","hrd_fig1.tiff"),plot_draws_placebo, height = 7, width = 12, units = "in", dpi = 400)
    #return
    return(list(plot_draws_placebo,plot_draws_all))
  }
  
  ############# RRST: plot 1
  plot1_rrst = function(){
    
    # only for placebo read data first PSI
    library(readr)
    RRST_threshold <- read_csv(here::here("data","sub_3025","visit_0003","InteroceptionTasks","RRST","RRST_threshold.csv"))
    thresholds = as.numeric(colnames(RRST_threshold))
    library(readr)
    RRST_thresholdse <- read_csv(here::here("data","sub_3025","visit_0003","InteroceptionTasks","RRST","RRST_thresholdse.csv"))
    thresholdse = as.numeric(colnames(RRST_thresholdse))
    
    # then the raw data
    df = read_csv(here::here("data","RRST_trial_level_data.csv")) %>% filter(subject == "sub_3025" & visit == "visit_0003")
    
    #renaming for the x-axis:
    stimms = c(0,8.5,17)
    labelss = round((stimms/ 17) * 100,0)
    
    #color for the plot
    color = "#0B96BC"
    
    # first plot
    placebo_plot1 = df %>% mutate(mean = thresholds, se = thresholdse, Response = ifelse(Resp == 1, "Correct","Incorrect")) %>% 
      mutate(trial = 1:n()) %>% 
      ggplot(aes(fill = drugs))+
      geom_point(aes(x = trial, y = Stim, shape = Response), size = 3) +
      # shapes' shape
      scale_shape_manual(values = c(16, 4)) +
      # 0 line
      # geom_hline(yintercept = 0, linetype = "dashed") +
      # line that is inside Confidence interval
      geom_line(aes(x = trial, y = mean),size = 1.1) +
      # confidence interval
      geom_ribbon(aes(x = trial, y = mean, ymin = mean-2*se, ymax = mean+2*se), alpha = 0.3)+
      # themes and text
      guides(color = "none", alpha = "none", fill = "none",
             shape = guide_legend(direction = "horizontal")) +
      scale_x_continuous(name = "#Trials", limits = c(0, nrow(df)), breaks = seq(0, nrow(df), by = 20)) +
      scale_y_continuous(name = expression(paste("Intensity (% resistance)")), limits = c(-2.5,22.5), breaks = stimms, labels = labelss)+
      theme_classic()+
      theme(legend.position = c(0.6,0.8),
            text = element_text(size = 28),         # Global text size
            axis.title = element_text(size = 28),   # Axis titles
            axis.text = element_text(size = 28),    # Axis tick labels
            legend.text = element_text(size = 28),  # Legend text
            legend.title = element_text(size = 28)  # Legend title (if applicable)
      )+
      scale_fill_manual(values = color)
    
    #get the subjects fitted psychometric:
    model <- readRDS(here::here("STAN models","RRST.RDS"))
    draw_id = sample(1:4000,50)
    individual_means = get_indi_rrst(model,draw_id)
    
    # plot 2
    placebo_plot2 = individual_means %>% filter(status == "placebo") %>% ggplot() +
      geom_line(aes(x = Alpha, y = prob, group = interaction(status,draw), col = status),
                linewidth = 0.75, alpha = 0.25)+coord_flip()+
      scale_x_continuous(name = expression(paste("Intensity")), limits = c(-2.5,22.5), breaks = seq(-2.5,22.5,10))+
      scale_y_continuous(name = "P(Correct)", limits = c(0.3,1), breaks = c(0.5,1))+
      theme_classic()+
      scale_color_manual(values = color)+
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            text = element_text(size = 28),         # Global text size
            axis.title = element_text(size = 28),   # Axis titles
            axis.text = element_text(size = 28),    # Axis tick labels
            legend.text = element_text(size = 28),  # Legend text
            legend.title = element_text(size = 28),  # Legend title (if applicable)
            legend.position = "none") +  # Remove the legend
      guides(col = "none")
    
    #combine
    placebo = placebo_plot1+placebo_plot2+ plot_layout(widths = c(4, 1))
    #save
    ggsave(here::here("Figures","rrst_fig1.tiff"),placebo, height = 7, width = 12, units = "in", dpi = 400)
    #return
    return(list(placebo))
    
  }
  
  hrd = plot1_hrd()
  rrst = plot1_rrst()
  
  #return the two plots (Only placebo for manuscript)
  return(list(hrd[[1]], rrst[[1]]))
  
}


# plot 2 for the manuscript
plot_2 = function(){
  
  ########## HRD:
  make_hrd = function(){
    #load model
    
    # model <- readRDS(here::here("STAN models","revisions","HRD_doublecontrol.RDS"))
    model <- readRDS(here::here("STAN models","revisions","revisions_v2","HRD_final.RDS"))
    line_width = 1.3
    fontsize = 30
    #psychometric equation HRD
    psycho_fit = function(x, alpha,lapse,beta){
      
      return(lapse + (1 - 2 * lapse) * (0.5 + 0.5 *pracma::erf((x-alpha)/(beta*sqrt(2)))))
    }
    
    box_color = "#F2F2F2"
    
    # data for the histogram on the plot (Slopes and threshold)
    
    hist = as_draws_df(model$draws(paste0("gm[",1:13,"]"))) %>% 
      mutate(beta_placebo = exp(`gm[6]`),
             beta_biso = exp(`gm[6]`+`gm[7]`),
             beta_prop = exp(`gm[6]`+`gm[8]`),
             control_beta = `gm[8]`,
             lapse_placebo = brms::inv_logit_scaled(`gm[11]`) / 2,
             lapse_biso = brms::inv_logit_scaled(`gm[11]`+`gm[12]`) / 2,
             lapse_prop = brms::inv_logit_scaled(`gm[11]`+`gm[13]`) / 2,
             alpha_placebo = `gm[1]`,
             alpha_biso = `gm[1]`+`gm[2]`,
             alpha_prop = `gm[1]`+`gm[3]`) %>% 
      mutate(draw = 1:10000)  %>% 
      select(-contains("."))%>%
      select(-contains("gm[")) %>% 
      pivot_longer(cols = starts_with(c("alpha", "beta", "lapse")), 
                   names_to = c(".value", "status"), 
                   names_sep = "_") %>% 
      pivot_longer(cols = c("alpha","beta","lapse"), values_to = "value",names_to = "parameters") %>% 
      mutate(status = ifelse(status == "placebo","placebo",ifelse(status == "biso","bisoprolol","propranolol"))) %>% 
      mutate(status = factor(status, levels = c("placebo","propranolol","bisoprolol")))
    
    hdi = hist %>% group_by(status,parameters) %>% dplyr::summarize(mean_alpha = mean(value),
                                                                    qhigh = quantile2(value, probs = c(0.05,0.95))[[1]],
                                                                    qlow = quantile2(value, probs = c(0.05,0.95))[[2]]
    )
    
    histo_data = hist %>% filter(parameters == "alpha")
    slopes = hist %>% filter(parameters == "beta")
    
    colors = c("#1200A8","#B00089","#009F73")
    
    # colors = c("#00A068","#EF9c00","#7E009E")
    
    
    # Slope histogram
    
    slopes_hist = slopes %>% ggplot(aes(x = value, fill = status))+
      # geom_histogram(aes(y = after_stat(density) / 4),col = "black", position = "identity", alpha = 0.75, bins = 25)+
      # geom_density(aes(y = after_stat(density) / 4),col = "black", position = "identity", alpha = 0.75)+
      geom_histogram(aes(y = after_stat(density) / 4), position = "identity", alpha = 0.75, bins = 25)+
      theme_void()+
      scale_fill_manual(values = colors)+
      ylab("Slope")+
      scale_color_manual(values = colors)+
      theme(
        legend.position = "none",
        axis.line.y = element_line(),          # Keep the x-axis line
        axis.ticks.y = element_line(),         # Keep x-axis ticks
        axis.text.y = element_text(size = fontsize/1.5, margin = margin(r = 3)),  # Rotate x-axis labels
        axis.title.x = element_text(size = fontsize/1.5, margin = margin(r = 3)),  # Rotate x-axis labels
        axis.line.x = element_blank(),         # Remove y-axis line
        axis.ticks.x = element_blank(),        # Remove y-axis ticks
        axis.text.x = element_blank()
      ) +
      coord_flip()+
      scale_x_continuous(labels = c(4,10), breaks = c(4,10))+
      theme(axis.line.y = element_line(linewidth = line_width))+
      scale_y_continuous(expand = c(0, 0))
    
    
    # get psychometric functions for the group means
    population_means = as_draws_df(model$draws(paste0("gm[",1:13,"]"))) %>% 
      mutate(beta_placebo = exp(`gm[6]`),
             beta_biso = exp(`gm[6]`+`gm[7]`),
             beta_prop = exp(`gm[6]`+`gm[8]`),
             control_beta = `gm[8]`,
             lapse_placebo = brms::inv_logit_scaled(`gm[11]`) / 2,
             lapse_biso = brms::inv_logit_scaled(`gm[11]`+`gm[12]`) / 2,
             lapse_prop = brms::inv_logit_scaled(`gm[11]`+`gm[13]`) / 2,
             alpha_placebo = `gm[1]`,
             alpha_biso = `gm[1]`+`gm[2]`,
             alpha_prop = `gm[1]`+`gm[3]`) %>% 
      mutate(draw = 1:10000)  %>% 
      select(-contains("."))%>%
      select(-contains("gm[")) %>% 
      pivot_longer(cols = starts_with(c("alpha", "beta", "lapse")), 
                   names_to = c(".value", "status"), 
                   names_sep = "_") %>% 
      drop_na() %>% 
      rowwise() %>%
      mutate(prob = list(psycho_fit(seq(-30,30,by = 0.1),alpha, lapse,beta)),
             Alpha = list(seq(-30,30,by = 0.1))) %>% 
      unnest() %>% 
      mutate(ids = 1) %>% 
      group_by(Alpha,status) %>%
      summarize(mean_prob = mean(prob),
                hdi_lower = HDInterval::hdi(prob, credMass = 0.95)[1],
                hdi_upper = HDInterval::hdi(prob, credMass = 0.95)[2],
                hdi_20 = HDInterval::hdi(prob, credMass = 0.80)[1],  # Example for 60% HDI
                hdi_80 = HDInterval::hdi(prob, credMass = 0.80)[2]
      ) %>% group_by(status) %>%
      mutate(
        smooth_prob = smooth.spline(Alpha, mean_prob)$y,
        smooth_hdi_lower = smooth.spline(Alpha, hdi_lower)$y,
        smooth_hdi_upper = smooth.spline(Alpha, hdi_upper)$y,
        smooth_hdi_20 = smooth.spline(Alpha, hdi_20)$y,
        smooth_hdi_80 = smooth.spline(Alpha, hdi_80)$y
      )
    
    
    ##### HRD psychometric with threshold hisogram plot:
    
    q = population_means %>% 
      mutate(status = ifelse(status == "placebo","placebo",ifelse(status == "biso","bisoprolol","propranolol"))) %>% 
      mutate(status = factor(status, levels = c("placebo","propranolol","bisoprolol"))) %>% 
      ggplot()+
      # geom_ribbon(aes(x = Alpha, y = mean_prob, ymin = hdi_lower, ymax = hdi_upper,fill = status),
      #             alpha = 0.3)+
      geom_line(aes(x = Alpha, y = mean_prob,group = status, col = status),
                # col = "black",
                linewidth = 1.2, alpha = 1)+
      geom_ribbon(aes(x = Alpha, y = smooth_prob, ymin = smooth_hdi_20, ymax = smooth_hdi_80,fill = status),
                  alpha = 0.6)+
      # geom_segment(data = histo_data %>% group_by(status) %>% mutate(median = median(value)), aes(x = median, xend = median, yend = 0.5, y = 0), linetype = 2)+
      # geom_segment(data = histo_data %>% group_by(status) %>% mutate(median = median(value)), aes(x = median, xend = 10, y = 0.5, yend = 0.5), linetype = 2)+
      
      # geom_histogram(data = histo_data, aes(x = value,y = after_stat(density) / 2, fill = status),col = "black", position = "identity", alpha = 0.75, bins = 125)+
      # geom_density(data = histo_data, aes(x = value,y = after_stat(density) / 2, fill = status),col = "black", position = "identity", alpha = 0.75)+
      geom_histogram(data = histo_data, aes(x = value,y = after_stat(density) / 2, fill = status), position = "identity", alpha = 0.75, bins = 125)+
      theme_classic()+
      
      # scale_fill_brewer(palette = "Dark2")+
      # scale_color_brewer(palette = "Dark2")+
      scale_fill_manual(values = colors)+
      scale_color_manual(values = colors)+
      labs(fill = "Drug", color = "Drug")+
      xlab("Stimulus intensity (ΔBPM)")+
      ylab("P(faster | ΔBPM)")+
      scale_x_continuous(limits = c(-31,31.5),breaks = seq(-30,30,by = 15), labels = seq(-30,30,by = 15),expand = c(0, 0))+
      scale_y_continuous(limits = c(0,1.05), expand = c(0, 0), breaks = c(0,0.50,1.00), labels = c("0.00","0.50","1.00")) +  # Ensures no padding below x-axis
      theme(legend.position = c(0.15,0.8),
            text = element_text(size = fontsize),           # All text
            axis.title = element_text(size = fontsize),     # Axis titles
            axis.text.x = element_text(size = fontsize),      # Axis tick labels
            axis.text.y = element_text(size = fontsize),      # Axis tick labels
            legend.text = element_text(size = fontsize),    # Legend text
            legend.title = element_text(size = fontsize),   # Legend title
            plot.title = element_text(size = fontsize),
            axis.line.y = element_line(linewidth = line_width),
            axis.line.x = element_line(linewidth = line_width)
      )
    
    # add the Slopes histogram
    slopes_hist <- ggplotGrob(slopes_hist)
    qq = q +
      annotation_custom(
        grob = slopes_hist,
        xmin = 4, xmax = 15,  # Adjust position
        ymin = 0.35, ymax = 0.65   # Adjust position
      )
    
    #save it
    ggsave(here::here("Figures","plot2.tiff"),qq, height = 7, width = 12, units = "in", dpi = 400)
    
    # include marginal differences below!
    
    marginal_dif = as_draws_df(model$draws(paste0("gm[",1:13,"]"))) %>% 
      mutate(slope_biso_v_placebo = `gm[7]`,
             slope_prop_v_placebo = `gm[8]`,
             lapse_biso_v_placebo = `gm[12]`,
             lapse_prop_v_placebo = `gm[13]`,
             alpha_biso_v_placebo = `gm[2]`,
             alpha_prop_v_placebo = `gm[3]`) %>% 
      mutate(draw = 1:10000)  %>% 
      select(-contains("."))%>%
      select(-contains("gm[")) %>% 
      pivot_longer(-draw) %>% 
      mutate(paramter = sub("_.*", "", name),
             name = sub("^[^_]*_", "", name)        # Remove prefix, keeping condition
      ) %>% # Extract the prefix before the first underscore
      pivot_wider(names_from = paramter, values_from = value) %>% 
      pivot_longer(cols = c(slope,lapse,alpha),names_to = "Parameters") %>% 
      rename(Contrast = name)
    
    hdi = marginal_dif %>% group_by(Contrast,Parameters) %>% dplyr::summarize(mean = mean(value),
                                                                              q_95h = quantile2(value, probs = c(0.05,0.95))[[1]],
                                                                              q_5l = quantile2(value, probs = c(0.05,0.95))[[2]],
                                                                              q_80h = quantile2(value, probs = c(0.20,0.80))[[1]],
                                                                              q_20l = quantile2(value, probs = c(0.20,0.80))[[2]]
    )
    
    colors = c("#009F73","white","#B00089")
    
    #Marginal difference for the threshold
    
    alpha = hdi %>% 
      mutate(Contrast = ifelse(Contrast == "biso_v_placebo", "Bisoprolol - Placebo", ifelse(Contrast == "prop_v_placebo", "Propanalol - Placebo"))) %>% 
      filter(Parameters == "alpha") %>% ungroup() %>% 
      mutate(Parameters = "Threshold") %>% add_row(Contrast = "", Parameters = "Threshold", mean = NA, q_95h = NA, q_5l = NA, q_80h = NA, q_20l = NA) %>% 
      mutate(Contrast = factor(Contrast, levels = c("Bisoprolol - Placebo", "", "Propanalol - Placebo"))) %>% 
      mutate(y_val = as.numeric(as.factor(Contrast))) %>% 
      ggplot(aes(col = Contrast))+
      geom_pointrange(aes(y = y_val, x = mean, xmin = q_5l, xmax = q_95h), linewidth = 1.5, show.legend = FALSE)+
      geom_pointrange(aes(y = y_val, x = mean, xmin = q_80h, xmax = q_20l),linewidth = 3, show.legend = FALSE)+
      geom_point(aes(y = y_val, x = mean),size = 5, show.legend = FALSE)+
      # geom_vline(xintercept = 0, linetype = 2, linewidth = line_width, col = "#1200A8", show.legend = FALSE)+
      facet_wrap(~Parameters, scales = "free")+
      theme_classic()+
      xlab("")+ylab("")+
      scale_x_continuous(breaks = c(0,2.5,5.5), labels = c("0","2.5","5"))+
      coord_cartesian(xlim = c(0,5), ylim = c(-1,5))+
      geom_segment(aes(x = 0, xend = 0, y = -2.3,yend = 8.3), linetype = 2, linewidth = line_width, col = "#1200A8", show.legend = FALSE)+
      theme(legend.position = "none",
            axis.text.y = element_blank(),
            axis.line.y = element_blank(),        # Remove the y-axis line
            strip.background = element_rect(fill = box_color, color = "white"),  # Default box style
            axis.ticks.y = element_blank())+scale_color_manual(values = colors)+  # Ensures no padding below x-axis
      theme(text = element_text(size = fontsize),           # All text
            axis.title = element_text(size = fontsize),     # Axis titles
            axis.text.x = element_text(size = fontsize),      # Axis tick labels
            legend.text = element_text(size = fontsize),    # Legend text
            legend.title = element_text(size = fontsize),   # Legend title
            plot.title = element_text(size = fontsize),      # Plot title
            axis.line.x = element_line(linewidth = line_width) 
      )
    
    #for the slope
    
    beta = hdi %>% 
      filter(Parameters == "slope") %>% 
      mutate(Parameters = "Slope") %>% 
      ungroup() %>% 
      mutate(Parameters = "Slope") %>% 
      add_row(Contrast = "", Parameters = "Slope", mean = NA, q_95h = NA, q_5l = NA, q_80h = NA, q_20l = NA) %>% 
      mutate(Contrast = factor(Contrast, levels = c("biso_v_placebo", "", "prop_v_placebo"))) %>% 
      mutate(y_val = as.numeric(as.factor(Contrast))) %>% 
      ggplot(aes(col = Contrast))+
      geom_pointrange(aes(y = y_val, x = mean, xmin = q_5l, xmax = q_95h), linewidth = 1.5, show.legend = FALSE)+
      geom_pointrange(aes(y = y_val, x = mean, xmin = q_80h, xmax = q_20l),linewidth = 3, show.legend = FALSE)+
      geom_point(aes(y = y_val, x = mean),size = 5, show.legend = FALSE)+
      # geom_vline(xintercept = 0, linetype = 2, linewidth = line_width, col = "#1200A8", show.legend = FALSE)+
      coord_cartesian(xlim = c(-0.25,0.25), ylim = c(-1,5))+
      geom_segment(aes(x = 0, xend = 0, y = -2.3,yend = 8.3), linetype = 2, linewidth = line_width, col = "#1200A8", show.legend = FALSE)+
      facet_wrap(~Parameters, scales = "free")+
      theme_classic()+ theme(
        axis.title.y = element_blank(),       # Remove y-axis title
        axis.text.y = element_blank(),        # Remove y-axis text labels
        strip.background = element_rect(fill = box_color, color = "white"),  # Default box style
        axis.ticks.y = element_blank(),       # Remove y-axis ticks
        axis.line.y = element_blank(),        # Remove the y-axis line
        panel.spacing = unit(1, "lines")      # Optional: Adjust space between facets
      )+
      xlab("Differences in parameter estimates")+
      xlab("")+
      scale_x_continuous(breaks = c(-0.2,0,0.2), labels = c("-0.2","0","0.2"))+
      theme(legend.position = "none")+scale_color_manual(values = colors)+
      theme(
        text = element_text(size = fontsize),           # All text
        axis.title = element_text(size = fontsize),     # Axis titles
        axis.text.x = element_text(size = fontsize),      # Axis tick labels
        legend.text = element_text(size = fontsize),    # Legend text
        legend.title = element_text(size = fontsize),   # Legend title
        plot.title = element_text(size = fontsize),      # Plot title
        axis.line.x = element_line(linewidth = line_width) 
      )
    
    
    # for the lapse rate:
    
    lapse = hdi %>% 
      filter(Parameters == "lapse") %>% 
      mutate(Parameters = "Lapse") %>% 
      ungroup() %>% 
      mutate(Parameters = "Lapse") %>% 
      add_row(Contrast = "", Parameters = "Lapse", mean = NA, q_95h = NA, q_5l = NA, q_80h = NA, q_20l = NA) %>% 
      mutate(Contrast = factor(Contrast, levels = c("biso_v_placebo", "", "prop_v_placebo"))) %>% 
      mutate(y_val = as.numeric(as.factor(Contrast))) %>% 
      ggplot(aes(col = Contrast))+
      geom_pointrange(aes(y = y_val, x = mean, xmin = q_5l, xmax = q_95h), linewidth = 1.5, show.legend = FALSE)+
      geom_pointrange(aes(y = y_val, x = mean, xmin = q_80h, xmax = q_20l),linewidth = 3, show.legend = FALSE)+
      geom_point(aes(y = y_val, x = mean),size = 5, show.legend = FALSE)+
      # geom_vline(xintercept = 0, linetype = 2, linewidth = line_width, col = "#1200A8", show.legend = FALSE)+
      coord_cartesian(xlim = c(-4.5,4.5), ylim = c(-1,5))+
      geom_segment(aes(x = 0, xend = 0, y = -2.3,yend = 8.3), linetype = 2, linewidth = line_width, col = "#1200A8", show.legend = FALSE)+
      facet_wrap(~Parameters, scales = "free")+
      theme_classic()  + theme(
        axis.title.y = element_blank(),       # Remove y-axis title
        axis.text.y = element_blank(),        # Remove y-axis text labels
        strip.background = element_rect(fill = box_color, color = "white"),  # Default box style
        axis.ticks.y = element_blank(),       # Remove y-axis ticks
        axis.line.y = element_blank(),        # Remove the y-axis line
        panel.spacing = unit(1, "lines")      # Optional: Adjust space between facets
      )+
      xlab("")+
      theme(legend.position = "none")+scale_color_manual(values = colors)+
      scale_x_continuous(breaks = c(-4,0,4), labels = c("-4","0","4"))+
      theme(
        text = element_text(size = fontsize),           # All text
        axis.title = element_text(size = fontsize),     # Axis titles
        axis.text.x = element_text(size = fontsize),      # Axis tick labels
        legend.text = element_text(size = fontsize),    # Legend text
        legend.title = element_text(size = fontsize),   # Legend title
        plot.title = element_text(size = fontsize),      # Plot title
        axis.line.x = element_line(linewidth = line_width) 
      )
    #return it all:
    return(list(alpha,beta,lapse,qq))
  }
  
  
  ########## RRST:
  make_rrst = function(){
    box_color = "#F2F2F2"
    
    #rrst psychometric equation
    psycho_fit_rrst_05 = function(x, alpha,lapse,beta){
      return(0.5 + (1 - 0.5 - lapse) * (1-exp(-10^(beta * (x-alpha)))))
    }
    
    # plot metrics and load model
    line_width = 1.3
    fontsize = 30
    colors = c("#1200A8","#B00089","#009F73")
    # model_rrst <- readRDS(here::here("STAN models","revisions","RRST_doublecontrol.RDS"))
    model_rrst <- readRDS(here::here("STAN models","revisions","revisions_v2","RRST_final.RDS"))
    # histograms slopes and threshold data
    hist_rrst = as_draws_df(model_rrst$draws(paste0("gm[",1:13,"]"))) %>% 
      mutate(beta_placebo = exp(`gm[1]`),
             beta_biso = exp(`gm[1]`+`gm[2]`),
             beta_prop = exp(`gm[1]`+`gm[3]`),
             lapse_placebo = brms::inv_logit_scaled(`gm[6]`) / 2,
             lapse_biso = brms::inv_logit_scaled(`gm[6]`+`gm[7]`) / 2,
             lapse_prop = brms::inv_logit_scaled(`gm[6]`+`gm[8]`) / 2,
             alpha_placebo = `gm[9]`,
             alpha_biso = `gm[9]`+`gm[10]`,
             alpha_prop = `gm[9]`+`gm[11]`) %>% 
      mutate(draw = 1:10000)  %>% 
      select(-contains("."))%>%
      select(-contains("gm[")) %>% 
      pivot_longer(cols = starts_with(c("alpha", "beta", "lapse")), 
                   names_to = c(".value", "status"), 
                   names_sep = "_")%>% 
      pivot_longer(cols = c("alpha","beta","lapse"), values_to = "value",names_to = "parameters") %>% 
      mutate(status = ifelse(status == "placebo","placebo",ifelse(status == "biso","bisoprolol","propranolol"))) %>% 
      mutate(status = factor(status, levels = c("placebo","propranolol","bisoprolol")))  
    
    hdi = hist_rrst %>% group_by(status,parameters) %>% dplyr::summarize(mean_alpha = mean(value),
                                                                         qhigh = quantile2(value, probs = c(0.05,0.95))[[1]],
                                                                         qlow = quantile2(value, probs = c(0.05,0.95))[[2]]
    )
    
    #data for threshold and slopes
    histo_data_rrst = hist_rrst %>% filter(parameters == "alpha")
    slopes_rrst = hist_rrst %>% filter(parameters == "beta")
    
    
    #Histogram for the slopes:
    
    slopes_hist_rrst = slopes_rrst %>% ggplot(aes(x = value, fill = status))+
      # geom_histogram(aes(y = after_stat(density) / 4),col = "black", position = "identity", alpha = 0.75, bins = 25)+
      # geom_density(aes(y = after_stat(density) / 4),col = "black", position = "identity", alpha = 0.75)+
      geom_histogram(aes(y = after_stat(density) / 4), position = "identity", alpha = 0.75, bins = 25)+
      theme_void()+
      ylab("Slope")+
      scale_fill_manual(values = colors)+
      scale_color_manual(values = colors)+
      theme(
        legend.position = "none",
        axis.line.y = element_line(),          # Keep the x-axis line
        axis.ticks.y = element_line(),         # Keep x-axis ticks
        axis.text.y = element_text(size = fontsize/1.5, margin = margin(r = 3)),  # Rotate x-axis labels
        axis.title.x = element_text(size = fontsize/1.5, margin = margin(r = 3)),  # Rotate x-axis labels
        axis.line.x = element_blank(),         # Remove y-axis line
        axis.ticks.x = element_blank(),        # Remove y-axis ticks
        axis.text.x = element_blank()          # Remove y-axis text
      ) +
      coord_flip()+
      scale_x_continuous(labels = c(0.2,0.4), breaks = c(0.2,0.4))+
      theme(axis.line.y = element_line(linewidth = line_width))+
      scale_y_continuous(expand = c(0, 0))
    
    
    
    # group means for the first plot of psychometric functions
    population_means_rrst = as_draws_df(model_rrst$draws(paste0("gm[",1:13,"]"))) %>% 
      mutate(beta_placebo = exp(`gm[1]`),
             beta_biso = exp(`gm[1]`+`gm[2]`),
             beta_prop = exp(`gm[1]`+`gm[3]`),
             lapse_placebo = brms::inv_logit_scaled(`gm[6]`) / 2,
             lapse_biso = brms::inv_logit_scaled(`gm[6]`+`gm[7]`) / 2,
             lapse_prop = brms::inv_logit_scaled(`gm[6]`+`gm[8]`) / 2,
             alpha_placebo = `gm[9]`,
             alpha_biso = `gm[9]`+`gm[10]`,
             alpha_prop = `gm[9]`+`gm[11]`) %>% 
      mutate(draw = 1:10000)  %>% 
      # filter(draw %in% draw_iq) %>% 
      select(-contains("."))%>%
      select(-contains("gm[")) %>% 
      pivot_longer(cols = starts_with(c("alpha", "beta", "lapse")), 
                   names_to = c(".value", "status"), 
                   names_sep = "_") %>% 
      drop_na() %>% 
      rowwise() %>%
      mutate(prob = list(psycho_fit_rrst_05(seq(0,20,by = 0.1),alpha, lapse,beta)),
             Alpha = list(seq(0,20,by = 0.1))) %>% 
      unnest() %>% 
      mutate(ids = 1) %>% 
      group_by(Alpha,status) %>%
      summarize(mean_prob = mean(prob),
                hdi_lower = quantile2(prob)[[1]],
                hdi_upper = quantile2(prob)[[2]],
                hdi_20 = quantile2(prob, probs = c(0.2,0.8))[[1]],
                hdi_80 = quantile2(prob, probs = c(0.2,0.8))[[2]],
      )%>% group_by(status) %>%
      mutate(
        smooth_prob = smooth.spline(Alpha, mean_prob)$y,
        smooth_hdi_lower = smooth.spline(Alpha, hdi_lower)$y,
        smooth_hdi_upper = smooth.spline(Alpha, hdi_upper)$y,
        smooth_hdi_20 = smooth.spline(Alpha, hdi_20)$y,
        smooth_hdi_80 = smooth.spline(Alpha, hdi_80)$y
      )
    
    #recoding for the x-axis
    stimms = c(0,4.25,8.5,12.75,17)
    labelss = round((stimms/ 17) * 100,0)
    
    
    
    
    population_means_rrst$status = factor(population_means_rrst$status, levels = c("placebo", "biso", "prop"))
    
    #first plot of group level psychometric functions
    q_rrst = population_means_rrst %>% 
      mutate(status = ifelse(status == "placebo","placebo",ifelse(status == "biso","bisoprolol","propranolol"))) %>% 
      mutate(status = factor(status, levels = c("placebo","propranolol","bisoprolol"))) %>% 
      ggplot()+
      # geom_ribbon(aes(x = Alpha, y = mean_prob, ymin = hdi_lower, ymax = hdi_upper,fill = status),
      #             alpha = 0.3)+
      geom_line(aes(x = Alpha, y = mean_prob,group = status, col = status),
                # col = "black",
                linewidth = 1, alpha = 1, show.legend = FALSE)+
      geom_ribbon(aes(x = Alpha, y = smooth_prob, ymin = smooth_hdi_20, ymax = smooth_hdi_80,fill = status),
                  alpha = 0.6, show.legend = FALSE)+
      # geom_segment(data = histo_data %>% group_by(status) %>% mutate(median = median(value)), aes(x = median, xend = median, yend = 0.5, y = 0), linetype = 2)+
      # geom_segment(data = histo_data %>% group_by(status) %>% mutate(median = median(value)), aes(x = median, xend = 10, y = 0.5, yend = 0.5), linetype = 2)+
      
      # geom_histogram(data = histo_data_rrst, aes(x = value,y = (after_stat(density) / 12)+0.5, fill = status),col = "black", position = "identity", alpha = 0.75, bins = 125)+
      # geom_density(data = histo_data_rrst, aes(x = value,y = (after_stat(density) / 12)+0.5, fill = status),col = "black", position = "identity", alpha = 0.75)+
      geom_histogram(data = histo_data_rrst, aes(x = value,y = (after_stat(density) / 12)+0.5, fill = status), position = "identity", alpha = 0.75, bins = 125, show.legend = FALSE)+
      #geom_pointrange(data = raw_data, aes(x = Alpha, y = mean, ymin = mean-se_resp, ymax = mean+se_resp, col = status))+
      theme_classic()+
      scale_fill_manual(values = colors)+
      scale_color_manual(values = colors)+
      labs(fill = "Drug", color = "Drug")+
      xlab("Stimulus intensity (% RRes)")+
      ylab("P(correct | % RRes)")+
      # theme(legend.position = "none")+
      scale_x_continuous(breaks = seq(0,20.2,by = 5), labels = labelss,expand = c(0, 0))+
      coord_cartesian(ylim = c(0.5,1.05), xlim = c(4.5,20.5))+
      scale_y_continuous(expand = c(0, 0), breaks = c(0.50,0.75,1.00), labels = c("0.50","0.75","1.00")) +  # Ensures no padding below x-axis
      theme(text = element_text(size = fontsize),           # All text
            axis.title = element_text(size = fontsize),     # Axis titles
            axis.text.x = element_text(size = fontsize),      # Axis tick labels
            axis.text.y = element_text(size = fontsize),      # Axis tick labels
            legend.text = element_text(size = fontsize),    # Legend text
            legend.title = element_text(size = fontsize),   # Legend title
            plot.title = element_text(size = fontsize),
            axis.line.y = element_line(linewidth = line_width),
            axis.line.x = element_line(linewidth = line_width)
      )
    
    # add the Slopes histogram
    slopes_hist_rrst <- ggplotGrob(slopes_hist_rrst)
    
    qq_rrst = q_rrst +
      annotation_custom(
        grob = slopes_hist_rrst,
        xmin = 17, xmax = 20,  # Adjust position
        ymin = 0.7, ymax = 0.9   # Adjust position
      )
    
    #Marginal difference (data)
    
    marginal_dif = as_draws_df(model_rrst$draws(paste0("gm[",1:13,"]"))) %>% 
      mutate(slope_biso_v_placebo = `gm[2]`,
             slope_prop_v_placebo = `gm[3]`,
             lapse_biso_v_placebo = `gm[7]`,
             lapse_prop_v_placebo = `gm[8]`,
             alpha_biso_v_placebo = `gm[10]`,
             alpha_prop_v_placebo = `gm[11]`) %>% 
      mutate(draw = 1:10000)  %>% 
      # filter(draw %in% draw_iq) %>% 
      select(-contains("."))%>%
      select(-contains("gm[")) %>% 
      pivot_longer(-draw) %>% 
      mutate(paramter = sub("_.*", "", name),
             name = sub("^[^_]*_", "", name)        # Remove prefix, keeping condition
      ) %>% # Extract the prefix before the first underscore
      pivot_wider(names_from = paramter, values_from = value) %>% 
      pivot_longer(cols = c(slope,lapse,alpha),names_to = "Parameters") %>% 
      rename(Contrast = name)
    
    hdi = marginal_dif %>% group_by(Contrast,Parameters) %>% dplyr::summarize(mean = mean(value),
                                                                              q_95h = quantile2(value, probs = c(0.05,0.95))[[1]],
                                                                              q_5l = quantile2(value, probs = c(0.05,0.95))[[2]],
                                                                              q_80h = quantile2(value, probs = c(0.20,0.80))[[1]],
                                                                              q_20l = quantile2(value, probs = c(0.20,0.80))[[2]]
    )
    
    # making the marginal difference below the plot. for threshold, slope and lapse
    
    colors = c("#009F73","white","#B00089")
    
    alpha_rrst = hdi %>% 
      mutate(Contrast = ifelse(Contrast == "biso_v_placebo", "Bisoprolol - Placebo", ifelse(Contrast == "prop_v_placebo", "Propanalol - Placebo"))) %>% 
      filter(Parameters == "alpha") %>% 
      mutate(Parameters = "Threshold") %>% 
      ungroup() %>% 
      mutate(Parameters = "Threshold") %>% add_row(Contrast = "", Parameters = "Threshold", mean = NA, q_95h = NA, q_5l = NA, q_80h = NA, q_20l = NA) %>% 
      mutate(Contrast = factor(Contrast, levels = c("Bisoprolol - Placebo", "", "Propanalol - Placebo"))) %>% 
      mutate(y_val = as.numeric(as.factor(Contrast))) %>% 
      ggplot(aes(col = Contrast, show.legend = FALSE))+
      geom_pointrange(aes(y = y_val, x = mean, xmin = q_5l, xmax = q_95h), linewidth = 1.5, show.legend = FALSE)+
      geom_pointrange(aes(y = y_val, x = mean, xmin = q_80h, xmax = q_20l),linewidth = 3, show.legend = FALSE)+
      geom_point(aes(y = y_val, x = mean),size = 5, show.legend = FALSE)+
      # geom_vline(xintercept = 0, linetype = 2, linewidth = line_width, col = "#1200A8", show.legend = FALSE)+
      coord_cartesian(xlim = c(-0.9,0.45), ylim = c(-1,5))+
      geom_segment(aes(x = 0, xend = 0, y = -2.3,yend = 8.3), linetype = 2, linewidth = line_width, col = "#1200A8", show.legend = FALSE)+
      facet_wrap(~Parameters, scales = "free")+
      theme_classic()+
      xlab("")+ylab("")+
      scale_x_continuous(breaks = c(-0.8,0,0.4), labels = c("-0.8","0","0.4"))+
      theme(axis.text.y = element_blank(),
            strip.background = element_rect(fill = box_color, color = "white"),  # Default box style
            axis.line.y = element_blank(),        # Remove the y-axis line
            axis.ticks.y = element_blank())+scale_color_manual(values = colors)+  # Ensures no padding below x-axis
      theme(text = element_text(size = fontsize),           # All text
            axis.title = element_text(size = fontsize),     # Axis titles
            axis.text.x = element_text(size = fontsize),      # Axis tick labels
            legend.text = element_text(size = fontsize),    # Legend text
            legend.title = element_text(size = fontsize),   # Legend title
            plot.title = element_text(size = fontsize),      # Plot title
            axis.line.x = element_line(linewidth = line_width) 
      )
    
    # for the slope
    
    beta_rrst = hdi %>% 
      filter(Parameters == "slope") %>% 
      mutate(Parameters = "Slope") %>% 
      ungroup() %>% 
      mutate(Parameters = "Slope") %>% 
      add_row(Contrast = "", Parameters = "Slope", mean = NA, q_95h = NA, q_5l = NA, q_80h = NA, q_20l = NA) %>% 
      mutate(Contrast = factor(Contrast, levels = c("biso_v_placebo", "", "prop_v_placebo"))) %>%
      mutate(y_val = as.numeric(as.factor(Contrast))) %>% 
      ggplot(aes(col = Contrast))+
      geom_pointrange(aes(y = y_val, x = mean, xmin = q_5l, xmax = q_95h), linewidth = 1.5, show.legend = FALSE)+
      geom_pointrange(aes(y = y_val, x = mean, xmin = q_80h, xmax = q_20l),linewidth = 3, show.legend = FALSE)+
      geom_point(aes(y = y_val, x = mean),size = 5, show.legend = FALSE)+
      # geom_vline(xintercept = 0, linetype = 2, linewidth = line_width, col = "#1200A8", show.legend = FALSE)+
      coord_cartesian(xlim = c(0,0.7), ylim = c(-1,5))+
      geom_segment(aes(x = 0, xend = 0, y = -2.3,yend = 8.3), linetype = 2, linewidth = line_width, col = "#1200A8", show.legend = FALSE)+
      facet_wrap(~Parameters, scales = "free")+
      theme_classic()+ theme(
        axis.title.y = element_blank(),       # Remove y-axis title
        axis.text.y = element_blank(),        # Remove y-axis text labels
        axis.ticks.y = element_blank(),       # Remove y-axis ticks
        axis.line.y = element_blank(),        # Remove the y-axis line
        panel.spacing = unit(1, "lines")      # Optional: Adjust space between facets
      )+
      #xlab("Differences in parameter estimates")+
      xlab("")+
      theme(legend.position = "none")+scale_color_manual(values = colors)+
      scale_x_continuous(breaks = c(0,0.3,0.6), labels = c("0","0.3","0.6"))+
      theme(
        text = element_text(size = fontsize),           # All text
        axis.title = element_text(size = fontsize),     # Axis titles
        axis.text.x = element_text(size = fontsize),      # Axis tick labels
        strip.background = element_rect(fill = box_color, color = "white"),  # Default box style
        legend.text = element_text(size = fontsize),    # Legend text
        legend.title = element_text(size = fontsize),   # Legend title
        plot.title = element_text(size = fontsize),      # Plot title
        axis.line.x = element_line(linewidth = line_width) 
      )
    
    # for the lapse rate
    
    lapse_rrst = hdi %>% 
      filter(Parameters == "lapse") %>% 
      mutate(Parameters = "Lapse") %>% 
      ungroup() %>% 
      mutate(Parameters = "Lapse") %>% 
      add_row(Contrast = "", Parameters = "Lapse", mean = NA, q_95h = NA, q_5l = NA, q_80h = NA, q_20l = NA) %>% 
      mutate(Contrast = factor(Contrast, levels = c("biso_v_placebo", "", "prop_v_placebo"))) %>% 
      mutate(y_val = as.numeric(as.factor(Contrast))) %>% 
      ggplot(aes(col = Contrast))+
      geom_pointrange(aes(y = y_val, x = mean, xmin = q_5l, xmax = q_95h), linewidth = 1.5, show.legend = FALSE)+
      geom_pointrange(aes(y = y_val, x = mean, xmin = q_80h, xmax = q_20l),linewidth = 3, show.legend = FALSE)+
      geom_point(aes(y = y_val, x = mean),size = 5, show.legend = FALSE)+
      # geom_vline(xintercept = 0, linetype = 2, linewidth = line_width, col = "#1200A8", show.legend = FALSE)+
      coord_cartesian(xlim = c(-2.25,1.5), ylim = c(-1,5))+
      geom_segment(aes(x = 0, xend = 0, y = -2.3,yend = 8.3), linetype = 2, linewidth = line_width, col = "#1200A8", show.legend = FALSE)+
      facet_wrap(~Parameters, scales = "free")+
      theme_classic()  + theme(
        axis.title.y = element_blank(),       # Remove y-axis title
        axis.text.y = element_blank(),        # Remove y-axis text labels
        strip.background = element_rect(fill = box_color, color = "white"),  # Default box style
        axis.ticks.y = element_blank(),       # Remove y-axis ticks
        axis.line.y = element_blank(),        # Remove the y-axis line
        panel.spacing = unit(1, "lines")      # Optional: Adjust space between facets
      )+
      xlab("")+
      scale_x_continuous(breaks = c(-2,-1,0,1), labels = c("-2","-1","0","1"))+
      theme(legend.position = "none")+
      scale_color_manual(values = colors)+
      theme(
        text = element_text(size = fontsize),           # All text
        axis.title = element_text(size = fontsize),     # Axis titles
        axis.text.x = element_text(size = fontsize),      # Axis tick labels
        legend.text = element_text(size = fontsize),    # Legend text
        legend.title = element_text(size = fontsize),   # Legend title
        plot.title = element_text(size = fontsize),      # Plot title
        axis.line.x = element_line(linewidth = line_width) 
      )
    
    
    #return it all
    
    return(list(alpha_rrst,beta_rrst,lapse_rrst,qq_rrst))
    
  }
  
  ##### Combine for full plot!
  
  hrds = make_hrd()
  rrsts = make_rrst()
  
  # combine plots:
  lower = hrds[[1]]|hrds[[2]]|hrds[[3]]|plot_spacer()|rrsts[[1]]|rrsts[[2]]|rrsts[[3]]|plot_layout(nrow = 1, widths = c(1,1,1,0.3,1,1,1))
  upper = (hrds[[4]]|rrsts[[4]])#+ plot_layout(heights = c(4,0.1, 4),widths = c(4,0.1, 4))
  
  plot = (upper/plot_spacer()/lower + plot_layout(heights = c(4,0.1, 1)))+ 
    plot_layout(guides = "collect")&
    theme(legend.position = "bottom", legend.title = element_blank())
  #watch it
  plot
  #save it
  ggsave(here::here("Figures","plot2_revision_visions2.tiff"),plot, height = 10, width = 18, units = "in", dpi = 400)
  #return it
  return(plot)
}

# plot_2()

# Plot 3:

# helper function:

clean_and_report <- function(data, dataset_name) {
  # Count the number of rows before removing NAs
  total_rows_before <- nrow(data)
  
  # Remove rows with any NAs
  cleaned_data <- data %>%
    drop_na()
  
  # Count the number of rows after removing NAs
  total_rows_after <- nrow(cleaned_data)
  
  # Calculate the percentage of dropped trials
  dropped_trials_percentage <- ((total_rows_before - total_rows_after) / total_rows_before) * 100
  
  # Print a summary report
  cat(paste0(
    "Dataset: ", dataset_name, "\n",
    "Total rows before cleaning: ", total_rows_before, "\n",
    "Total rows after cleaning: ", total_rows_after, "\n",
    "Percentage of dropped trials due to NAs: ", round(dropped_trials_percentage, 2), "%\n\n"
  ))
  
  # Return cleaned data
  return(cleaned_data)
}


## HRD plots
plot_3 = function(){
  
  
  # Get metacognition data for HRD
  HRD_trl_data = read.csv(here::here("data","cleaned","HRD.csv")) %>% 
    mutate(Conf = Confidence/100,
           drugs = as.factor(drugs),
           Drug = factor(drugs, levels = c("PLACEBO", "PROP", "BISO"), labels = c("placebo", "propranolol", "bisoprolol")),
           Drug = relevel(Drug, ref = "placebo"),
           ResponseCorrect = factor(ResponseCorrect, levels = c("True", "False"), labels = c("Correct", "Incorrect")),
           BPM_scaled = scale(listenBPM),
           Condition = drugs
    ) %>% 
    filter(subject != "sub_4049")
  
  # Cleaning HRD data
  HRD_trl_data <- clean_and_report(HRD_trl_data, "HRD_trl_data")
  
  # HRD - Ordered beta regression confidence model
  HRD_conf_model <- glmmTMB(
    Conf ~ Drug * ResponseCorrect + BPM_scaled +  (1  + ResponseCorrect + BPM_scaled | subject),
    data = HRD_trl_data,
    family = ordbeta(),
    start = list(psi = c(0, 1)),
    control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000))
  )
  
  
  # Get metacognition data for RRST
  RRST_trial_data <- read_csv(here::here("data","cleaned","RRST.csv")) %>% 
    filter(ConfResp >= 0 & ConfResp <= 100) %>%
    mutate(Conf = ConfResp/100,
           Stimulus = Stim,
           Drug = factor(drugs, levels = c("PLACEBO", "PROP", "BISO"), labels = c("placebo", "propranolol", "bisoprolol")),
           Drug = relevel(Drug, ref = "placebo"),
           Condition = drugs,
           Accuracy = factor(Resp, levels = c(1, 0), labels = c("Correct", "Incorrect"))
    )
  
  # Cleaning RRST data
  RRST_trial_data <- clean_and_report(RRST_trial_data, "RRST_trial_data")
  
  # RRST - Ordered beta regression confidence model
  RRST_conf_model <- glmmTMB(
    Conf ~ Stimulus + Drug * Accuracy + (1 + Stimulus + Accuracy | subject),
    data = RRST_trial_data,
    family = ordbeta(),
    control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000))
  )
  
  preds <- ggpredict(HRD_conf_model, terms = c("Drug", "ResponseCorrect"))
  
  ribbon_col <- "#B0B0B0"
  binary_colourmap <- c("#18699D", "#E8BB1B")
  line_width = 1.3
  fontsize = 18
  
  # Create the plot for HRD predicted confidence
  p1 <-ggplot(preds, aes(x = x, y = predicted, group = group, color = group)) +
    geom_line(size = 1.3, show.legend = FALSE) +
    geom_point(size = 4, show.legend = FALSE) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = ribbon_col), alpha = 0.15, outline.type = "both", show.legend = FALSE) +
    scale_color_manual(values = binary_colourmap) +
    scale_fill_manual(values = ribbon_col) +
    labs(
      x = "",
      y = "Predicted Confidence",
      color = "Response Correct",
      fill = ""
    ) +
    theme_classic(base_size = 12) +
    coord_cartesian(ylim = c(0.43,0.71))+
    scale_x_discrete(expand = c(0.1,0))+
    scale_y_continuous(breaks = c(0.5,0.6,0.7), labels = c("0.5","0.6","0.7"))+
    theme(plot.margin = margin(t = 6, r = 6, b = 6, l = 6), "pt")+
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = fontsize),
      legend.text = element_text(size = fontsize),
      axis.text = element_text(size = 12),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = fontsize),
      axis.text.x = element_text(size = fontsize), #, hjust = 0.6),      # Axis tick labels
      axis.text.y = element_text(size = fontsize),      # Axis tick labels
      plot.title = element_text(size = fontsize), #, face = "bold")
      axis.line.y = element_line(linewidth = 1),
      axis.line.x = element_line(linewidth = 1)
    )
  
  # Obtain estimated marginal means for each Drug and ResponseCorrect combination
  emm <- emmeans(HRD_conf_model, ~ Drug | ResponseCorrect)
  
  # Perform pairwise comparisons between Drugs within each ResponseCorrect level
  contrast_results <- contrast(emm, method = "pairwise")
  
  # Convert the contrast results to a data frame
  contrast_df <- as.data.frame(contrast_results)
  
  # Filter contrasts for ResponseCorrect == Correct
  contrast_correct <- subset(contrast_df, ResponseCorrect == "Correct")
  
  # Calculate 95% confidence intervals
  # Since df is infinite (df = Inf), we'll use the normal distribution (z-distribution)
  z_value <- qnorm(0.975)  # For 95% CI
  
  contrast_correct$lower.CL <- contrast_correct$estimate - z_value * contrast_correct$SE
  contrast_correct$upper.CL <- contrast_correct$estimate + z_value * contrast_correct$SE
  
  
  # Create the HRD contrasts plot
  colors = c("#404040","#009F73","#B00089")
  box_color = "#F2F2F2"
  
  p2 <- contrast_correct %>%
    mutate(contrast = factor(contrast, levels = c("propranolol - bisoprolol", "placebo - bisoprolol", "placebo - propranolol"))) %>% 
    mutate(y_val = as.numeric(as.factor(contrast))) %>% 
    ggplot(aes(col = contrast, show.legend = FALSE)) +
    geom_point(aes(y = y_val, x = estimate),size = 5, show.legend = FALSE)+
    geom_pointrange(aes(y = y_val, x = estimate, xmin = lower.CL, xmax = upper.CL), linewidth = 1.5, show.legend = FALSE)+
    labs(
      title = "",
      y = "",
      x = "Difference in Predicted Confidence"
    ) +
    coord_cartesian(xlim = c(-0.325,0), ylim = c(0.5,3.4))+
    geom_segment(aes(x = 0, xend = 0, y = -2.3,yend = 8.3), linetype = 2, linewidth = line_width, col = "#1200A8", show.legend = FALSE)+
    theme_classic(base_size = 12) +
    scale_y_continuous(breaks = c(1,2,3), labels = c("prop - biso","plac - biso","plac - prop"))+
    theme(strip.background = element_rect(fill = box_color, color = "white"),  # Default box style
          axis.line.y = element_blank(), #)+        # Remove the y-axis line
          axis.ticks.y = element_blank())+
    scale_color_manual(values = colors)+  # Ensures no padding below x-axis
    theme(plot.margin = margin(t = 6, r = 6, b = 6, l = 6), "pt")+
    theme(
      plot.title = element_text(face = "bold", size = fontsize),
      axis.text.y = element_text(angle = 25, hjust = 0.5, size = fontsize),
      axis.text.x = element_text(size = fontsize),
      axis.title.x = element_text(size = fontsize),
      axis.line.x = element_line(linewidth = 1),
      legend.position = "none"
    )
  
  
  ## RRST PLOTS
  RSpreds <- ggpredict(RRST_conf_model, terms = c("Drug", "Accuracy"))
  RSpreds$group <- factor(RSpreds$group, levels=c('Correct', 'Incorrect'))
  
  # Create the plot for RRST predicted confidence
  p3 <-ggplot(RSpreds, aes(x = x, y = predicted, group = group, color = group)) +
    geom_line(size = 1.3) +
    geom_point(size = 4) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = ribbon_col), alpha = 0.15, outline.type = "both", show.legend = FALSE) +
    scale_color_manual(values = binary_colourmap) +
    scale_fill_manual(values = ribbon_col) +
    labs(
      x = "",
      y = "",
      color = "Accuracy",
      fill = ""
    ) +
    theme_classic(base_size = 12) +
    coord_cartesian(ylim = c(0.43,0.71))+
    scale_x_discrete(expand = c(0.1,0))+
    scale_y_continuous(breaks = c(0.5,0.6,0.7), labels = c("0.5","0.6","0.7"))+
    theme(plot.margin = margin(t = 6, r = 6, b = 6, l = 6), "pt")+
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = fontsize),
      legend.text = element_text(size = fontsize),
      axis.text = element_text(size = fontsize),
      axis.text.x = element_text(size = fontsize), #, hjust = 0.6), 
      axis.title = element_blank(),
      plot.title = element_text(size = 12, face = "bold"),
      axis.line.y = element_line(linewidth = 1),
      axis.line.x = element_line(linewidth = 1)
    )
  
  # Obtain estimated marginal means for each Drug and Accuracy combination
  RSemm <- emmeans(RRST_conf_model, ~ Drug | Accuracy)
  
  # Perform pairwise comparisons between Drugs within each Accuracy level
  RScontrast_results <- contrast(RSemm, method = "pairwise")
  
  # Convert the contrast results to a data frame
  RScontrast_df <- as.data.frame(RScontrast_results)
  
  # Filter contrasts for Accuracy == Correct
  RScontrast_correct <- subset(RScontrast_df, Accuracy == "Correct")
  
  # Calculate 95% confidence intervals
  # Since df is infinite (df = Inf), we'll use the normal distribution (z-distribution)
  z_value <- qnorm(0.975)  # For 95% CI
  
  RScontrast_correct$lower.CL <- RScontrast_correct$estimate - z_value * RScontrast_correct$SE
  RScontrast_correct$upper.CL <- RScontrast_correct$estimate + z_value * RScontrast_correct$SE
  
  
  # Create the RRST contrasts plot
  p4 <- RScontrast_correct %>%
    mutate(contrast = factor(contrast, levels = c("propranolol - bisoprolol", "placebo - bisoprolol", "placebo - propranolol"))) %>% 
    mutate(y_val = as.numeric(as.factor(contrast))) %>% 
    ggplot(aes(col = contrast, show.legend = FALSE)) +
    geom_point(aes(y = y_val, x = estimate),size = 5, show.legend = FALSE)+
    geom_pointrange(aes(y = y_val, x = estimate, xmin = lower.CL, xmax = upper.CL), linewidth = 1.5, show.legend = FALSE)+
    labs(
      title = "",
      y = "",
      x = "Difference in Predicted Confidence"
    ) +
    coord_cartesian(xlim = c(-0.1,0.1), ylim = c(0.5,3.4))+
    geom_segment(aes(x = 0, xend = 0, y = -2.3,yend = 8.3), linetype = 2, linewidth = line_width, col = "#1200A8", show.legend = FALSE)+
    theme_classic(base_size = 12) +
    scale_y_continuous(breaks = c(1,2,3), labels = c("","",""))+
    theme(strip.background = element_rect(fill = box_color, color = "white"),  # Default box style
          axis.line.y = element_blank(),  # Remove the y-axis line
          axis.ticks.y = element_blank())+
    scale_color_manual(values = colors)+  # Ensures no padding below x-axis
    theme(plot.margin = margin(t = 6, r = 6, b = 6, l = 6), "pt")+
    theme(
      plot.title = element_text(face = "bold", size = fontsize),
      axis.text.y = element_text(angle = 0, hjust = 0.5, size = fontsize),
      axis.text.x = element_text(size = fontsize),
      axis.title.x = element_text(size = fontsize),
      axis.line.x = element_line(linewidth = 1),
      legend.position = "none"
    )
  
  combined_plot_top <- ggarrange(
    p1+p3,
    common.legend = TRUE, legend = "bottom",
    heights = c(1,1),
    widths = c(1,1)
  )
  
  combined_plot_bottom <- ggarrange(
    p2+p4,
    heights = c(1,1),
    widths = c(2,2)
  )
  
  combined_plot_bottom <- annotate_figure(combined_plot_bottom, 
                                          top = text_grob("Pairwise Differences Between Drugs (Correct Trials)", 
                                                          face = "bold", size = fontsize))
  
  combined_plot_full <- ggarrange(
    combined_plot_top, NULL, combined_plot_bottom,
    ncol = 1,
    nrow = 3,
    heights = c(3,0.1,2),
    align = "v"
  )
  
  # Display the combined plot
  print(combined_plot_full)
  
  ggsave(here::here("Figures","plot3.tiff"), combined_plot_full, device = "tiff", width = 12, height = 10, units = "in", dpi = 400)
  ggsave(here::here("Figures","plot3.png"), combined_plot_full, device = "png", width = 12, height = 10, units = "in", dpi = 400)
  
}

# Supplementary Figures

# Supp Fig 1: HRV metrics per drug:

plot_s1 = function(){
  
  HRV <- read_csv(here::here("data", "alldrugs_full_HRV.csv"))
  
  # Summarise the HRV by drug and subject
  subject_summary2 <- HRV %>%
    group_by(Drug, SubNo) %>%
    mutate(Drug = ifelse(Drug == "plac","placebo",ifelse(Drug == "biso","bisoprolol",ifelse(Drug == "prop","propranolol",NA))))
  
  # 1. Subplot: Boxplot + Scatter for Mean BPM (rsECG)
  meanHR_plot <- HRV_scatter(subject_summary2$meanHR, subject_summary2$Drug, "Mean Heart Rate")
  
  # 2. Subplot: Boxplot + Scatter for RMSSD (rsECG)
  RMSSD_plot <- HRV_scatter(subject_summary2$RMSSD, subject_summary2$Drug, "Heart Rate Variability (RMSSD)")
  
  # 3. Subplot: Boxplot + Scatter for LF/HF Ratio (rsECG)
  LFHF_plot <- HRV_scatter(subject_summary2$LFHF, subject_summary2$Drug, "Heart Rate Variability (LF/HF Ratio)")
  
  # Combine and display plot
  plots1 = meanHR_plot+RMSSD_plot+LFHF_plot
  print(plots1)
  
  # save plot
  ggsave(here::here("Figures","Supplementary1_HRV.tiff"),plots1,dpi = 400, width = 18,height = 10)
  
}

# Supp Figs 2 and 3: Individual subject plots of psychometrics with raw data overlay:
plot_s2 = function(){
  
  #load two models and their data
  rrst_data = read.csv(here::here("data","cleaned","RRST.csv"))
  hrd_data = read.csv(here::here("data","cleaned","HRD.csv"))
  
  HRD <- readRDS(here::here("STAN models","HRD.RDS"))
  RRST <- readRDS(here::here("STAN models","RRST.RDS"))
  
  
  # function to get and plot a single subjects' estimates as draws (HRD). This is the function that is looped over.
  
  n_draws = 50
  draws = sample(1:4000,n_draws)
  get_estimates_hrd = function(hrd_data,id,draws){
    
    fontsize = 14
    line_width = 1
    colors = c("#009F73","#1200A8","#B00089")
    
    subject_hrd = hrd_data %>% filter(sub == id)
    
    # making bins for the raw data such that one can see the increase in the probability of response and not just 0 and 1's.
    breaks <- unique(quantile(subject_hrd$Alpha, probs = seq(0, 1, length.out = 7)))
    
    # calculate mean and se of the raw data for each condition
    raw_data = subject_hrd %>% 
      mutate(
        bins = cut(Alpha, breaks = breaks, include.lowest = TRUE, labels = FALSE),
        Alpha = sapply(bins, function(bin) mean(Alpha[bins == bin])),
      )%>% group_by(drugs, Alpha,bins) %>% 
      summarize(prob = sum(resp)/n(), prob_se = mean(resp)*(1-mean(resp)) / sqrt(n()), n = n()) %>% 
      mutate(Alpha_low = sapply(bins, function(bin) breaks[bin]),
             Alpha_high = sapply(bins, function(bin) breaks[bin + 1])) %>% 
      mutate(draw = NA)
    
    ## Getting the subject level psychometric functions
    
    #threshold
    alphas = as_draws_df(HRD$draws(paste0("alpha[",subject_hrd$full_trials,"]"))) %>% 
      select(-contains(".")) %>% mutate(draw = 1:n()) %>% 
      pivot_longer(-draw, names_to = "parameter") %>% 
      mutate(full_trials = as.numeric(str_extract(parameter, "(?<=\\[)\\d+(?=\\])"))) %>% 
      mutate(parameter = "alpha")
    
    #slope
    betas = as_draws_df(HRD$draws(paste0("beta[",subject_hrd$full_trials,"]"))) %>% 
      select(-contains(".")) %>% mutate(draw = 1:n()) %>% 
      pivot_longer(-draw, names_to = "parameter") %>% 
      mutate(full_trials = as.numeric(str_extract(parameter, "(?<=\\[)\\d+(?=\\])"))) %>% 
      mutate(parameter = "beta")
    
    #lapserate
    lapses = as_draws_df(HRD$draws(paste0("lapse[",subject_hrd$full_trials,"]"))) %>% 
      select(-contains(".")) %>% mutate(draw = 1:n()) %>% 
      pivot_longer(-draw, names_to = "parameter") %>% 
      mutate(full_trials = as.numeric(str_extract(parameter, "(?<=\\[)\\d+(?=\\])"))) %>% 
      mutate(parameter = "lapse")
    
    # psychometric eqaution
    psycho_hrd = function(x,alpha,beta,lapse){
      return(lapse + (1 - 2 * lapse) * ((0.5+0.5 * pracma::erf((x-alpha) / (beta * sqrt(2))))))
    }
    
    #combine threshold, slope and lapse-rate and make them into a format that is workable
    dfq = inner_join(rbind(alphas,betas,lapses) %>% filter(draw %in% draws),
                     subject_hrd %>% mutate(listenBPM = NULL, X = NULL, Alpha = NULL, resp = NULL)) %>% mutate(trials = NULL) %>% 
      distinct() %>% select(draw,parameter,drugs,value) %>% pivot_wider(names_from = "drugs", values_from = "value")
    
    min_Alpha = min(subject_hrd$Alpha)
    max_Alpha = max(subject_hrd$Alpha)
    
    # now some subjects do not have all sessions so need to account for that:
    
    # all sessions:
    if(sum(c("BISO","PROP","PLACEBO") %in% colnames(dfq))== 3){
      
      plot = dfq %>% unnest(cols = c(BISO, PLACEBO, PROP)) %>% 
        pivot_wider(values_from = c("PLACEBO","BISO","PROP"),names_from = "parameter") %>% unnest() %>% 
        group_by(draw) %>% 
        summarize(p_placebo = list(psycho_hrd(seq(min_Alpha,max_Alpha, by = 0.1),PLACEBO_alpha,PLACEBO_beta,PLACEBO_lapse)),
                  p_biso = list(psycho_hrd(seq(min_Alpha,max_Alpha,by = 0.1),BISO_alpha,BISO_beta,BISO_lapse)),
                  p_prop = list(psycho_hrd(seq(min_Alpha,max_Alpha,by = 0.1),PROP_alpha,PROP_beta,PROP_lapse)),
                  x = list(seq(min_Alpha,max_Alpha,by = 0.1))) %>% unnest() %>% pivot_longer(cols = c("p_placebo","p_biso","p_prop"),names_to = "drugs",values_to = "prob") %>% 
        mutate(drugs = ifelse(drugs == "p_biso","bisoprolol",ifelse(drugs == "p_placebo","placebo","propranolol"))) %>% 
        ggplot(aes(x = x, y = prob, col = drugs, group = interaction(draw,drugs)))+
        geom_line(alpha = 0.2)+
        geom_pointrange(data = raw_data %>% mutate(drugs = ifelse(drugs == "BISO","bisoprolol",ifelse(drugs == "PLACEBO","placebo","propranolol"))), aes(x = Alpha, y = prob, ymin = prob-2*prob_se, ymax = prob+2*prob_se, fill = drugs), col = "black", shape = 21,size = 0.7, linewidth = line_width)+
        geom_pointrange(data = raw_data %>% mutate(drugs = ifelse(drugs == "BISO","bisoprolol",ifelse(drugs == "PLACEBO","placebo","propranolol"))), aes(x = Alpha, y = prob, xmin = Alpha_low, xmax = Alpha_high, fill = drugs), col = "black", shape = 21,size = 0.7, linewidth = line_width)+
        # xlab("Stimulus intensity (ΔBPM)")+
        scale_y_continuous("P(Response = faster | ΔBPM)",breaks = c(0,0.5,1),labels = c("0.0","0.5","1.0"))+
        scale_x_continuous(" ",breaks = c(-20,0,20),labels = c("-20","0","20"))+
        # ylab("")+
        theme_classic()+
        scale_color_manual(values = colors)+
        scale_fill_manual(values = colors)+
        theme(text = element_text(size = fontsize),           # All text
              axis.title = element_text(size = fontsize),     # Axis titles
              axis.text.x = element_text(size = fontsize),      # Axis tick labels
              axis.text.y = element_text(size = fontsize),      # Axis tick labels
              legend.text = element_text(size = fontsize),    # Legend text
              legend.title = element_text(size = fontsize),   # Legend title
              plot.title = element_text(size = fontsize),
              axis.line.y = element_line(linewidth = line_width),
              axis.line.x = element_line(linewidth = line_width)
        )+
        scale_color_manual(values = c("placebo" = "#1200A8", "propranolol" = "#B00089", "bisoprolol" = "#009F73"),
                           breaks = c("placebo", "propranolol", "bisoprolol"))+
        scale_fill_manual(values = c("placebo" = "#1200A8", "propranolol" = "#B00089", "bisoprolol" = "#009F73"),
                          breaks = c("placebo", "propranolol", "bisoprolol"))
      
      return(plot)
      
      #only Propanolol and placebo
    }else if(sum(c("PROP","PLACEBO") %in% colnames(dfq)) == 2){
      
      
      colors = c("#1200A8","#B00089")
      
      plot = dfq %>% pivot_wider(values_from = c("PLACEBO","PROP"),names_from = "parameter") %>% unnest() %>% 
        group_by(draw) %>% 
        summarize(p_placebo = list(psycho_hrd(seq(min_Alpha,max_Alpha, by = 0.1),PLACEBO_alpha,PLACEBO_beta,PLACEBO_lapse)),
                  p_prop = list(psycho_hrd(seq(min_Alpha,max_Alpha, by = 0.1),PROP_alpha,PROP_beta,PROP_lapse)),
                  x = list(seq(min_Alpha,max_Alpha, by = 0.1))) %>% unnest() %>% pivot_longer(cols = c("p_placebo","p_prop"),names_to = "drugs",values_to = "prob") %>% 
        mutate(drugs = ifelse(drugs == "p_biso","bisoprolol",ifelse(drugs == "p_placebo","placebo","propranolol"))) %>% 
        ggplot(aes(x = x, y = prob, col = drugs, group = interaction(draw,drugs)))+
        geom_line(alpha = 0.2, show.legend = F)+
        geom_pointrange(data = raw_data %>% mutate(drugs = ifelse(drugs == "BISO","bisoprolol",ifelse(drugs == "PLACEBO","placebo","propranolol"))), aes(x = Alpha, y = prob, ymin = prob-2*prob_se, ymax = prob+2*prob_se, fill = drugs), col = "black", shape = 21,size = 0.7, linewidth = line_width, show.legend = F)+
        geom_pointrange(data = raw_data %>% mutate(drugs = ifelse(drugs == "BISO","bisoprolol",ifelse(drugs == "PLACEBO","placebo","propranolol"))), aes(x = Alpha, y = prob, xmin = Alpha_low, xmax = Alpha_high, fill = drugs), col = "black", shape = 21,size = 0.7, linewidth = line_width, show.legend = F)+
        # xlab("Stimulus intensity (ΔBPM)")+
        scale_y_continuous("P(Response = faster | ΔBPM)",breaks = c(0,0.5,1),labels = c("0.0","0.5","1.0"))+
        scale_x_continuous(" ",breaks = c(-20,0,20),labels = c("-20","0","20"))+
        # ylab("")+
        theme_classic()+
        scale_color_manual(values = colors)+
        scale_fill_manual(values = colors)+
        theme(text = element_text(size = fontsize),           # All text
              axis.title = element_text(size = fontsize),     # Axis titles
              axis.text.x = element_text(size = fontsize),      # Axis tick labels
              axis.text.y = element_text(size = fontsize),      # Axis tick labels
              legend.text = element_text(size = fontsize),    # Legend text
              legend.title = element_text(size = fontsize),   # Legend title
              plot.title = element_text(size = fontsize),
              axis.line.y = element_line(linewidth = line_width),
              axis.line.x = element_line(linewidth = line_width))
      
      return(plot)
      #only Bisoprolol and placebo
    }else if(sum(c("BISO","PLACEBO") %in% colnames(dfq)) == 2){
      
      
      colors = c("#1200A8","#009F73")
      
      plot = dfq %>% pivot_wider(values_from = c("PLACEBO","BISO"),names_from = "parameter") %>% unnest() %>% 
        group_by(draw) %>% 
        summarize(p_placebo = list(psycho_hrd(seq(min_Alpha,max_Alpha, by = 0.1),PLACEBO_alpha,PLACEBO_beta,PLACEBO_lapse)),
                  p_biso = list(psycho_hrd(seq(min_Alpha,max_Alpha, by = 0.1),BISO_alpha,BISO_beta,BISO_lapse)),
                  x = list(seq(min_Alpha,max_Alpha, by = 0.1))) %>% unnest() %>% pivot_longer(cols = c("p_placebo","p_biso"),names_to = "drugs",values_to = "prob") %>% 
        mutate(drugs = ifelse(drugs == "p_biso","bisoprolol",ifelse(drugs == "p_placebo","placebo","propranolol"))) %>% 
        ggplot(aes(x = x, y = prob, col = drugs, group = interaction(draw,drugs)))+
        geom_line(alpha = 0.2, show.legend = F)+
        geom_pointrange(data = raw_data %>% mutate(drugs = ifelse(drugs == "BISO","bisoprolol",ifelse(drugs == "PLACEBO","placebo","propranolol"))), aes(x = Alpha, y = prob, ymin = prob-2*prob_se, ymax = prob+2*prob_se, fill = drugs), col = "black", shape = 21,size = 0.7, linewidth = line_width, show.legend = F)+
        geom_pointrange(data = raw_data %>% mutate(drugs = ifelse(drugs == "BISO","bisoprolol",ifelse(drugs == "PLACEBO","placebo","propranolol"))), aes(x = Alpha, y = prob, xmin = Alpha_low, xmax = Alpha_high, fill = drugs), col = "black", shape = 21,size = 0.7, linewidth = line_width, show.legend = F)+
        # xlab("Stimulus intensity (ΔBPM)")+
        scale_y_continuous("P(Response = faster | ΔBPM)",breaks = c(0,0.5,1),labels = c("0.0","0.5","1.0"))+
        scale_x_continuous(" ",breaks = c(-20,0,20),labels = c("-20","0","20"))+
        theme_classic()+
        scale_color_manual(values = colors)+
        scale_fill_manual(values = colors)+
        theme(text = element_text(size = fontsize),           # All text
              axis.title = element_text(size = fontsize),     # Axis titles
              axis.text.x = element_text(size = fontsize),      # Axis tick labels
              axis.text.y = element_text(size = fontsize),      # Axis tick labels
              legend.text = element_text(size = fontsize),    # Legend text
              legend.title = element_text(size = fontsize),   # Legend title
              plot.title = element_text(size = fontsize),
              axis.line.y = element_line(linewidth = line_width),
              axis.line.x = element_line(linewidth = line_width))
      
      return(plot)
      #only Bisoprolol and Propanolol
    }else if(sum(c("BISO","PROP") %in% colnames(dfq)) == 2){
      
      
      colors = c("#009F73","#B00089")
      
      plot = dfq %>% pivot_wider(values_from = c("PROP","BISO"),names_from = "parameter") %>% unnest() %>% 
        group_by(draw) %>% 
        summarize(p_prop = list(psycho_hrd(seq(min_Alpha,max_Alpha, by = 0.1),PROP_alpha,PROP_beta,PROP_lapse)),
                  p_biso = list(psycho_hrd(seq(min_Alpha,max_Alpha, by = 0.1),BISO_alpha,BISO_beta,BISO_lapse)),
                  x = list(seq(min_Alpha,max_Alpha, by = 0.1))) %>% unnest() %>% pivot_longer(cols = c("p_prop","p_biso"),names_to = "drugs",values_to = "prob") %>% 
        mutate(drugs = ifelse(drugs == "p_biso","bisoprolol",ifelse(drugs == "p_placebo","placebo","propranolol"))) %>% 
        ggplot(aes(x = x, y = prob, col = drugs, group = interaction(draw,drugs)))+
        geom_line(alpha = 0.2, show.legend = F)+
        geom_pointrange(data = raw_data %>% mutate(drugs = ifelse(drugs == "BISO","bisoprolol",ifelse(drugs == "PLACEBO","placebo","propranolol"))), aes(x = Alpha, y = prob, ymin = prob-2*prob_se, ymax = prob+2*prob_se, fill = drugs), col = "black", shape = 21,size = 0.7, linewidth = line_width, show.legend = F)+
        geom_pointrange(data = raw_data %>% mutate(drugs = ifelse(drugs == "BISO","bisoprolol",ifelse(drugs == "PLACEBO","placebo","propranolol"))), aes(x = Alpha, y = prob, xmin = Alpha_low, xmax = Alpha_high, fill = drugs), col = "black", shape = 21,size = 0.7, linewidth = line_width, show.legend = F)+
        # xlab("Stimulus intensity (ΔBPM)")+
        scale_y_continuous("P(Response = faster | ΔBPM)",breaks = c(0,0.5,1),labels = c("0.0","0.5","1.0"))+
        scale_x_continuous(" ",breaks = c(-20,0,20),labels = c("-20","0","20"))+
        theme_classic()+
        scale_color_manual(values = colors)+
        scale_fill_manual(values = colors)+
        theme(text = element_text(size = fontsize),           # All text
              axis.title = element_text(size = fontsize),     # Axis titles
              axis.text.x = element_text(size = fontsize),      # Axis tick labels
              axis.text.y = element_text(size = fontsize),      # Axis tick labels
              legend.text = element_text(size = fontsize),    # Legend text
              legend.title = element_text(size = fontsize),   # Legend title
              plot.title = element_text(size = fontsize),
              axis.line.y = element_line(linewidth = line_width),
              axis.line.x = element_line(linewidth = line_width))
      
      return(plot)
    }
    
    
  }
  
  # Same for RRST as for HRD
  get_estimates_rrst = function(rrst_data,id,draws){
    
    
    fontsize = 14
    line_width = 1
    colors = c("#009F73","#1200A8","#B00089")
    subject_rrst = rrst_data %>% filter(sub == id)
    
    #again binning the stimulus value:
    breaks <- unique(quantile(subject_rrst$Stim, probs = seq(0, 1, length.out = 5)))
    
    #raw datapoints dataframe
    raw_data = subject_rrst %>% 
      mutate(
        bins = cut(Stim, breaks = breaks, include.lowest = TRUE, labels = FALSE),
        Stim = sapply(bins, function(bin) mean(Stim[bins == bin])),
      )%>% group_by(drugs, Stim,bins) %>% 
      summarize(prob = sum(Resp)/n(), prob_se = mean(Resp)*(1-mean(Resp)) / sqrt(n()), n = n()) %>% 
      mutate(Alpha_low = sapply(bins, function(bin) breaks[bin]),
             Alpha_high = sapply(bins, function(bin) breaks[bin + 1])) %>% 
      mutate(draw = NA)
    
    
    #threshold psychometric
    alphas = as_draws_df(RRST$draws(paste0("alpha[",subject_rrst$full_trials,"]"))) %>% 
      select(-contains(".")) %>% mutate(draw = 1:n()) %>% 
      pivot_longer(-draw, names_to = "parameter") %>% 
      mutate(full_trials = as.numeric(str_extract(parameter, "(?<=\\[)\\d+(?=\\])"))) %>% 
      mutate(parameter = "alpha")
    
    #slope psychometric
    betas = as_draws_df(RRST$draws(paste0("beta[",subject_rrst$full_trials,"]"))) %>% 
      select(-contains(".")) %>% mutate(draw = 1:n()) %>% 
      pivot_longer(-draw, names_to = "parameter") %>% 
      mutate(full_trials = as.numeric(str_extract(parameter, "(?<=\\[)\\d+(?=\\])"))) %>% 
      mutate(parameter = "beta")
    
    #lapse-rate psychometric
    lapses = as_draws_df(RRST$draws(paste0("lapse[",subject_rrst$full_trials,"]"))) %>% 
      select(-contains(".")) %>% mutate(draw = 1:n()) %>% 
      pivot_longer(-draw, names_to = "parameter") %>% 
      mutate(full_trials = as.numeric(str_extract(parameter, "(?<=\\[)\\d+(?=\\])"))) %>% 
      mutate(parameter = "lapse")
    
    # psychometric equation
    psycho_rrst = function(x,alpha,beta,lapse){
      return(0.5 + (1 - 0.5 - lapse) * (1-exp(-10^(beta * (x-alpha)))))
    }
    
    min_stim = min(subject_rrst$Stim)
    max_stim = max(subject_rrst$Stim)
    
    ## ALL subjects had all sessions so only need one plotting function
    
    #combine the parameters and plot
    inner_join(rbind(alphas,betas,lapses) %>% filter(draw %in% draws),
               subject_rrst %>% mutate(listenBPM = NULL, X = NULL, Stim = NULL, Resp = NULL, ConfResp = NULL)) %>% mutate(trials = NULL) %>% 
      distinct() %>% select(draw,parameter,drugs,value) %>% distinct() %>% 
      pivot_wider(names_from = "drugs", values_from = "value") %>% unnest() %>% 
      pivot_wider(values_from = c("PLACEBO","BISO","PROP"),names_from = "parameter") %>% unnest() %>% 
      group_by(draw) %>% 
      summarize(p_placebo = list(psycho_rrst(seq(min_stim,max_stim,by = 0.1),PLACEBO_alpha,PLACEBO_beta,PLACEBO_lapse)),
                p_biso = list(psycho_rrst(seq(min_stim,max_stim,by = 0.1),BISO_alpha,BISO_beta,BISO_lapse)),
                p_prop = list(psycho_rrst(seq(min_stim,max_stim,by = 0.1),PROP_alpha,PROP_beta,PROP_lapse)),
                x = list(seq(min_stim,max_stim,by = 0.1))) %>% 
      unnest() %>% 
      pivot_longer(cols = c("p_placebo","p_biso","p_prop"),names_to = "drugs",values_to = "prob") %>% 
      mutate(drugs = as.factor(ifelse(drugs == "p_biso","bisoprolol",ifelse(drugs == "p_placebo","placebo","propranolol")))) %>% 
      ggplot(aes(x = x, y = prob, col = drugs, group = interaction(draw,drugs)))+
      geom_line(alpha = 0.2)+
      geom_pointrange(data = raw_data%>% mutate(drugs = ifelse(drugs == "BISO","bisoprolol",ifelse(drugs == "PLACEBO","placebo","propranolol"))), aes(x = Stim, y = prob, ymin = prob-2*prob_se, ymax = prob+2*prob_se, fill = drugs), col = "black", shape = 21,size = 0.7, linewidth = line_width)+
      geom_pointrange(data = raw_data%>% mutate(drugs = ifelse(drugs == "BISO","bisoprolol",ifelse(drugs == "PLACEBO","placebo","propranolol"))), aes(x = Stim, y = prob, xmin = Alpha_low, xmax = Alpha_high, fill = drugs), col = "black", shape = 21,size = 0.7, linewidth = line_width)+
      coord_cartesian(ylim = c(0.5,1))+
      ylab("P(correct | % RRes)")+
      scale_y_continuous("P(correct | % RRes)", breaks = c(0.5,0.75,1.0), labels =  c("0.5","0.75","1.0"))+
      scale_x_continuous(" ", breaks = c(8,12,16), labels =  c("8","12","16"))+
      theme_classic()+
      theme_classic()+
      scale_color_manual(values = colors)+
      scale_fill_manual(values = colors)+
      theme(text = element_text(size = fontsize),           # All text
            axis.title = element_text(size = fontsize),     # Axis titles
            axis.text.x = element_text(size = fontsize),      # Axis tick labels
            axis.text.y = element_text(size = fontsize),      # Axis tick labels
            legend.text = element_text(size = fontsize),    # Legend text
            legend.title = element_text(size = fontsize),   # Legend title
            plot.title = element_text(size = fontsize),
            axis.line.y = element_line(linewidth = line_width),
            axis.line.x = element_line(linewidth = line_width)
      )+
      scale_color_manual(values = c("placebo" = "#1200A8", "propranolol" = "#B00089", "bisoprolol" = "#009F73"),
                         breaks = c("placebo", "propranolol", "bisoprolol"))+
      scale_fill_manual(values = c("placebo" = "#1200A8", "propranolol" = "#B00089", "bisoprolol" = "#009F73"),
                        breaks = c("placebo", "propranolol", "bisoprolol"))
    
  }
  
  
  ## Loop through the HRD
  plot_hrd = list()
  for(i in (unique(hrd_data$sub))){
    print(paste0("number ", i, " in: ", length((unique(hrd_data$sub)))))
    plot_hrd[[i]] = get_estimates_hrd(hrd_data,i,draws)
  }
  
  #combine the plots
  combined_plot <- wrap_plots(plot_hrd, ncol = 10, nrow = 5)+
    plot_layout(axes = "collect_y",
                guides = "collect")+
    labs(tag = "Stimulus intensity (ΔBPM)") &
    theme(
      plot.tag.position = c(-3.3,0.29),
      legend.position = "bottom",
      axis.title = element_text(size = 22),  # Ensuring consistency
      legend.text = element_text(size = 22),  # Match legend text size to axis title
      plot.tag = element_text(size = 22),  # Match tag size to axis title
      legend.title = element_blank()
    )
  combined_plot
  #save it
  ggsave(here::here("Figures","Supplementary2_HRD.tiff"), combined_plot, width = 20, height = 10, units = "in", dpi = 400)
  
  ## Loop through the RRST
  plot_rrst = list()
  for(i in (unique(rrst_data$sub))){
    print(paste0("number ", i, " in: ", length((unique(rrst_data$sub)))))
    plot_rrst[[i]] = get_estimates_rrst(rrst_data,i,draws)
  }
  
  # combine it
  combined_plot_rrst <- wrap_plots(plot_rrst, ncol = 10, nrow = 5)+
    plot_layout(guides = "collect",
                axes = "collect_y")+
    labs(tag = "Stimulus intensity (% RRes)") &
    theme(
      plot.tag.position = c(-3.3,0.29),
      legend.position = "bottom",
      axis.title = element_text(size = 22),  # Ensuring consistency
      legend.text = element_text(size = 22),  # Match legend text size to axis title
      plot.tag = element_text(size = 22),  # Match tag size to axis title
      legend.title = element_blank()
    )
  combined_plot_rrst
  #save it
  ggsave(here::here("Figures","Supplementary3_RRST.tiff"), combined_plot_rrst, width = 20, height = 10, units = "in", dpi = 400)
  
  
}


# Supp Figs 4 & 5 : group level plots of meta cognition model with raw data overlay:
plot_s3 = function(){
  
  RRST_trial_data <- read_csv(here::here("data","cleaned","RRST.csv")) %>% 
    filter(ConfResp >= 0 & ConfResp <= 100) %>%
    mutate(Conf = ConfResp/100,
           Stimulus = Stim,
           Drug = factor(drugs, levels = c("PLACEBO", "PROP", "BISO"), labels = c("placebo", "propranolol", "bisoprolol")),
           Drug = relevel(Drug, ref = "placebo"),
           Condition = drugs,
           Accuracy = factor(Resp, levels = c(1, 0), labels = c("Correct", "Incorrect"))
    )
  
  
  # Cleaning RRST data
  RRST_trial_data <- clean_and_report(RRST_trial_data, "RRST_trial_data")
  
  
  ## Models
  
  # RRST - Ordered beta regression confidence model
  RRST_conf_model <- glmmTMB(
    Conf ~ Stimulus + Drug * Accuracy + (1 + Stimulus + Accuracy | subject),
    data = RRST_trial_data,
    family = ordbeta(),
    control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000))
  )
  
  # summarise effects
  summary(RRST_conf_model)
  
  
  ## plotting the estimated means and overlaying a grand average data.
  P1 = (ggemmeans(RRST_conf_model, terms = c("Stimulus","Accuracy" ,"Drug")))
  
  
  subj_data <- RRST_trial_data %>%
    group_by(subject, Accuracy, Stimulus, Drug) %>%
    summarize(mean_conf = mean(Conf), .groups = "drop") %>%
    group_by(Accuracy, Stimulus, Drug) %>%
    summarize(
      mean = mean(mean_conf),
      se = sd(mean_conf) / sqrt(n()),
      q5 = mean - 2 * se,
      q95 = mean + 2 * se,
      .groups = "drop"
    ) %>%
    mutate(facet = Drug)
  
  
  combined_plot_rrst = data.frame(P1) %>% 
    mutate(Accuracy = as.factor(group)) %>% 
    ggplot()+
    geom_line(aes(x = x, y = predicted, ymin = conf.low,, ymax = conf.high, fill = Accuracy))+
    geom_ribbon(aes(x = x, y = predicted, ymin = conf.low,, ymax = conf.high, fill = Accuracy), alpha = 0.25)+
    geom_pointrange(data = subj_data, aes(x = Stimulus, y = mean, ymin = q5,ymax = q95, fill = as.factor(Accuracy)),
                    shape = 21,
                    col = "black", 
                    position = position_dodge(width = 1))+
    scale_color_manual(values = c("green","red"))+
    scale_fill_manual(values = c("green","red"))+
    facet_wrap(~facet)+
    labs(x = "Stimulus intensity (% RRes)",
         y = "P(Response = faster | % RRes)")+
    scale_x_continuous(breaks = c(4.25,8.5,12.75,17), labels = c(25,50,75,100))+
    theme_classic(base_size = 16)
  
  
  ggsave(here::here("figures","revisions","Supplementary4_RRST.tiff"), combined_plot_rrst, width = 20, height = 10, units = "in", dpi = 400)
  
  
  HRD_trl_data = read.csv(here::here("data","cleaned","HRD.csv")) %>% 
    mutate(Conf = Confidence/100,
           drugs = as.factor(drugs),
           Drug = factor(drugs, levels = c("PLACEBO", "PROP", "BISO"), labels = c("placebo", "propranolol", "bisoprolol")),
           Drug = relevel(Drug, ref = "placebo"),
           ResponseCorrect = factor(ResponseCorrect, levels = c("True", "False"), labels = c("Correct", "Incorrect")),
           BPM_scaled = scale(listenBPM),
           Condition = drugs
    ) %>% 
    filter(subject != "sub_4049")
  
  # Cleaning HRD data
  HRD_trl_data <- clean_and_report(HRD_trl_data, "HRD_trl_data")
  
  # HRD - Ordered beta regression confidence model
  HRD_conf_model <- glmmTMB(
    Conf ~ Drug * ResponseCorrect  + BPM_scaled +  (1  + ResponseCorrect + BPM_scaled | subject),
    data = HRD_trl_data,
    family = ordbeta(),
    start = list(psi = c(0, 1)),
    control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000))
  )
  
  summary(HRD_conf_model)
  
  p1  = (ggemmeans(HRD_conf_model, terms = c("ResponseCorrect" ,"Drug")))
  
  
  breaks <- unique(quantile(HRD_trl_data$Alpha, probs = seq(0, 1, length.out = 9)))
  
  
  
  subj_data <-  HRD_trl_data %>% 
    group_by(Drug) %>%
    mutate(
      bins = cut(Alpha, breaks = breaks, include.lowest = TRUE, labels = FALSE),
      Alpha = sapply(bins, function(bin) mean(Alpha[bins == bin])),
    )%>% 
    # group_by(subject, ResponseCorrect, Drug, Alpha,bins) %>% 
    mutate(Alpha_low = sapply(bins, function(bin) breaks[bin]),
           Alpha_high = sapply(bins, function(bin) breaks[bin + 1])) %>% 
    # mutate(Alpha = cut(Alpha,11))%>% 
    group_by(subject, ResponseCorrect, Alpha, Drug) %>%
    summarize(mean_conf = mean(Conf, na.rm = T), .groups = "drop") %>%
    group_by(ResponseCorrect, Alpha, Drug) %>%
    summarize(
      mean = mean(mean_conf, na.rm = T),
      se = sd(mean_conf, na.rm = T) / sqrt(n()),
      q5 = mean - 2 * se,
      q95 = mean + 2 * se,
      .groups = "drop"
    ) %>%
    mutate(facet = Drug)
  
  
  preds <- data.frame(p1) %>%
    rename(
      Accuracy = x,
      facet = group
    ) %>%
    tidyr::crossing(x = seq(min(subj_data$Alpha)-1, max(subj_data$Alpha)+1, by = 5))
  
  # preds = data.frame(facet = c("placebo","placebo","propranolol","propranolol","bisoprolol","bisoprolol"),
  #                    Accuracy = c("Correct","Incorrect","Correct","Incorrect","Correct","Incorrect"),
  #                    predicted = c(0.61,0.51,0.64,0.52,0.66,0.53),
  #                    conf.low = c(0.56,0.45,0.59,0.45,0.61,0.46),
  #                    conf.high = c(0.65,0.58,0.68,0.58,0.70,0.60)) %>% mutate(x = list(seq(0,6,by = 1))) %>% unnest()
  
  
  
  combined_plot_hrd = subj_data %>% 
    mutate(Accuracy = as.factor(ResponseCorrect)) %>% 
    ggplot()+
    geom_pointrange(aes(x = Alpha, y = mean, ymin = q5,ymax = q95, fill = Accuracy), col = "black", shape = 21, position = position_dodge(width = 5))+
    geom_line(data = preds,aes(x = x, y = predicted, col = Accuracy))+
    geom_ribbon(data = preds,aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = Accuracy), alpha = 0.25)+
    scale_color_manual(values = c("green","red"))+
    scale_fill_manual(values = c("green","red"))+
    facet_wrap(~facet,scales = "free")+
    theme_classic(base_size = 16)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    labs(x =  "Binned Stimulus intensity (ΔBPM)",
         y = "P(Response = faster | ΔBPM)")
  
  
  ggsave(here::here("figures","revisions","Supplementary5_HRD.tiff"), combined_plot_hrd, width = 20, height = 10, units = "in", dpi = 400)
  
  
}

#plotting all marginal parameter estimates for all models:
plot_rev = function(){
  
  
  library(tidyverse)
  library(posterior)
  line_width = 1.3
  fontsize = 30
  names = c("Threshold_placebo_vs_biso","Threshold_placebo_vs_prop","Slope_placebo_vs_biso","Slope_placebo_vs_prop","Lapse_placebo_vs_biso","Lapse_placebo_vs_prop")
  
  HRD_control <- readRDS(here::here("STAN models","HRD_control.RDS"))
  HRD_none <- readRDS(here::here("STAN models","revisions","HRD.RDS"))
  HRD_full <- readRDS(here::here("STAN models","revisions","HRD_doublecontrol.RDS"))
  HRD_final <- readRDS(here::here("STAN models","revisions","revisions_v2","HRD_final.RDS"))
  
  
  dd = rbind(
    HRD_control$draws(c("gm[2]",
                        "gm[3]",
                        "gm[6]",
                        "gm[7]",
                        "gm[10]",
                        "gm[11]")) %>% 
      as_draws_df() %>% 
      select(-contains(".")) %>% 
      rename_with(~names) %>% 
      mutate(model = "avgHR_control"),
    HRD_none$draws(c("gm[2]",
                     "gm[3]",
                     "gm[5]",
                     "gm[6]",
                     "gm[8]",
                     "gm[9]")) %>% as_draws_df() %>% 
      select(-contains(".")) %>% 
      rename_with(~names) %>% 
      mutate(model = "no_control"),
    HRD_full$draws(c("gm[2]",
                     "gm[3]",
                     "gm[7]",
                     "gm[8]",
                     "gm[12]",
                     "gm[13]")) %>% as_draws_df() %>% 
      select(-contains(".")) %>% 
      rename_with(~names)%>% 
      mutate(model = "avgHR&avgSDHR_control"),
    HRD_final$draws(c("gm[2]",
                      "gm[3]",
                      "gm[7]",
                      "gm[8]",
                      "gm[12]",
                      "gm[13]")) %>% as_draws_df() %>% 
      select(-contains(".")) %>% 
      rename_with(~names)%>% 
      mutate(model = "Final")
  ) %>% 
    pivot_longer(
      cols = -model,
      names_to = c("parameter", "comparison"),
      names_sep = "_placebo_vs_",
      values_to = "value"
    )
  
  hdi =dd %>% group_by(model,parameter,comparison) %>% summarize(mean = mean(value),
                                                                 q_95h = quantile2(value, probs = c(0.05,0.95))[[1]],
                                                                 q_5l = quantile2(value, probs = c(0.05,0.95))[[2]],
                                                                 q_80h = quantile2(value, probs = c(0.20,0.80))[[1]],
                                                                 q_20l = quantile2(value, probs = c(0.20,0.80))[[2]])
  
  
  colors = c(
    "black",
    "#66C9B2",
    "#009F73",
    "#007A58",
    "#8FE3D1",
    "#C2189B",
    "#B00089",
    "#870066",
    
    "#D66BB8"
  )
  box_color = "#F2F2F2"
  #Marginal difference for the threshold
  
  alpha = hdi %>% 
    filter(parameter == "Threshold") %>% ungroup() %>% 
    add_row(model = "",comparison = "", parameter = "Threshold", mean = NA, q_95h = NA, q_5l = NA, q_80h = NA, q_20l = NA) %>% 
    # mutate(Contrast = factor(Contrast, levels = c("Bisoprolol - Placebo", "", "Propanalol - Placebo"))) %>% 
    mutate(y_val = as.numeric(as.factor(comparison)),
           comparison = as.factor(comparison)) %>% 
    mutate(model_comparison = interaction(model, comparison)) %>% 
    mutate(model_comparison = factor(model_comparison, 
                                     levels = c(
                                       "no_control.biso", 
                                       "avgHR_control.biso", 
                                       "avgHR&avgSDHR_control.biso",
                                       "Final.biso",
                                       "no_control.prop", 
                                       "avgHR_control.prop", 
                                       "avgHR&avgSDHR_control.prop",
                                       "Final.prop"
                                     ))) %>% 
    
    mutate(y_val = ifelse(comparison == "prop",y_val+4,y_val)) %>% 
    ggplot(aes(col = interaction(model,comparison), group = model_comparison))+
    # geom_pointrange(aes(y = y_val, x = mean, xmin = q_5l, xmax = q_95h), linewidth = 1.5, show.legend = FALSE)+
    geom_pointrange(aes(y = y_val, x = mean, xmin = q_5l, xmax = q_95h), linewidth = 1.5, position = position_dodge(width = 1))+
    geom_pointrange(aes(y = y_val, x = mean, xmin = q_80h, xmax = q_20l),linewidth = 3, position = position_dodge(width = 1))+
    geom_pointrange(aes(y = y_val, x = mean, xmin = q_80h, xmax = q_20l),linewidth = 0, size = 1.2, position = position_dodge(width = 1))+
    
    facet_wrap(~parameter, scales = "free")+
    theme_classic()+
    xlab("")+
    ylab("")+
    scale_x_continuous(breaks = c(0,2.5,5.5), labels = c("0","2.5","5"))+
    coord_cartesian(xlim = c(0,6), ylim = c(0,9))+
    geom_segment(aes(x = 0, xend = 0, y = -2.3,yend = 12), linetype = 2, linewidth = line_width, col = "#1200A8", show.legend = FALSE)+
    theme(legend.position = "none",
          axis.text.y = element_blank(),
          axis.line.y = element_blank(),        # Remove the y-axis line
          strip.background = element_rect(fill = box_color, color = "white"),  # Default box style
          axis.ticks.y = element_blank())+scale_color_manual(values = colors)+  # Ensures no padding below x-axis
    theme(text = element_text(size = fontsize),           # All text
          axis.title = element_text(size = fontsize),     # Axis titles
          axis.text.x = element_text(size = fontsize),      # Axis tick labels
          legend.text = element_text(size = fontsize),    # Legend text
          legend.title = element_text(size = fontsize),   # Legend title
          plot.title = element_text(size = fontsize),      # Plot title
          axis.line.x = element_line(linewidth = line_width) 
    )
  
  alpha
  
  #for the slope
  
  beta = hdi %>% 
    filter(parameter == "Slope") %>% ungroup() %>% 
    add_row(model = "",comparison = "", parameter = "Slope", mean = NA, q_95h = NA, q_5l = NA, q_80h = NA, q_20l = NA) %>% 
    # mutate(Contrast = factor(Contrast, levels = c("Bisoprolol - Placebo", "", "Propanalol - Placebo"))) %>% 
    mutate(y_val = as.numeric(as.factor(comparison)),
           comparison = as.factor(comparison)) %>% 
    mutate(y_val = ifelse(comparison == "prop",y_val+4,y_val)) %>% 
    mutate(model_comparison = interaction(model, comparison)) %>% 
    mutate(model_comparison = factor(model_comparison, 
                                     levels = c(
                                       "no_control.biso", 
                                       "avgHR_control.biso", 
                                       "avgHR&avgSDHR_control.biso",
                                       "Final.biso",
                                       "no_control.prop", 
                                       "avgHR_control.prop", 
                                       "avgHR&avgSDHR_control.prop",
                                       "Final.prop"
                                     ))) %>% 
    ggplot(aes(col = interaction(model,comparison), group = model_comparison))+
    # geom_pointrange(aes(y = y_val, x = mean, xmin = q_5l, xmax = q_95h), linewidth = 1.5, show.legend = FALSE)+
    geom_pointrange(aes(y = y_val, x = mean, xmin = q_5l, xmax = q_95h), linewidth = 1.5, position = position_dodge(width = 1))+
    geom_pointrange(aes(y = y_val, x = mean, xmin = q_80h, xmax = q_20l),linewidth = 3, position = position_dodge(width = 1))+
    geom_pointrange(aes(y = y_val, x = mean, xmin = q_80h, xmax = q_20l),linewidth = 0, size = 1.2, position = position_dodge(width = 1))+
    geom_segment(aes(x = 0, xend = 0, y = -2.3,yend = 12), linetype = 2, linewidth = line_width, col = "#1200A8", show.legend = FALSE)+
    facet_wrap(~parameter, scales = "free")+
    theme_classic()+ theme(
      axis.title.y = element_blank(),       # Remove y-axis title
      axis.text.y = element_blank(),        # Remove y-axis text labels
      strip.background = element_rect(fill = box_color, color = "white"),  # Default box style
      axis.ticks.y = element_blank(),       # Remove y-axis ticks
      axis.line.y = element_blank(),        # Remove the y-axis line
      panel.spacing = unit(1, "lines")      # Optional: Adjust space between facets
    )+
    xlab("Differences in parameter estimates")+
    xlab("")+
    scale_x_continuous(breaks = c(-0.2,0,0.2), labels = c("-0.2","0","0.2"))+
    coord_cartesian(xlim = c(-0.25,0.25), ylim = c(0,9))+
    theme(legend.position = "none")+
    scale_color_manual(values = colors)+
    theme(
      text = element_text(size = fontsize),           # All text
      axis.title = element_text(size = fontsize),     # Axis titles
      axis.text.x = element_text(size = fontsize),      # Axis tick labels
      legend.text = element_text(size = fontsize),    # Legend text
      legend.title = element_text(size = fontsize),   # Legend title
      plot.title = element_text(size = fontsize),      # Plot title
      axis.line.x = element_line(linewidth = line_width) 
    )
  
  
  beta
  
  # for the lapse rate:
  
  lapse = hdi %>% 
    filter(parameter == "Lapse") %>% ungroup() %>% 
    add_row(model = "",comparison = "", parameter = "Lapse", mean = NA, q_95h = NA, q_5l = NA, q_80h = NA, q_20l = NA) %>% 
    # mutate(Contrast = factor(Contrast, levels = c("Bisoprolol - Placebo", "", "Propanalol - Placebo"))) %>% 
    mutate(y_val = as.numeric(as.factor(comparison)),
           comparison = as.factor(comparison)) %>% 
    mutate(model_comparison = interaction(model, comparison)) %>% 
    mutate(model_comparison = factor(model_comparison, 
                                     levels = c(
                                       "no_control.biso", 
                                       "avgHR_control.biso", 
                                       "avgHR&avgSDHR_control.biso",
                                       "Final.biso",
                                       "no_control.prop", 
                                       "avgHR_control.prop", 
                                       "avgHR&avgSDHR_control.prop",
                                       "Final.prop"
                                     ))) %>% 
    mutate(y_val = ifelse(comparison == "prop",y_val+4,y_val)) %>% 
    ggplot(aes(col = interaction(model,comparison), group = model_comparison))+
    # geom_pointrange(aes(y = y_val, x = mean, xmin = q_5l, xmax = q_95h), linewidth = 1.5, show.legend = FALSE)+
    geom_pointrange(aes(y = y_val, x = mean, xmin = q_5l, xmax = q_95h), linewidth = 1.5, position = position_dodge(width = 1))+
    geom_pointrange(aes(y = y_val, x = mean, xmin = q_80h, xmax = q_20l),linewidth = 3, position = position_dodge(width = 1))+
    geom_pointrange(aes(y = y_val, x = mean, xmin = q_80h, xmax = q_20l),linewidth = 0, size = 1.2, position = position_dodge(width = 1))+
    geom_segment(aes(x = 0, xend = 0, y = -2.3,yend = 12), linetype = 2, linewidth = line_width, col = "#1200A8", show.legend = FALSE)+
    facet_wrap(~parameter, scales = "free")+
    theme_classic()  + 
    theme(
      axis.title.y = element_blank(),       # Remove y-axis title
      axis.text.y = element_blank(),        # Remove y-axis text labels
      strip.background = element_rect(fill = box_color, color = "white"),  # Default box style
      axis.ticks.y = element_blank(),       # Remove y-axis ticks
      axis.line.y = element_blank(),        # Remove the y-axis line
      panel.spacing = unit(1, "lines")      # Optional: Adjust space between facets
    )+
    xlab("")+
    theme(legend.position = "none")+scale_color_manual(values = colors)+
    scale_x_continuous(breaks = c(-4,0,4), labels = c("-4","0","4"))+
    coord_cartesian(xlim = c(-4.5,4.5), ylim = c(0,9))+
    theme(
      text = element_text(size = fontsize),           # All text
      axis.title = element_text(size = fontsize),     # Axis titles
      axis.text.x = element_text(size = fontsize),      # Axis tick labels
      legend.text = element_text(size = fontsize),    # Legend text
      legend.title = element_text(size = fontsize),   # Legend title
      plot.title = element_text(size = fontsize),      # Plot title
      axis.line.x = element_line(linewidth = line_width) 
    )
  
  lapse
  
  library(patchwork)
  p1 = alpha|beta|lapse
  
  
  ggsave(here::here("figures","revisions","plot2_revision_v2.tiff"),p1, height = 8, width = 9, units = "in", dpi = 400)
  
  
  
  
  ## rrst
  
  RRST_nocontrol <- readRDS(here::here("STAN models","RRST.RDS"))
  RRST_doublecontrol <-readRDS(here::here("STAN models","revisions","RRST_doublecontrol.RDS")) 
  RRST_control<-readRDS(here::here("STAN models","revisions","RRST_control.RDS")) 
  RRST_final<-readRDS(here::here("STAN models","revisions","revisions_v2","RRST_final.RDS")) 
  
  
  
  library(tidyverse)
  library(posterior)
  
  names = c("Slope_placebo_vs_biso","Slope_placebo_vs_prop","Lapse_placebo_vs_biso","Lapse_placebo_vs_prop","Threshold_placebo_vs_biso","Threshold_placebo_vs_prop")
  
  dd = rbind(
    RRST_control$draws(c("gm[2]","gm[3]","gm[6]","gm[7]","gm[9]","gm[10]")) %>% as_draws_df() %>% 
      select(-contains(".")) %>% 
      rename_with(~names) %>% 
      mutate(model = "avgHR_control"),
    RRST_nocontrol$draws(c("gm[2]","gm[3]","gm[5]","gm[6]","gm[8]","gm[9]")) %>% as_draws_df() %>% 
      select(-contains(".")) %>% 
      rename_with(~names) %>% 
      mutate(model = "no_control"),
    RRST_doublecontrol$draws(c("gm[2]","gm[3]","gm[7]","gm[8]","gm[10]","gm[11]")) %>% as_draws_df() %>% 
      select(-contains(".")) %>% 
      rename_with(~names)%>% 
      mutate(model = "avgHR&avgSDHR_control"),
    RRST_final$draws(c("gm[2]","gm[3]","gm[7]","gm[8]","gm[10]","gm[11]")) %>% as_draws_df() %>% 
      select(-contains(".")) %>% 
      rename_with(~names)%>% 
      mutate(model = "Final")
  ) %>% 
    pivot_longer(
      cols = -model,
      names_to = c("parameter", "comparison"),
      names_sep = "_placebo_vs_",
      values_to = "value"
    )
  
  
  hdi =dd %>% group_by(model,parameter,comparison) %>% summarize(mean = mean(value),
                                                                 q_95h = quantile2(value, probs = c(0.05,0.95))[[1]],
                                                                 q_5l = quantile2(value, probs = c(0.05,0.95))[[2]],
                                                                 q_80h = quantile2(value, probs = c(0.20,0.80))[[1]],
                                                                 q_20l = quantile2(value, probs = c(0.20,0.80))[[2]])
  
  
  colors = c(
    "black",
    "#66C9B2",
    "#009F73",
    "#007A58",
    "#8FE3D1",
    "#C2189B",
    "#B00089",
    "#870066",
    
    "#D66BB8"
  )
  box_color = "#F2F2F2"
  #Marginal difference for the threshold
  
  
  
  alpha2 = hdi %>% 
    filter(parameter == "Threshold") %>% ungroup() %>% 
    add_row(model = "",comparison = "", parameter = "Threshold", mean = NA, q_95h = NA, q_5l = NA, q_80h = NA, q_20l = NA) %>% 
    # mutate(Contrast = factor(Contrast, levels = c("Bisoprolol - Placebo", "", "Propanalol - Placebo"))) %>% 
    mutate(y_val = as.numeric(as.factor(comparison)),
           comparison = as.factor(comparison)) %>% 
    mutate(y_val = ifelse(comparison == "prop",y_val+4,y_val)) %>% 
    mutate(model_comparison = interaction(model, comparison)) %>% 
    mutate(model_comparison = factor(model_comparison, 
                                     levels = c(
                                       "no_control.biso", 
                                       "avgHR_control.biso", 
                                       "avgHR&avgSDHR_control.biso",
                                       "Final.biso",
                                       "no_control.prop", 
                                       "avgHR_control.prop", 
                                       "avgHR&avgSDHR_control.prop",
                                       "Final.prop"
                                     ))) %>% 
    ggplot(aes(col = interaction(model,comparison), group = model_comparison))+
    # geom_pointrange(aes(y = y_val, x = mean, xmin = q_5l, xmax = q_95h), linewidth = 1.5, show.legend = FALSE)+
    geom_pointrange(aes(y = y_val, x = mean, xmin = q_5l, xmax = q_95h), linewidth = 1.5, position = position_dodge(width = 1))+
    geom_pointrange(aes(y = y_val, x = mean, xmin = q_80h, xmax = q_20l),linewidth = 3, position = position_dodge(width = 1))+
    geom_pointrange(aes(y = y_val, x = mean, xmin = q_80h, xmax = q_20l),linewidth = 0, size = 1.2, position = position_dodge(width = 1))+
    facet_wrap(~parameter, scales = "free")+
    theme_classic()+
    xlab("")+
    ylab("")+
    scale_x_continuous(breaks = c(-1,0,1), labels = c("-1","0","1"))+
    coord_cartesian(xlim = c(-1,1), ylim = c(0,9))+
    geom_segment(aes(x = 0, xend = 0, y = -2.3,yend = 12), linetype = 2, linewidth = line_width, col = "#1200A8", show.legend = FALSE)+
    theme(legend.position = "none",
          axis.text.y = element_blank(),
          axis.line.y = element_blank(),        # Remove the y-axis line
          strip.background = element_rect(fill = box_color, color = "white"),  # Default box style
          axis.ticks.y = element_blank())+scale_color_manual(values = colors)+  # Ensures no padding below x-axis
    theme(text = element_text(size = fontsize),           # All text
          axis.title = element_text(size = fontsize),     # Axis titles
          axis.text.x = element_text(size = fontsize),      # Axis tick labels
          legend.text = element_text(size = fontsize),    # Legend text
          legend.title = element_text(size = fontsize),   # Legend title
          plot.title = element_text(size = fontsize),      # Plot title
          axis.line.x = element_line(linewidth = line_width) 
    )
  
  alpha2
  
  #for the slope
  
  beta2 = hdi %>% 
    filter(parameter == "Slope") %>% ungroup() %>% 
    add_row(model = "",comparison = "", parameter = "Slope", mean = NA, q_95h = NA, q_5l = NA, q_80h = NA, q_20l = NA) %>% 
    # mutate(Contrast = factor(Contrast, levels = c("Bisoprolol - Placebo", "", "Propanalol - Placebo"))) %>% 
    mutate(y_val = as.numeric(as.factor(comparison)),
           comparison = as.factor(comparison)) %>% 
    mutate(y_val = ifelse(comparison == "prop",y_val+4,y_val)) %>% 
    mutate(model_comparison = interaction(model, comparison)) %>% 
    mutate(model_comparison = factor(model_comparison, 
                                     levels = c(
                                       "no_control.biso", 
                                       "avgHR_control.biso", 
                                       "avgHR&avgSDHR_control.biso",
                                       "Final.biso",
                                       "no_control.prop", 
                                       "avgHR_control.prop", 
                                       "avgHR&avgSDHR_control.prop",
                                       "Final.prop"
                                     ))) %>% 
    ggplot(aes(col = interaction(model,comparison), group = model_comparison))+
    # geom_pointrange(aes(y = y_val, x = mean, xmin = q_5l, xmax = q_95h), linewidth = 1.5, show.legend = FALSE)+
    geom_pointrange(aes(y = y_val, x = mean, xmin = q_5l, xmax = q_95h), linewidth = 1.5, position = position_dodge(width = 1))+
    geom_pointrange(aes(y = y_val, x = mean, xmin = q_80h, xmax = q_20l),linewidth = 3, position = position_dodge(width = 1))+
    geom_pointrange(aes(y = y_val, x = mean, xmin = q_80h, xmax = q_20l),linewidth = 0, size = 1.2, position = position_dodge(width = 1))+
    geom_segment(aes(x = 0, xend = 0, y = -2.3,yend = 12), linetype = 2, linewidth = line_width, col = "#1200A8", show.legend = FALSE)+
    facet_wrap(~parameter, scales = "free")+
    theme_classic()+ theme(
      axis.title.y = element_blank(),       # Remove y-axis title
      axis.text.y = element_blank(),        # Remove y-axis text labels
      strip.background = element_rect(fill = box_color, color = "white"),  # Default box style
      axis.ticks.y = element_blank(),       # Remove y-axis ticks
      axis.line.y = element_blank(),        # Remove the y-axis line
      panel.spacing = unit(1, "lines")      # Optional: Adjust space between facets
    )+
    xlab("Differences in parameter estimates")+
    xlab("")+
    scale_x_continuous(breaks = c(0,0.4,0.8), labels = c("0","0.4","0.8"))+
    coord_cartesian(xlim = c(0,0.8), ylim = c(0,9))+
    theme(legend.position = "none")+
    scale_color_manual(values = colors)+
    theme(
      text = element_text(size = fontsize),           # All text
      axis.title = element_text(size = fontsize),     # Axis titles
      axis.text.x = element_text(size = fontsize),      # Axis tick labels
      legend.text = element_text(size = fontsize),    # Legend text
      legend.title = element_text(size = fontsize),   # Legend title
      plot.title = element_text(size = fontsize),      # Plot title
      axis.line.x = element_line(linewidth = line_width) 
    )
  
  
  beta2
  
  # for the lapse rate:
  
  lapse2 = hdi %>% 
    filter(parameter == "Lapse") %>% ungroup() %>% 
    add_row(model = "",comparison = "", parameter = "Lapse", mean = NA, q_95h = NA, q_5l = NA, q_80h = NA, q_20l = NA) %>% 
    # mutate(Contrast = factor(Contrast, levels = c("Bisoprolol - Placebo", "", "Propanalol - Placebo"))) %>% 
    mutate(y_val = as.numeric(as.factor(comparison)),
           comparison = as.factor(comparison)) %>% 
    mutate(y_val = ifelse(comparison == "prop",y_val+4,y_val)) %>% 
    mutate(model_comparison = interaction(model, comparison)) %>% 
    mutate(model_comparison = factor(model_comparison, 
                                     levels = c(
                                       "no_control.biso", 
                                       "avgHR_control.biso", 
                                       "avgHR&avgSDHR_control.biso",
                                       "Final.biso",
                                       "no_control.prop", 
                                       "avgHR_control.prop", 
                                       "avgHR&avgSDHR_control.prop",
                                       "Final.prop"
                                     ))) %>% 
    ggplot(aes(col = interaction(model,comparison), group = model_comparison))+
    # geom_pointrange(aes(y = y_val, x = mean, xmin = q_5l, xmax = q_95h), linewidth = 1.5, show.legend = FALSE)+
    geom_pointrange(aes(y = y_val, x = mean, xmin = q_5l, xmax = q_95h), linewidth = 1.5, position = position_dodge(width = 1))+
    geom_pointrange(aes(y = y_val, x = mean, xmin = q_80h, xmax = q_20l),linewidth = 3, position = position_dodge(width = 1))+
    geom_pointrange(aes(y = y_val, x = mean, xmin = q_80h, xmax = q_20l),linewidth = 0, size = 1.2, position = position_dodge(width = 1))+
    geom_segment(aes(x = 0, xend = 0, y = -2.3,yend = 12), linetype = 2, linewidth = line_width, col = "#1200A8", show.legend = FALSE)+
    facet_wrap(~parameter, scales = "free")+
    theme_classic()  + 
    theme(
      axis.title.y = element_blank(),       # Remove y-axis title
      axis.text.y = element_blank(),        # Remove y-axis text labels
      strip.background = element_rect(fill = box_color, color = "white"),  # Default box style
      axis.ticks.y = element_blank(),       # Remove y-axis ticks
      axis.line.y = element_blank(),        # Remove the y-axis line
      panel.spacing = unit(1, "lines")      # Optional: Adjust space between facets
    )+
    xlab("")+
    theme(legend.position = "none")+scale_color_manual(values = colors)+
    scale_x_continuous(breaks = c(-2,0,1), labels = c("-2","0","1"))+
    coord_cartesian(xlim = c(-2.5,2), ylim = c(0,9))+
    theme(
      text = element_text(size = fontsize),           # All text
      axis.title = element_text(size = fontsize),     # Axis titles
      axis.text.x = element_text(size = fontsize),      # Axis tick labels
      legend.text = element_text(size = fontsize),    # Legend text
      legend.title = element_text(size = fontsize),   # Legend title
      plot.title = element_text(size = fontsize),      # Plot title
      axis.line.x = element_line(linewidth = line_width) 
    )
  
  lapse2
  
  library(patchwork)
  p2 = alpha2|beta2|lapse2
  
  ggsave(here::here("figures","revisions","plot2_revisionRRST_v2.tiff"),p1, height = 8, width = 9, units = "in", dpi = 400)
  
  # legend:
  
  legend = hdi %>% 
    filter(parameter == "Lapse") %>% ungroup() %>% 
    add_row(model = "",comparison = "", parameter = "Lapse", mean = NA, q_95h = NA, q_5l = NA, q_80h = NA, q_20l = NA) %>% 
    # mutate(Contrast = factor(Contrast, levels = c("Bisoprolol - Placebo", "", "Propanalol - Placebo"))) %>% 
    mutate(y_val = as.numeric(as.factor(comparison)),
           comparison = as.factor(comparison)) %>% 
    mutate(y_val = ifelse(comparison == "prop",y_val+4,y_val)) %>% 
    mutate(model_comparison = interaction(model, comparison)) %>% 
    mutate(model_comparison = factor(model_comparison, 
                                     levels = c(
                                       "no_control.biso", 
                                       "avgHR_control.biso", 
                                       "avgHR&avgSDHR_control.biso",
                                       "Final.biso",
                                       "no_control.prop", 
                                       "avgHR_control.prop", 
                                       "avgHR&avgSDHR_control.prop",
                                       "Final.prop"
                                     ))) %>% 
    ggplot(aes(col = interaction(model,comparison), group = model_comparison))+
    # geom_pointrange(aes(y = y_val, x = mean, xmin = q_5l, xmax = q_95h), linewidth = 1.5, show.legend = FALSE)+
    geom_pointrange(aes(y = y_val, x = mean, xmin = q_5l, xmax = q_95h), linewidth = 1.5, position = position_dodge(width = 1))+
    geom_pointrange(aes(y = y_val, x = mean, xmin = q_80h, xmax = q_20l),linewidth = 3, position = position_dodge(width = 1))+
    geom_pointrange(aes(y = y_val, x = mean, xmin = q_80h, xmax = q_20l),linewidth = 0, size = 1.2, position = position_dodge(width = 1))+
    geom_segment(aes(x = 0, xend = 0, y = -2.3,yend = 12), linetype = 2, linewidth = line_width, col = "#1200A8", show.legend = FALSE)+
    facet_wrap(~parameter, scales = "free")+
    theme_classic()  + 
    theme(
      axis.title.y = element_blank(),       # Remove y-axis title
      axis.text.y = element_blank(),        # Remove y-axis text labels
      strip.background = element_rect(fill = box_color, color = "white"),  # Default box style
      axis.ticks.y = element_blank(),       # Remove y-axis ticks
      axis.line.y = element_blank(),        # Remove the y-axis line
      panel.spacing = unit(1, "lines")      # Optional: Adjust space between facets
    )+
    xlab("")+
    scale_color_manual(values = colors)+
    scale_x_continuous(breaks = c(-2,0,1), labels = c("-2","0","1"))+
    coord_cartesian(xlim = c(-2.5,2), ylim = c(0,9))+
    theme(
      text = element_text(size = fontsize),           # All text
      axis.title = element_text(size = fontsize),     # Axis titles
      axis.text.x = element_text(size = fontsize),      # Axis tick labels
      legend.text = element_text(size = fontsize),    # Legend text
      legend.title = element_text(size = fontsize),   # Legend title
      plot.title = element_text(size = fontsize),      # Plot title
      axis.line.x = element_line(linewidth = line_width) 
    )
  
  legend
  ggsave(here::here("figures","revisions","plot2_revision_revisionsv2_legend.tiff"),legend, height = 8, width = 9, units = "in", dpi = 400)
  
  marginaleffects = alpha | beta | lapse | alpha2 | beta2 | lapse2
  
  ggsave(here::here("figures","revisions","plot2_revision_revisionsv2_marginaleffects.tiff"),marginaleffects, height = 8, width = 16, units = "in", dpi = 100)
  
  
}



HRV_scatter <- function(data_to_plot, Drug, metric_label){
  custom_colors <- c("placebo" = "#1200A8", "propranolol" = "#B00089", "bisoprolol" = "#009F73")
  level_order <- c('placebo', 'propranolol', 'bisoprolol')
  
  plotted_data <- ggplot()+
    geom_point(aes(x = factor(Drug, level = level_order), y = data_to_plot, fill = Drug), alpha = 0.4, width = 2, size = 5, col = "black", shape = 21,
               position = position_jitterdodge(seed = 123,jitter.width = 0.5, dodge.width = 0.75)) +
    geom_boxplot(aes(x = factor(Drug, level = level_order), y = data_to_plot, col = Drug), alpha = 0.80, size = 1.2,outliers = F)+
    scale_fill_manual(values = custom_colors) +  # Custom colors
    scale_color_manual(values = custom_colors)+
    labs(title = " ",
         x = " ",
         y = metric_label) +
    theme_classic() +
    theme(legend.position = "none", 
          plot.title = element_text(size = 18, face = "bold"),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18), # Center the title and refine axis labels
          axis.text = element_text(size = 18))
  
  return(plotted_data)
}