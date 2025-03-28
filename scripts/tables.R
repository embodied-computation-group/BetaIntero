# function for generating the demographics_table in the supplementary material

Demographics_table = function(){
  
  #read the data csv
  TraitSurveys_scores <- read.csv(here::here("data","TraitSurveys_scores.csv"))
  
  TraitSurveys_scores <- TraitSurveys_scores %>% 
    mutate(sex = str_trim(sex)) %>% 
    filter(!Subject_ID %in% c(4024,4049))
  
  # Summary of numeric variables (age, height, weight, BMI, ADHD, Anxiety, Depression)
  numeric_summary <- tibble(
    `Variable` = c("Age (Mean ± SD)", 
                   "Height (cm)", 
                   "Weight (kg)", 
                   "BMI", 
                   "ADHD Score", 
                   "Anxiety Score (BAI)", 
                   "Depression Score (BDI II)"),
    `Mean (± SD)` = c(
      paste0(round(mean(TraitSurveys_scores$age, na.rm = TRUE), 1), " ± ", round(sd(TraitSurveys_scores$age, na.rm = TRUE), 1)),
      paste0(round(mean(TraitSurveys_scores$height, na.rm = TRUE), 1), " ± ", round(sd(TraitSurveys_scores$height, na.rm = TRUE), 1)),
      paste0(round(mean(TraitSurveys_scores$weight, na.rm = TRUE), 1), " ± ", round(sd(TraitSurveys_scores$weight, na.rm = TRUE), 1)),
      paste0(round(mean(TraitSurveys_scores$bmi_calc, na.rm = TRUE), 1), " ± ", round(sd(TraitSurveys_scores$bmi_calc, na.rm = TRUE), 1)),
      paste0(round(mean(TraitSurveys_scores$ADHD_Main, na.rm = TRUE), 1), " ± ", round(sd(TraitSurveys_scores$ADHD_Main, na.rm = TRUE), 1)),
      paste0(round(mean(TraitSurveys_scores$BAI_Anxiety, na.rm = TRUE), 1), " ± ", round(sd(TraitSurveys_scores$BAI_Anxiety, na.rm = TRUE), 1)),
      paste0(round(mean(TraitSurveys_scores$BDI_II_Depression, na.rm = TRUE), 1), " ± ", round(sd(TraitSurveys_scores$BDI_II_Depression, na.rm = TRUE), 1))
    )
  )
  
  # Summary of categorical variable (Sex)
  sex_summary <- TraitSurveys_scores %>%
    group_by(sex) %>%
    summarize(N = n()) %>%
    mutate(Percent = round(100 * N / sum(N), 1)) %>%
    summarize(Sex = paste0(N[1], " (", Percent[1], "%) Female, ", N[2], " (", Percent[2], "%) Male"))
  
  # Combine the numeric and categorical summaries into a single tibble
  demographic_table <- bind_rows(
    numeric_summary,
    tibble(
      Variable = "Sex",
      `Mean (± SD)` = sex_summary$Sex
    )
  )
  
  # Create flextable for APA formatting
  flex_table <- flextable(demographic_table) %>%
    add_header_row(values = c("Demographics", "Mean (± SD)"), colwidths = c(1, 1)) %>%
    bold(i = 1, bold = TRUE) %>%
    autofit() %>%
    set_table_properties(layout = "autofit")
  
  # Save the table to a Word document
  doc <- read_docx() %>%
    body_add_flextable(flex_table) %>%
    body_add_par("Note: Standard deviations are presented in parentheses. This table follows APA formatting guidelines.", style = "Normal") %>%
    print(target = here::here("docs/Demographic_Table_APA.docx"))
  
}

# function for generating the IBQ in the supplementary material

IBQ_table = function(){
  
  # Data for the table with Domain/Modality
  item_data <- data.frame(
    `Domain/Modality` = c("Breathing Sensitivity", "Heartbeat Sensitivity", "Heartbeat Affective", "Breathing Affective"),
    `Item Text` = c(
      "I am generally sensitive to the sensations of my own breathing.",
      "I can generally tell how fast my heart is beating.",
      "When you focus on the sensation of your heartbeat, how does it generally make you feel?",
      "When you focus on the sensation of your breathing, how does it generally make you feel?"
    ),
    `Response Scale` = c(
      "1 (Strongly Disagree) - 10 (Strongly Agree)",
      "1 (Strongly Disagree) - 10 (Strongly Agree)",
      "1 (Extremely Disturbing) - 7 (Extremely Relaxing)",
      "1 (Extremely Disturbing) - 7 (Extremely Relaxing)"
    )
  )
  
  # Create a flextable
  ft <- flextable(item_data) %>%
    add_header_row(values = c("Interoceptive Belief Questionnaire: Domains/Modalities and Response Scales"), colwidths = 3) %>%
    bold(part = "header") %>%
    autofit() %>%
    set_table_properties(layout = "autofit")
  
  # Save the table to a Word document with the note below the table
  doc <- read_docx() %>%
    body_add_flextable(ft) %>%
    body_add_par("Note: Items 1 and 2 are belief-based questions rated on a 1 to 10 scale. Items 3 and 4 assess affective responses to interoceptive focus, rated on a 1 to 7 scale.", style = "Normal") %>%
    print(target = here::here("docs/Interoceptive_Belief_Questionnaire_Corrected.docx"))
  
}

# function for generating the table that hold information on the priors of the stanmodels

prior_tables = function(){
  
  # for the HRDT
  get_hrd_table_priors = function(){
    library(flextable)
    
    # Data for the flextable
    prior_data <- data.frame(
      parameter = c(
        "gm_Placebo_threshold", "gm_Biso_dif_threshold", "gm_Prop_dif_threshold", "gm_scaled_BPM_threshold", "gm_Placebo_slope", "gm_Biso_dif_slope", "gm_Prop_dif_slope", "gm_scaled_BPM_slope",
        "gm_Placebo_lapse", "gm_Biso_dif_lapse", "gm_Prop_dif_lapse", 
        "gsd_Placebo_threshold", "gsd_Biso_dif_threshold", "gsd_Prop_dif_threshold", "gsd_scaled_BPM_threshold", "gsd_Placebo_slope", "gsd_Biso_dif_slope", "gsd_Prop_dif_slope", "gsd_scaled_BPM_slope",
        "gsd_Placebo_lapse", "gsd_Biso_dif_lapse", "gsd_Prop_dif_lapse",  "R"
      ),
      mean = c(
        -10, 0, 0, 0, 2.25, 0, 0, 0,
        -5.5, 0, 0, 
        10, 0, 0, 0, 0.25, 0, 
        0, 0, 2.5, 0, 0, NA
      ),
      sd = c(
        5, 5, 5, 5, 0.5, 3, 3, 3, 
        1.5, 2, 2, 
        10, 10, 10, 3, 0.5, 2, 
        2, 2, 1, 3, 3, NA
      ),
      eta = c(
        rep(NA, 22), 2 # Only L_u has eta = 2
      ),
      distribution = c(
        rep("Normal", 11),
        rep("Normal+", 11),
        "LKJ"
      )
    ) %>% mutate(group = case_when(
      startsWith(parameter, "gm_") ~ "mean",
      startsWith(parameter, "gsd_") ~ "sd",
      parameter == "R" ~ "correlation",
      TRUE ~ NA_character_
    ),
    parameter = gsub("^gm_|^gsd_", "", parameter)) %>% 
      select(group,parameter,mean,sd,eta,distribution)
    
    
    # Create the flextable
    ft <- flextable(prior_data)
    ft <- set_header_labels(
      ft, parameter = "Parameter", mean = "Mean", sd = "SD",
      eta = "Eta", distribution = "Distribution"
    )
    ft <- fontsize(ft, size = 10, part = "all")
    ft <- autofit(ft)
    
    # Display the flextable
    ft
  } 
  
  # for the RRST
  get_rrst_table_priors = function(){
    library(flextable)
    
    # Data for the flextable
    prior_data <- data.frame(
      parameter = c(
        "gm_Placebo_threshold", "gm_Biso_dif_threshold", "gm_Prop_dif_threshold", "gm_Placebo_slope", "gm_Biso_dif_slope", "gm_Prop_dif_slope",
        "gm_Placebo_lapse", "gm_Biso_dif_lapse", "gm_Prop_dif_lapse", 
        "gsd_Placebo_threshold", "gsd_Biso_dif_threshold", "gsd_Prop_dif_threshold", "gsd_Placebo_slope", "gsd_Biso_dif_slope", "gsd_Prop_dif_slope",
        "gsd_Placebo_lapse", "gsd_Biso_dif_lapse", "gsd_Prop_dif_lapse",  "R"
      ),
      mean = c(
        10, 0, 0, 0, 0, 0,
        -4, 0, 0, 
        0, 0, 0, 0, 0, 
        0, 0, 0, 0, NA
      ),
      sd = c(
        10, 5, 5, 3, 3, 3, 3, 
        1.5, 2,
        5, 5, 5, 3, 3, 3, 
        3, 3, 3, NA
      ),
      eta = c(
        rep(NA, 18), 2 # Only L_u has eta = 2
      ),
      distribution = c(
        rep("Normal", 9),
        rep("Normal+", 9),
        "LKJ"
      )
    ) %>% mutate(group = case_when(
      startsWith(parameter, "gm_") ~ "mean",
      startsWith(parameter, "gsd_") ~ "sd",
      parameter == "R" ~ "correlation",
      TRUE ~ NA_character_
    ),
    parameter = gsub("^gm_|^gsd_", "", parameter)) %>% 
      select(group,parameter,mean,sd,eta,distribution)

    
    # Create the flextable
    ft <- flextable(prior_data)
    ft <- set_header_labels(
      ft, parameter = "Parameter", mean = "Mean", sd = "SD",
      eta = "Eta", distribution = "Distribution"
    )
    ft <- fontsize(ft, size = 10, part = "all")
    ft <- autofit(ft)
    
    # Display the flextable
    ft
    
  }
  
  # make both and return
  hrd = get_hrd_table_priors()
  rrst = get_rrst_table_priors()
  
  return(list(hrd,rrst))
}
