# Debug subject mapping
library(dplyr)
library(readr)

df_trial <- read_csv("data/HRD_trial_level_data.csv", show_col_types = FALSE)
df_hrv <- read_csv("data/alldrugs_full_HRV.csv", show_col_types = FALSE)

cat("Unique subjects in Trial Data:\n")
print(unique(df_trial$subject))

cat("\nUnique SubNo in HRV Data:\n")
print(unique(df_hrv$SubNo))

cat("\nUnique visits in Trial Data:\n")
print(unique(df_trial$visit))

# Test my mapping logic
df_hrv_mapped <- df_hrv %>%
  mutate(
    subject_mapped = sprintf("sub_3%03d", as.numeric(SubNo))
  )

cat("\nMapped subjects from HRV Data:\n")
print(unique(df_hrv_mapped$subject_mapped))

# Check intersection
common_subjects <- intersect(unique(df_trial$subject), unique(df_hrv_mapped$subject_mapped))
cat("\nCommon subjects:\n")
print(common_subjects)

cat("\nMissing subjects (in Trial but not in HRV):\n")
print(setdiff(unique(df_trial$subject), unique(df_hrv_mapped$subject_mapped)))
