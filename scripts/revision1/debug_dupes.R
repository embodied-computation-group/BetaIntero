# Check for duplicates in HRV data
library(dplyr)
library(readr)

df_hrv <- read_csv("data/alldrugs_full_HRV.csv", show_col_types = FALSE)

# Check if SubNo + Drug is unique
dupes <- df_hrv %>%
  group_by(SubNo, Drug) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n > 1)

print(dupes)

# Check specific case 0024
df_0024 <- df_hrv %>% filter(SubNo == "0024")
print(df_0024)
