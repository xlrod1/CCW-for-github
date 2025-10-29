library(dplyr)

# Function: clone_censor
# Description: This function generates four clones of the dataset with different censoring rules.
# It applies different rules for defining follow-up outcome, censoring, and event occurrence.
#
# Parameters:
# - data: The input dataset.
# - treatment_threshold: The threshold for early treatment (default is 12 hours).
#
# Steps:
# 1. Creates separate clones (1 to 4) with specific eligibility criteria.
# 2. Assigns `fup_outcome`, `censor`, and `outcome` based on conditions for each clone.
# 3. Combines the cloned datasets into one final output.

apply_censoring_clone1 <- function(data, treatment_threshold) {
  data %>%
    mutate(
      clone_id = 1,
      fup_outcome = case_when(
        is.na(time.prom.oxitocin.hour) & time.prom.to.delivery <= treatment_threshold ~ time.prom.to.delivery,
        is.na(time.prom.oxitocin.hour) & time.prom.to.delivery > treatment_threshold ~ treatment_threshold,
        time.prom.oxitocin.hour <= treatment_threshold & first_time_reach_c_o <= time.prom.oxitocin.hour ~ time.prom.to.delivery,
        time.prom.oxitocin.hour <= treatment_threshold & first_time_reach_c_o > time.prom.oxitocin.hour ~ time.prom.oxitocin.hour,
        time.prom.oxitocin.hour > treatment_threshold ~ treatment_threshold,
        TRUE ~ NA_real_
      ),
      censor = as.integer(fup_outcome != time.prom.to.delivery),
      outcome = as.integer(fup_outcome == time.prom.to.delivery)
    )
}

apply_censoring_clone2 <- function(data, treatment_threshold) {
  data %>%
    mutate(
      clone_id = 2,
      fup_outcome = case_when(
        is.na(time.prom.oxitocin.hour) & time.prom.to.delivery <= treatment_threshold & first_time_reach_c_o <= treatment_threshold ~ first_time_reach_c_o,
        is.na(time.prom.oxitocin.hour) & time.prom.to.delivery > treatment_threshold & first_time_reach_c_o > treatment_threshold ~ treatment_threshold,
        is.na(time.prom.oxitocin.hour) & time.prom.to.delivery > treatment_threshold & first_time_reach_c_o <= treatment_threshold ~ first_time_reach_c_o,
        time.prom.oxitocin.hour <= treatment_threshold & first_time_reach_c_o <= time.prom.oxitocin.hour ~ first_time_reach_c_o,
        time.prom.oxitocin.hour <= treatment_threshold & first_time_reach_c_o > time.prom.oxitocin.hour ~ time.prom.to.delivery,
        time.prom.oxitocin.hour > treatment_threshold & first_time_reach_c_o <= treatment_threshold ~ first_time_reach_c_o,
        time.prom.oxitocin.hour > treatment_threshold & first_time_reach_c_o > treatment_threshold ~ treatment_threshold,
        TRUE ~ NA_real_
      ),
      censor = as.integer(fup_outcome != time.prom.to.delivery),
      outcome = as.integer(fup_outcome == time.prom.to.delivery)
    )
}

apply_censoring_clone3 <- function(data, treatment_threshold) {
  data %>%
    mutate(
      clone_id = 3,
      fup_outcome = case_when(
        is.na(time.prom.oxitocin.hour) ~ time.prom.to.delivery,
        time.prom.oxitocin.hour <= treatment_threshold ~ time.prom.oxitocin.hour,
        time.prom.oxitocin.hour > treatment_threshold & first_time_reach_c_o <= time.prom.oxitocin.hour ~ time.prom.to.delivery,
        time.prom.oxitocin.hour > treatment_threshold & first_time_reach_c_o > time.prom.oxitocin.hour ~ time.prom.oxitocin.hour,
        TRUE ~ NA_real_
      ),
      censor = as.integer(fup_outcome != time.prom.to.delivery),
      outcome = as.integer(fup_outcome == time.prom.to.delivery)
    )
}

apply_censoring_clone4 <- function(data, treatment_threshold) {
  data %>%
    mutate(
      clone_id = 4,
      fup_outcome = case_when(
        is.na(time.prom.oxitocin.hour) ~ first_time_reach_c_o,
        time.prom.oxitocin.hour <= treatment_threshold & first_time_reach_c_o <= time.prom.oxitocin.hour ~ first_time_reach_c_o,
        time.prom.oxitocin.hour <= treatment_threshold & first_time_reach_c_o > time.prom.oxitocin.hour ~ time.prom.oxitocin.hour,
        time.prom.oxitocin.hour > treatment_threshold & first_time_reach_c_o > time.prom.oxitocin.hour ~ time.prom.to.delivery,
        time.prom.oxitocin.hour > treatment_threshold & first_time_reach_c_o <= time.prom.oxitocin.hour ~ first_time_reach_c_o,
        TRUE ~ NA_real_
      ),
      censor = as.integer(fup_outcome != time.prom.to.delivery),
      outcome = as.integer(fup_outcome == time.prom.to.delivery)
    )
}

clone_censor <- function(data, treatment_threshold = 12) {
  bind_rows(
    apply_censoring_clone1(data, treatment_threshold),
    apply_censoring_clone2(data, treatment_threshold),
    apply_censoring_clone3(data, treatment_threshold),
    apply_censoring_clone4(data, treatment_threshold)
  )
}
