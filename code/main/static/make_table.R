library(dplyr)
library(tidyr)

# Mapping clone names to readable labels
clone_map <- c(
  "window_0_8" = "a_1",
  "window_8_16" = "a_2",
  "window_16_24" = "a_3",
  "window_24_32" = "a_4",
  "window_32_Inf" = "a_5"
)

#table for cox
#lead the data
cif_plot_df_TV<-read.csv( "G:\\My Drive\\Phd4\\R code\\final_code\\SW\\final\\results\\main\\cif_TV.csv")
cif_plot_df_not_TV<-read.csv( "G:\\My Drive\\Phd4\\R code\\final_code\\SW\\final\\results\\main\\cif_no_TV.csv")

cif_plot_df_TV<-cif_plot_df_TV %>% 
  mutate(TV="yes")

cif_plot_df_not_TV<-cif_plot_df_not_TV %>% 
  mutate(TV="no") %>% 
  filter(method!="No weights")


cif_plot_df<-rbind(cif_plot_df_TV,cif_plot_df_not_TV)
highlight_pts <- cif_plot_df %>% filter(time %in% c(12, 24, 48,72))

h.cox.SW<-highlight_pts %>% filter(TV=="yes",
                                   method %in% c( "Cox.SW‐model weights")) %>% 
  arrange(event,clone,time) 

h.log.SW<-highlight_pts %>% filter(TV=="yes",
                                   method %in% c( "Logistic.SW‐model weights")) %>% 
  arrange(event,clone,time) 

h.not.SW<-highlight_pts %>% filter(TV=="yes",
                                   method %in% c( "No weights")) %>% 
  arrange(event,clone,time) 
h.log.SW<-highlight_pts %>% filter(TV=="no",
                                   method %in% c( "Logistic.SW‐model weights")) %>% 
  arrange(event,clone,time) 

h.cox.SW<-highlight_pts %>% filter(TV=="no",
                                   method %in% c( "Cox.SW‐model weights")) %>% 
  arrange(event,clone,time) 


df_summary <- h.log.SW%>%
  # Keep only time points of interest
  filter(time %in% c(12, 24, 48, 72)) %>%
  
  # Recode clone to strategy label
  mutate(
    strategy = recode(clone, !!!clone_map),
    
    # Create formatted summary string
    summary_str = sprintf("%.2f (%.2f, %.2f)", prob, cif_lower, cif_upper),
    
    # Create column identifier like: "Assisted_12"
    col_name = paste0(event, "_", time)
  ) %>%
  
  # Pivot wider: one row per strategy
  select(strategy, col_name, summary_str) %>%
  pivot_wider(
    names_from = col_name,
    values_from = summary_str
  ) %>%
  arrange(strategy)
# Optional: arrange column order (e.g., a_1_12, a_2_24, etc.)
# df_summary <- df_summary %>% select(event, sort(tidyselect::peek_vars()))
write.csv(df_summary,
          file = "G:/My Drive/Phd4/R code/final_code/SW/final/results/main/tables/log_not_tv.csv",
          row.names = FALSE)
