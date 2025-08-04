library(dplyr)
library(purrr)
library(data.table)
library(cmprsk)
library(tidyr)
library(ggplot2)
source("G:\\My Drive\\Phd4\\R code\\final_code\\f_compute_survival_for_clone_new_event_time.R")
df.only.base<- read.csv("G:/My Drive/Phd4/data/final_data/data_base_line.csv") 
measur_df<-fread("G:/My Drive/Phd4/data/final_data/opening_height.csv") 
measur_df<-measur_df %>% 
  mutate(opening_clean_cf=ifelse(time.prom.measurement==0&is.na(opening_clean_cf),0,opening_clean_cf),
         height_clean_cf=ifelse(time.prom.measurement==0&is.na(height_clean_cf),-3,height_clean_cf)) %>% 
  filter(!is.na(height_clean_cf),!is.na(opening_clean_cf))
                                      
#select only the women with measurmnt information
df.only.base<-df.only.base %>% right_join(.,measur_df %>% select(id) %>% distinct())

#Opening cut-off X=3 canchange here
c_o<-4
first.time.open<-measur_df%>%
  filter(opening_clean_cf >=c_o) %>%
  group_by(id) %>%
  summarize(first_time_reach_c_o = min(time.prom.measurement)) %>%
  ungroup()
#left_join with df.only.base
df.only.base<- df.only.base %>% left_join(first.time.open)


round(colMeans(is.na(df.only.base)) * 100, 1)

#for now we will do complete casews anayslsis-without GBS+proteinuria
#precenteg of NA
round(colMeans(is.na(df.only.base)) * 100, 1)

list.of.confounder<-age ~ bmi+parity+ART_pregnancy+
  pre_HTN+Gestational_diabetes+Gestational_HTN+
anemia+thrombocytopenia
  

#table(df.only.base$parity)
vars <- all.vars(list.of.confounder)
df.only.base<-df.only.base %>%
  # select only the vars in the formula, then drop any rows with NA
  filter(
    if_all(all_of(vars), ~ !is.na(.))
  )

brks<- c(0,12,Inf)
brks.e <- c(0,12)

df.only.base <- df.only.base %>%
  mutate(
    woman_type = case_when(
      first_time_reach_c_o <= 12 & time.prom.oxitocin.hour <= 12 ~ "A",
      first_time_reach_c_o <= 12 & (time.prom.oxitocin.hour > 12|is.na(time.prom.oxitocin.hour) ) ~ "B",
      first_time_reach_c_o > 12  & time.prom.oxitocin.hour <= 12 ~ "C",
      first_time_reach_c_o > 12  & (time.prom.oxitocin.hour > 12|is.na(time.prom.oxitocin.hour) )  ~ "D",
      TRUE ~ NA_character_
    )
  )


#Here for each women I create the variable fup_outcome- which the the 
#time i follow-up the women util censoring
make_clone_strategy <- function(data, strategy) {
  data %>%
    mutate(
      clone = paste0("strategy_", strategy),
      fup_outcome = time.prom.to.delivery,
      
      # STRATEGY 1: Induce ≤ 12h if dilation < 4 cm
      fup_outcome = if_else(strategy == 1 & !is.na(first_time_reach_c_o) & first_time_reach_c_o <= 12,
                            first_time_reach_c_o, fup_outcome),
      fup_outcome = if_else(strategy == 1 & is.na(first_time_reach_c_o) & 
                              is.na(time.prom.oxitocin.hour) & time.prom.to.delivery > 12,
                            12, fup_outcome),
      fup_outcome = if_else(strategy == 1 & !is.na(time.prom.oxitocin.hour) & time.prom.oxitocin.hour > 12,
                            time.prom.oxitocin.hour, fup_outcome),
      
      # STRATEGY 2: Induce >12h if dilation < 4 cm
      fup_outcome = if_else(strategy == 2 & !is.na(first_time_reach_c_o) & first_time_reach_c_o <= 12,
                            first_time_reach_c_o, fup_outcome),
      fup_outcome = if_else(strategy == 2 & !is.na(time.prom.oxitocin.hour) & time.prom.oxitocin.hour <= 12,
                            time.prom.oxitocin.hour, fup_outcome),
      fup_outcome = if_else(strategy == 2 & is.na(time.prom.oxitocin.hour) & 
                              time.prom.to.delivery <= 12,
                            pmin(first_time_reach_c_o, time.prom.to.delivery, na.rm = TRUE),
                            fup_outcome),
      fup_outcome = if_else(strategy == 2 & is.na(time.prom.oxitocin.hour) &
                              !is.na(first_time_reach_c_o),
                            first_time_reach_c_o, fup_outcome),
      
      # STRATEGY 3: Induce ≤12h if dilation ≥ 4 cm
      fup_outcome = if_else(strategy == 3 & (is.na(first_time_reach_c_o) | first_time_reach_c_o > 12),
                            12, fup_outcome),
      fup_outcome = if_else(strategy == 3 & !is.na(first_time_reach_c_o) & first_time_reach_c_o <= 12 &
                              (is.na(time.prom.oxitocin.hour) | time.prom.oxitocin.hour > 12),
                            12, fup_outcome),
      fup_outcome = if_else(strategy == 3 & !is.na(time.prom.oxitocin.hour) & 
                              time.prom.oxitocin.hour < first_time_reach_c_o,
                            time.prom.oxitocin.hour, fup_outcome),
      
      # STRATEGY 4: Induce >12h if dilation ≥ 4 cm
      fup_outcome = if_else(strategy == 4 & !is.na(time.prom.oxitocin.hour) & time.prom.oxitocin.hour <= 12,
                            time.prom.oxitocin.hour, fup_outcome),
      fup_outcome = if_else(strategy == 4 & !is.na(time.prom.oxitocin.hour) & !is.na(first_time_reach_c_o) &
                              time.prom.oxitocin.hour > 12 & time.prom.oxitocin.hour < first_time_reach_c_o,
                            time.prom.oxitocin.hour, fup_outcome),
      fup_outcome = if_else(strategy == 4 & is.na(time.prom.oxitocin.hour) & time.prom.to.delivery <= 12,
                            pmin(first_time_reach_c_o, time.prom.to.delivery, na.rm = TRUE),
                            fup_outcome),
      
      outcome = as.integer(fup_outcome == time.prom.to.delivery),
      censor  = as.integer(fup_outcome <  time.prom.to.delivery)
    )
}


# build all the clones in one data‐frame

df_clones <- map_dfr(1:4, ~ make_clone_strategy(df.only.base, .x))

epsilon <- 1e-4
df_clones<- df_clones %>%
  mutate(fup_outcome = if_else(fup_outcome == 0, epsilon, fup_outcome))
#add the variabel foe competenf risks:
df_clones <- df_clones%>%
  mutate(
    status = factor(
      case_when(
        censor == 1                   ~ 0L,
        outcome == 1 & NVD == 1       ~ 1L,
        outcome == 1 & NVD == 0       ~ 2L
      ),
      levels = 0:2,
      labels = c("censored","vaginal","other")
    )
  )

#The denuminator formulas

glm_formula<- (censor == 0) ~ ns(t, df =2) +opening_clean_cf+height_clean_cf +age +bmi+parity+ART_pregnancy+
  pre_HTN+Gestational_diabetes+Gestational_HTN+
  anemia+thrombocytopenia
cox_formula <- Surv(t_start, t_stop, censor)~opening_clean_cf +height_clean_cf+age +bmi+parity+ART_pregnancy+
  pre_HTN+Gestational_diabetes+Gestational_HTN+
  anemia+thrombocytopenia

strategy_ids <- 1:4
results_list <- vector("list", length(strategy_ids))
#i=2
for (i in seq_along(strategy_ids)) {
  strategy <- strategy_ids[i]
  df_sub   <- df_clones %>% filter(clone == paste0("strategy_", strategy))
  
  data = df_sub
  only.measur.full.time = measur_df
  glm_censor_formula = glm_formula
  cox_censor_formula = cox_formula
  time_grid = seq(0, 72, by = 0.5)
  competing.risks = TRUE
  
  
  surv_dat <- compute_survival_for_clone_dynamic(
    data = df_sub,
    only.measur.full.time = measur_df,
    glm_censor_formula = glm_formula,
    cox_censor_formula = cox_formula,
    time_grid = seq(0, 72, by = 0.5),
    competing.risks = TRUE
  )
  
  surv_dat$strategy <- paste0("strategy_", strategy)
  results_list[[i]] <- surv_dat
  print(i)
}

# row-bind them all into one data frame
cif_by_clone <- bind_rows(results_list)
# now surv_by_clone has columns: time, logit_w, cox_w, unweighted, and clone

# h1<-surv_by_clone%>% filter(clone=="window_0_12") %>% select(time,unw_risk)
# h2<-surv_by_clone%>% filter(clone!="window_0_12") %>% select(unw_risk)
# hh<-cbind(h1,h2)
# Now I can make KM plot for all kind of delivery
# assume your data.frame is called `surv_tbl`
# with columns: time, logit_w, cox_w, unweighted, clone

cif_clean <- cif_by_clone %>%
  rename(
    coxSW_vaginal     = cox_vaginal_SW,
    coxSW_competing   = cox_competing_SW,
    logSW_vaginal     = log_vaginal_SW,
    logSW_competing   = log_competing_SW
  )

cif_plot_df <- cif_clean %>%
  pivot_longer(
    cols = -c(time, clone),
    names_to = c("method", "event"),
    names_sep = "_",
    values_to = "cif"
  ) %>%
  mutate(
    method = recode(method,
                    "unw"      = "No weights",
                    "cox"      = "Cox‐model weights",
                    "cox_SW"   = "Cox‐model S weights",
                    "log_SW"   = "Logistic‐model S weights"
    ),
    event = recode(event,
                   "vaginal"   = "Vaginal delivery",
                   "competing" = "Other modes"
    )
  )
# optional: pick off highlights
highlight_pts <- cif_plot_df %>% filter(time %in% c(12, 24, 48,72))


cif_plot_df %>% 
  filter(method=="Cox‐model weights",clone%in%c("window_0_8","window_16_24","window_32_Inf")) %>% 
ggplot(., aes(x = time, y = cif, color = clone, linetype = event)) +
  geom_step(size = 1) +
  geom_point(data = highlight_pts %>% filter(method=="Cox‐model weights",clone%in%c("window_0_8","window_16_24","window_32_Inf")), size = 2) +
  geom_text(
    data    = highlight_pts%>% filter(method=="Cox‐model weights",clone%in%c("window_0_8","window_16_24","window_32_Inf")),
    aes(label = sprintf("%.2f", cif)),
    vjust   = -0.5, size = 3, show.legend = FALSE
  ) +
  #facet_wrap(~ method, ncol = 1) +
  scale_colour_manual(
    name   = "Treatemt arm",
    # match old levels to new labels
    values = c("window_0_8"    = "#E41A1C",
               "window_16_24"  = "#377EB8",
               "window_32_Inf" = "#4DAF4A"),
    labels = c("0–8", "16–24", "32+")
  ) +
  scale_linetype_manual(
    name   = "Delivery mode",
    # give each old level a line-type,
    # and relabel “Other modes” to “Instrumental delivery + CS”
    values = c("Other modes"       = "dashed",
               "Vaginal delivery"  = "solid"),
    labels = c("Operative delivery", "Vaginal delivery")
  ) +
  labs(
  #  title    = "Cumulative‐Incidence of Vaginal vs. Other Delivery by Method & Clone",
    x        = "Time (hours)",
    y        = "Cumulative incidence",
    color    = "Treatemt arm",
    linetype = "Delivery mode"
  ) +
  theme_minimal(base_size = 8) +
  theme(
    legend.position = c(0.85, 0.6), 
    strip.text      = element_text(face = "bold")
  )


cif_plot_df %>% 
  filter(method=="Logistic‐model weights",clone%in%c("window_0_8","window_16_24","window_32_Inf")) %>% 
  ggplot(., aes(x = time, y = cif, color = clone, linetype = event)) +
  geom_step(size = 1) +
  geom_point(data = highlight_pts %>% filter(method=="Logistic‐model weights",clone%in%c("window_0_8","window_16_24","window_32_Inf")), size = 2) +
  geom_text(
    data    = highlight_pts%>% filter(method=="Logistic‐model weights",clone%in%c("window_0_8","window_16_24","window_32_Inf")),
    aes(label = sprintf("%.2f", cif)),
    vjust   = -0.5, size = 3, show.legend = FALSE
  ) +
  #facet_wrap(~ method, ncol = 1) +
  scale_colour_manual(
    name   = "Treatemt arm",
    # match old levels to new labels
    values = c("window_0_8"    = "#E41A1C",
               "window_16_24"  = "#377EB8",
               "window_32_Inf" = "#4DAF4A"),
    labels = c("0–8", "16–24", "32+")
  ) +
  scale_linetype_manual(
    name   = "Delivery mode",
    # give each old level a line-type,
    # and relabel “Other modes” to “Instrumental delivery + CS”
    values = c("Other modes"       = "dashed",
               "Vaginal delivery"  = "solid"),
    labels = c("Operative delivery", "Vaginal delivery")
  ) +
  labs(
    #  title    = "Cumulative‐Incidence of Vaginal vs. Other Delivery by Method & Clone",
    x        = "Time (hours)",
    y        = "Cumulative incidence",
    color    = "Treatemt arm",
    linetype = "Delivery mode"
  ) +
  theme_minimal(base_size = 8) +
  theme(
    legend.position = c(0.85, 0.6), 
    strip.text      = element_text(face = "bold")
  )

cif_plot_df %>% 
  filter(method=="No weights",clone%in%c("window_0_8","window_16_24","window_32_Inf")) %>% 
  ggplot(., aes(x = time, y = cif, color = clone, linetype = event)) +
  geom_step(size = 1) +
  geom_point(data = highlight_pts %>% filter(method=="No weights",clone%in%c("window_0_8","window_16_24","window_32_Inf")), size = 2) +
  geom_text(
    data    = highlight_pts%>% filter(method=="No weights",clone%in%c("window_0_8","window_16_24","window_32_Inf")),
    aes(label = sprintf("%.2f", cif)),
    vjust   = -0.5, size = 3, show.legend = FALSE
  ) +
  #facet_wrap(~ method, ncol = 1) +
  scale_colour_manual(
    name   = "Treatemt arm",
    # match old levels to new labels
    values = c("window_0_8"    = "#E41A1C",
               "window_16_24"  = "#377EB8",
               "window_32_Inf" = "#4DAF4A"),
    labels = c("0–8", "16–24", "32+")
  ) +
  scale_linetype_manual(
    name   = "Delivery mode",
    # give each old level a line-type,
    # and relabel “Other modes” to “Instrumental delivery + CS”
    values = c("Other modes"       = "dashed",
               "Vaginal delivery"  = "solid"),
    labels = c("Operative delivery", "Vaginal delivery")
  ) +
  labs(
    #  title    = "Cumulative‐Incidence of Vaginal vs. Other Delivery by Method & Clone",
    x        = "Time (hours)",
    y        = "Cumulative incidence",
    color    = "Treatemt arm",
    linetype = "Delivery mode"
  ) +
  theme_minimal(base_size = 8) +
  theme(
    legend.position = c(0.85, 0.6), 
    strip.text      = element_text(face = "bold")
  )




#gc()
ggsave(
  filename = "G:\\My Drive\\Phd4\\overleaf\\plots\\CIF\\SW\\cif_no_weights.pdf",
  plot     = p,
  device   = cairo_pdf,        # or "pdf"
  width    = 7,              # inches
  height   = 4,                # inches
  units    = "in",
  dpi      = 300               # only matters for raster layers
)

p<-cif_plot_df %>% 
#  filter(method=="NO weights",clone%in%c("window_0_8","window_16_24","window_32_Inf")) %>% 
  ggplot(., aes(x = time, y = cif, color = clone, linetype = event)) +
  geom_step(size = 0.5) +
  geom_point(data = highlight_pts , size = 2) +
  geom_text(
    data    = highlight_pts,
    aes(label = sprintf("%.2f", cif)),
    vjust   = -0.5, size = 2, show.legend = FALSE
  ) +
  facet_wrap(~ method, ncol = 1) +
  scale_colour_manual(
    name   = "Treatemt arm",
    # match old levels to new labels
    values = c("window_0_8"    = "#E41A1C",
               "window_16_24"  = "#377EB8",
               "window_32_Inf" = "#4DAF4A"),
    labels = c("0–8", "16–24", "32+")
  ) +
  scale_linetype_manual(
    name   = "Delivery mode",
    # give each old level a line-type,
    # and relabel “Other modes” to “Instrumental delivery + CS”
    values = c("Other modes"       = "dashed",
               "Vaginal delivery"  = "solid"),
    labels = c("Operative delivery", "Vaginal delivery")
  ) +
  labs(
    #  title    = "Cumulative‐Incidence of Vaginal vs. Other Delivery by Method & Clone",
    x        = "Time (hours)",
    y        = "Cumulative incidence",
    color    = "Treatemt arm",
    linetype = "Delivery mode"
  ) +
  theme_minimal(base_size = 8) +
  theme(
    legend.position = "above", 
    strip.text      = element_text(face = "bold")
  )
#gc()
ggsave(
  filename = "G:\\My Drive\\Phd4\\overleaf\\plots\\CIF\\SW\\cif_together.pdf",
  plot     = p,
  device   = cairo_pdf,        # or "pdf"
  width    = 7,              # inches
  height   = 4,                # inches
  units    = "in",
  dpi      = 300               # only matters for raster layers
)
