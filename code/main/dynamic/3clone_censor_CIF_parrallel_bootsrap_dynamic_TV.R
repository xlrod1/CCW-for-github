library(dplyr)
library(purrr)
library(data.table)
library(cmprsk)
library(tidyr)
library(ggplot2)
library(survival)
library(foreach)
library(doParallel)
library(progressr)
library(slider)
library(zoo)


#num_cores <- parallel::detectCores() - 5


source("C:\\Users\\Rachel\\Desktop\\SW\\final\\code\\main\\dynamic\\boot_SW_CIF_TV.R")
source("C:\\Users\\Rachel\\Desktop\\SW\\final\\code\\main\\dynamic\\f_compute_survival_for_clone_SW_TV_no_glm.R")
df.only.base<- read.csv("C:\\Users\\Rachel\\Downloads\\data_base_line.csv") 
measur_df<-fread("C:\\Users\\Rachel\\Downloads\\opening_height2.csv") 
measur_df<-measur_df %>% select(-height_clean_cf) %>% ###take only the opening data
  mutate(opening_clean_cf=ifelse(time.prom.measurement==0&is.na(opening_clean_cf),0,opening_clean_cf)
  ) %>% 
  filter(!is.na(opening_clean_cf))


#jj<-measur_df %>% filter(id==2)


#measur_df1<-measur_df %>% filter(bishop.higher.5==1)
#round(colMeans(is.na(df.only.base)) * 100, 1)
df.only.base<-df.only.base %>% right_join(.,measur_df %>% select(id) %>% distinct()) %>% 
  mutate(GBS.pos=ifelse(is.na(GBS)|GBS==0,0,1))

#Opening cut-off X=3 canchange here
c_o<-4
first.time.open<-measur_df%>%
  filter(opening_clean_cf >=c_o) %>%
  group_by(id) %>%
  summarize(first_time_reach_c_o = min(time.prom.measurement)) %>%
  ungroup()



df.only.base<- df.only.base %>% left_join(first.time.open)

#sample.id<-data.frame(id=sample(unique(df.only.base$id),1000))
#df.only.base<-left_join(sample.id,df.only.base)


#for now we will do complete casews anayslsis-without GBS+proteinuria
list.of.confounder<-age ~parity+GWG+smoking+GBS.pos

# list.of.confounder2<-age ~ bmi+parity+ART_pregnancy+
#   Gestational_diabetes+anemia+thrombocytopenia 

# list.of.confounder3<-age ~ bmi+parity
#table(df.only.base$thrombocytopenia)


vars <- all.vars(list.of.confounder)
df.only.base<-df.only.base %>%
  # select only the vars in the formula, then drop any rows with NA
  filter(
    if_all(all_of(vars), ~ !is.na(.))
  )

brks<- c(0,12,Inf)
brks.e <- c(0,12)


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


  


glm_formula <- (censor == 0) ~ ns(t, df =3) +opening_clean_cf +age+parity+GWG+smoking+GBS.pos
cox_formula <- Surv(t_start, t_stop, censor)~opening_clean_cf+age+parity+GWG+smoking+GBS.pos



handlers(global = TRUE)  # Use default progress bar


clone_info_list <- unique(df_clones$clone)


#jj<-df_clones %>% filter(id==2)
# 

start=Sys.time()
cl <- makeCluster(10)
registerDoParallel(cl)
# Not in parallel — just loop over windows
results_list <- list()
#i=3
for (i in seq_along(clone_info_list)) {
  ci <- clone_info_list[[i]]
  #clone.n <- ci$clone.n[[1]]
  df_sub  <- df_clones[df_clones$clone == ci, , drop = FALSE]
  
  message("Running bootstrap for clone: ", ci)
  
  df_clone      =      df_sub
  only.measur.full.time=measur_df
  
  results_list[[i]] <- bootstrap_cif_TV_with_CI(
    df_clone      =      df_sub,
    ci                 = ci,
    only.measur.full.time=measur_df,
   # glm_formula_first  = glm_formula_first,
    glm_formula= glm_formula,
    cox_formula        = cox_formula,
    B                  = 250,
    seed               = 2024 + i  # different seed per clone
  )
}



stop=Sys.time()
stopCluster(cl)

print(min(stop-start))
# row-bind them all into one data frame
cif_by_clone <- bind_rows(results_list)


cif_plot_df <- cif_by_clone %>%
  mutate(
    event = recode(event,
                   survival="Overall survival",
                   vaginal   = "Vaginal delivery",
                   competing = "Other modes")
  )


#table(cif_plot_df$method)



write.csv(cif_plot_df, "C:\\Users\\Rachel\\Desktop\\SW\\final\\results\\main\\cif_dynamic_TV.csv", row.names = FALSE)
# optional: pick off highlights
highlight_pts <- cif_plot_df %>% filter(time %in% c(12, 24, 48,72))

old_cif<-read.csv("C:\\Users\\Rachel\\Desktop\\SW\\results_SW_TV\\cif_SW_TV_CI.csv")
# old_cif<-old_cif %>% 
#   rename("old_cif"="cif",
#          "cif_lower_old"="cif_lower",
#          "cif_upper_old"="cif_upper")
# jj<-cif_plot_df %>% left_join(.,old_cif) %>% 
#   filter(old_cif!=cif)

#jj<-cif_plot_df %>% filter(clone=="window_0_8",method=="Cox‐model weights")















#Maybe for the main text, only part of teh windows
cif_plot_df %>% 
  filter(method=="Cox‐model weights") %>% 
  filter(clone%in%c("window_0_8","window_16_24","window_32_Inf")) %>% 
ggplot(., aes(x = time, y = cif, color = clone, linetype = event)) +
  geom_step(size = 1) +
  geom_point(data = highlight_pts %>% filter(method=="Cox‐model weights",clone%in%c("window_0_8","window_16_24","window_32_Inf")), size = 2) +
  geom_ribbon(aes(ymin = cif_lower, ymax = cif_upper), fill = "blue", alpha = 0.2) +
  geom_text(
    data    = highlight_pts%>% filter(method=="Cox‐model weights",clone%in%c("window_0_8","window_16_24","window_32_Inf")),
    aes(label = sprintf("%.2f", cif)),
    vjust   = -0.5, size = 3, show.legend = FALSE
  ) +
  #facet_wrap(~ method, ncol = 1) +
  scale_color_brewer(palette = "Set1") +
  labs(
  #  title    = "Cumulative‐Incidence of Vaginal vs. Other Delivery by Method & Clone",
    x        = "Time (hours)",
    y        = "Cumulative incidence",
    color    = "Clone",
    linetype = "Event type"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text      = element_text(face = "bold")
  )+ coord_cartesian(xlim = c(0, 73))
ggsave("C:\\Users\\Rachel\\Desktop\\SW\\results_SW_TV\\CIF_cox_not_SW_TV.pdf", width = 6, height = 4, units = "in")


cif_plot_df %>% 
  filter(method=="Cox.SW‐model weights") %>% 
  filter(clone%in%c("window_0_8","window_16_24","window_32_Inf")) %>% 
  ggplot(., aes(x = time, y = cif, color = clone, linetype = event)) +
  geom_step(size = 1) +
  geom_point(data = highlight_pts %>% filter(method=="Cox.SW‐model weights",clone%in%c("window_0_8","window_16_24","window_32_Inf")), size = 2) +
  geom_ribbon(aes(ymin = cif_lower, ymax = cif_upper), fill = "blue", alpha = 0.2) +
  geom_text(
    data    = highlight_pts%>% filter(method=="Cox.SW‐model weights",clone%in%c("window_0_8","window_16_24","window_32_Inf")),
    aes(label = sprintf("%.2f", cif)),
    vjust   = -0.5, size = 3, show.legend = FALSE
  ) +
  #facet_wrap(~ method, ncol = 1) +
  scale_color_brewer(palette = "Set1") +
  labs(
    #  title    = "Cumulative‐Incidence of Vaginal vs. Other Delivery by Method & Clone",
    x        = "Time (hours)",
    y        = "Cumulative incidence",
    color    = "Clone",
    linetype = "Event type"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text      = element_text(face = "bold")
  )+ coord_cartesian(xlim = c(0, 73))


ggsave("C:\\Users\\Rachel\\Desktop\\SW\\results_SW_TV\\CIF_cox_SW_TV.pdf", width = 6, height = 4, units = "in")



#Maybe for the main text, only part of teh windows
cif_plot_df %>% 
  filter(method=="No weights") %>% 
  filter(clone%in%c("window_0_8","window_16_24","window_32_Inf")) %>% 
  ggplot(., aes(x = time, y = cif, color = clone, linetype = event)) +
  geom_step(size = 1) +
  geom_point(data = highlight_pts %>% filter(method=="No weights",clone%in%c("window_0_8","window_16_24","window_32_Inf")), size = 2) +
  geom_ribbon(aes(ymin = cif_lower, ymax = cif_upper), fill = "blue", alpha = 0.2) +
  geom_text(
    data    = highlight_pts%>% filter(method=="No weights",clone%in%c("window_0_8","window_16_24","window_32_Inf")),
    aes(label = sprintf("%.2f", cif)),
    vjust   = -0.5, size = 3, show.legend = FALSE
  ) +
  #facet_wrap(~ method, ncol = 1) +
  scale_color_brewer(palette = "Set1") +
  labs(
  #  title    = "Cumulative‐Incidence of Vaginal vs. Other Delivery by Method & Clone",
    x        = "Time (hours)",
    y        = "Cumulative incidence",
    color    = "Clone",
    linetype = "Event type"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text      = element_text(face = "bold")
  )+
  coord_cartesian(xlim = c(0, 72))
ggsave("C:\\Users\\Rachel\\Desktop\\SW\\results_SW_TV\\CIF_no_weights_SW_TV.pdf", width = 6, height = 4, units = "in")

#Maybe for the main text, only part of teh windows
cif_plot_df %>% 
  filter(method=="Logistic‐model weights") %>% 
  filter(clone%in%c("window_0_8","window_16_24","window_32_Inf")) %>% 
  ggplot(., aes(x = time, y = cif, color = clone, linetype = event)) +
  geom_step(size = 1) +
  geom_point(data = highlight_pts %>% filter(method=="Logistic‐model weights",clone%in%c("window_0_8","window_16_24","window_32_Inf")), size = 2) +
  geom_ribbon(aes(ymin = cif_lower, ymax = cif_upper), fill = "blue", alpha = 0.2) +
  geom_text(
    data    = highlight_pts%>% filter(method=="Logistic‐model weights",clone%in%c("window_0_8","window_16_24","window_32_Inf")),
    aes(label = sprintf("%.2f", cif)),
    vjust   = -0.5, size = 3, show.legend = FALSE
  ) +
  #facet_wrap(~ method, ncol = 1) +
  scale_color_brewer(palette = "Set1") +
  labs(
   # title    = "Cumulative‐Incidence of Vaginal vs. Other Delivery by Method & Clone",
    x        = "Time (hours)",
    y        = "Cumulative incidence",
    color    = "Clone",
    linetype = "Event type"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    strip.text      = element_text(face = "bold")
  )+
  coord_cartesian(xlim = c(0, 72))
ggsave("C:\\Users\\Rachel\\Desktop\\SW\\results_SW_TV\\CIF_logistic_weights_no_SW_TV.pdf", width = 6, height = 4, units = "in")


#Maybe for the main text, only part of teh windows
cif_plot_df %>% 
  filter(method=="Logistic.SW‐model weights") %>% 
  filter(clone%in%c("window_0_8","window_16_24","window_32_Inf")) %>% 
  ggplot(., aes(x = time, y = cif, color = clone, linetype = event)) +
  geom_step(size = 1) +
  geom_point(data = highlight_pts %>% filter(method=="Logistic.SW‐model weights",clone%in%c("window_0_8","window_16_24","window_32_Inf")), size = 2) +
  geom_ribbon(aes(ymin = cif_lower, ymax = cif_upper), fill = "blue", alpha = 0.2) +
  geom_text(
    data    = highlight_pts%>% filter(method=="Logistic.SW‐model weights",clone%in%c("window_0_8","window_16_24","window_32_Inf")),
    aes(label = sprintf("%.2f", cif)),
    vjust   = -0.5, size = 3, show.legend = FALSE
  ) +
  #facet_wrap(~ method, ncol = 1) +
  scale_color_brewer(palette = "Set1") +
  labs(
    # title    = "Cumulative‐Incidence of Vaginal vs. Other Delivery by Method & Clone",
    x        = "Time (hours)",
    y        = "Cumulative incidence",
    color    = "Clone",
    linetype = "Event type"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    strip.text      = element_text(face = "bold")
  )+
  coord_cartesian(xlim = c(0, 72))
ggsave("C:\\Users\\Rachel\\Desktop\\SW\\results_SW_TV\\CIF_logistic_weights_SW_TV.pdf", width = 6, height = 4, units = "in")




#Foe the apedix, all the medhod, all the windows
#Maybe for the main text, only part of teh windows
cif_plot_df %>% 
  #filter(method=="Cox‐model weights") %>% 
 # filter(clone%in%c("window_0_8","window_16_24","window_32_Inf")) %>% 
  ggplot(., aes(x = time, y = cif, color = clone, linetype = event)) +
  geom_step(size = 1) +
  geom_point(data = highlight_pts , size = 2) +
  geom_ribbon(aes(ymin = cif_lower, ymax = cif_upper), fill = "blue", alpha = 0.2) +
  geom_text(
    data    = highlight_pts,
    aes(label = sprintf("%.2f", cif)),
    vjust   = -0.5, size = 3, show.legend = FALSE
  ) +
  facet_wrap(~ method, ncol = 1) +
  scale_color_brewer(palette = "Set1") +
  labs(
   # title    = "Cumulative‐Incidence of Vaginal vs. Other Delivery by Method & Clone",
    x        = "Time (hours)",
    y        = "Cumulative incidence",
    color    = "Clone",
    linetype = "Event type"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    strip.text      = element_text(face = "bold")
  )+
  coord_cartesian(xlim = c(0, 72))

ggsave("C:\\Users\\Rachel\\Downloads\\results_SW_no_TV\\CIF_everthing_SW_TV.pdf", width = 6, height = 12, units = "in")

#gc()


