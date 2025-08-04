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

60*250

#num_cores <- parallel::detectCores() - 5


source("C:\\Users\\Rachel\\Desktop\\SW\\final\\boot_SW_CIF_TV.R")
source("C:\\Users\\Rachel\\Desktop\\SW\\tv\\f_compute_survival_for_clone_SW_TV.R")
df.only.base<- read.csv("C:\\Users\\Rachel\\Downloads\\data_base_line.csv") 
measur_df<-fread("C:\\Users\\Rachel\\Downloads\\opening_height2.csv") 
measur_df<-measur_df %>% select(-height_clean_cf) %>% ###take only the opening data
  mutate(opening_clean_cf=ifelse(time.prom.measurement==0&is.na(opening_clean_cf),0,opening_clean_cf)
  ) %>% 
  filter(!is.na(opening_clean_cf))


df.only.base<-df.only.base %>% right_join(.,measur_df %>% select(id) %>% distinct()) %>% 
  mutate(GBS.pos=ifelse(is.na(GBS)|GBS==0,0,1))


round(colMeans(is.na(df.only.base)) * 100, 1)

#for now we will do complete casews anayslsis-without GBS+proteinuria
list.of.confounder<-age ~parity+GWG+smoking+GBS.pos



vars <- all.vars(list.of.confounder)
df.only.base<-df.only.base %>%
  # select only the vars in the formula, then drop any rows with NA
  filter(
    if_all(all_of(vars), ~ !is.na(.))
  )

brks <- c(0,8,16,24,32,Inf)
brks.e <- c(0,8,16,24,32)


# event_times <- sort(unique(c(0,df.only.base$time.prom.to.delivery, df.only.base$time.prom.oxitocin.hour, brks.e)))
# event_times_df <- data.frame(
#   t_event = event_times,
#   time_id = 1:length(event_times)
# )
make_clone <- function(data, LB, UB = Inf) {
  # build the clone name, e.g. "window_0_12" or "window_12_Inf"
  clone_name <- paste0("window_", LB, "_", ifelse(is.infinite(UB), "Inf", UB))
  
  data %>%
    mutate(
      clone       = clone_name,
      # 1) start everyone at their true delivery time
      fup_outcome = time.prom.to.delivery,
      
      # 2) too‐early induction: censor at oxytocin time if < LB
      fup_outcome = if_else(
        !is.na(time.prom.oxitocin.hour) & time.prom.oxitocin.hour < LB,
        time.prom.oxitocin.hour,
        fup_outcome
      ),
      
      # 3) too‐late induction: censor at oxytocin time if > UB
      fup_outcome = if_else(
        !is.na(time.prom.oxitocin.hour) & time.prom.oxitocin.hour > UB,
        UB,
        fup_outcome
      ),
      
      # 4) no‐oxytocin & delivery after UB: censor at UB
      fup_outcome = if_else(
        is.na(time.prom.oxitocin.hour) & time.prom.to.delivery > UB,
        UB,
        fup_outcome
      ),
      
      # 5) flag event vs. censor
      outcome = as.integer(fup_outcome == time.prom.to.delivery),
      censor  = as.integer(fup_outcome <  time.prom.to.delivery)
    )
}
#threshold<-12
# df.only.base<- df.only.base%>% 
#   mutate(
#     type=case_when(
#       is.na(time.prom.oxitocin.hour)&time.prom.to.delivery <= threshold ~ "A",
#       !is.na(time.prom.oxitocin.hour)&time.prom.oxitocin.hour <= threshold&time.prom.to.delivery <= threshold ~ "B",
#       !is.na(time.prom.oxitocin.hour)&time.prom.oxitocin.hour <= threshold&time.prom.to.delivery > threshold ~ "C",
#       !is.na(time.prom.oxitocin.hour)&time.prom.oxitocin.hour >threshold&time.prom.to.delivery > threshold ~"D",
#       is.na(time.prom.oxitocin.hour)&time.prom.to.delivery > threshold ~ "E",
#       TRUE ~ NA_character_
#     )
#   )
# build all the clones in one data‐frame
df_clones <- map2_dfr(
  brks[-length(brks)],    # LB
  brks[-1],               # UB
  ~ make_clone(df.only.base, .x, .y)
)

# E<-df_clones %>% filter(type=="E")
#I want to explore teh censoring time. types E+D are censored in the first arm.
#B+C are censoerd from the later arm and A is nit censored from both arms
# table(df.only.base$type)
# hist(df.only.base$time.prom.oxitocin.hour)
# median(df.only.base$time.prom.oxitocin.hour,na.rm = TRUE)

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


# glm_formula_first <- (censor == 0) ~ age +bmi+parity
# glm_formula_rest <- (censor == 0) ~ ns(t, df =3) +age +bmi+parity
# cox_formula <- Surv(t_start, t_stop, censor)~ age +bmi+parity
#  


glm_formula_first <- (censor == 0) ~ age +parity+GWG+smoking+GBS.pos

glm_formula_rest <- (censor == 0) ~ ns(t, df =3) +age +parity+GWG+smoking+GBS.pos

cox_formula <- Surv(t_start, t_stop, censor)~ age +parity+GWG+smoking+GBS.pos



clone_info <- tibble(
  LB    = brks[-length(brks)],
  
  
  UB    = brks[-1],
  clone.n = ifelse(is.finite(UB),
                 paste0("window_", LB, "_", UB),
                 paste0("window_", LB, "_Inf"))
)

handlers(global = TRUE)  # Use default progress bar


clone_info_list <- split(clone_info, seq_len(nrow(clone_info)))


start=Sys.time()
cl <- makeCluster(10)
registerDoParallel(cl)
# Not in parallel — just loop over windows
results_list <- list()
#i=1
for (i in seq_along(clone_info_list)) {
  ci <- clone_info_list[[i]]
  clone.n <- ci$clone.n[[1]]
  df_sub  <- df_clones[df_clones$clone == clone.n, , drop = FALSE]
  
  message("Running bootstrap for clone: ", clone.n)
  
  results_list[[i]] <- bootstrap_cif_TV_with_CI(
    df_clone      =      df_sub,
    ci                 = ci,
    only.measur.full.time=NULL,
    glm_formula_first  = glm_formula_first,
    glm_formula_rest   = glm_formula_rest,
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
# now surv_by_clone has columns: time, logit_w, cox_w, unweighted, and clone

# h1<-surv_by_clone%>% filter(clone=="window_0_12") %>% select(time,unw_risk)
# h2<-surv_by_clone%>% filter(clone!="window_0_12") %>% select(unw_risk)
# hh<-cbind(h1,h2)
# Now I can make KM plot for all kind of delivery
# assume your data.frame is called `surv_tbl`
# with columns: time, logit_w, cox_w, unweighted, clone
# 
# 
# cif_plot_df <- cif_by_clone %>%
#   pivot_longer(
#     cols      = c(unw_vaginal, unw_competing,
#                   cox_vaginal, cox_competing,
#                   log_vaginal, log_competing),
#     names_to  = c("method", "event"),
#     names_sep = "_",
#     values_to = "cif"
#   ) %>%
#   mutate(
#     method = recode(method,
#                     unw  = "No weights",
#                     cox  = "Cox‐model weights",
#                     log  = "Logistic‐model weights"
#     ),
#     event = recode(event,
#                    vaginal   = "Vaginal delivery",
#                    competing = "Other modes"
#     )
#   )


cif_plot_df <- cif_by_clone %>%
  mutate(
    event = recode(event,
                   survival="Overall survival",
                   vaginal   = "Vaginal delivery",
                   competing = "Other modes")
  )

write.csv(cif_plot_df, "C:\\Users\\Rachel\\Desktop\\SW\\final\\results\\main\\cif_no_TV.csv", row.names = FALSE)
# optional: pick off highlights
highlight_pts <- cif_plot_df %>% filter(time %in% c(12, 24, 48,72))


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
ggsave("C:\\Users\\Rachel\\Desktop\\SW\\results_SW_no_TV\\CIF_cox_no_TV.pdf", width = 6, height = 4, units = "in")



#Maybe for the main text, only part of teh windows
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
ggsave("C:\\Users\\Rachel\\Desktop\\SW\\results_SW_no_TV\\CIF_cox_SW_no_TV.pdf", width = 6, height = 4, units = "in")


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
ggsave("C:\\Users\\Rachel\\Desktop\\SW\\results_SW_no_TV\\CIF_no_weights_SW_no_TV.pdf", width = 6, height = 4, units = "in")

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
ggsave("C:\\Users\\Rachel\\Desktop\\SW\\results_SW_no_TV\\CIF_logistic_weights_no_TV.pdf", width = 6, height = 4, units = "in")


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
ggsave("C:\\Users\\Rachel\\Desktop\\SW\\results_SW_no_TV\\CIF_logistic_weights_SW_no_TV.pdf", width = 6, height = 4, units = "in")


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

ggsave("C:\\Users\\Rachel\\Downloads\\results_SW_no_TV\\CIF_everthing_SW_no_TV.pdf", width = 6, height = 12, units = "in")

#gc()


