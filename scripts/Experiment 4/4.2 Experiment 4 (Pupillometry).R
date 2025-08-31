
#
# Windowed analysis ----------------------------------------------------------------
# 

# Preprocessing
# Early = 500 - 3500, Late = 3500 - 6000
e4_study_pd %>%
  filter(bin >= 500 & bin <= 6000) %>%
  mutate(window = ifelse(bin <= 3500, 'early', 'late')) %>%
  group_by(sid, condition, window) %>%
  summarize(pupil_z = mean(pupil_z)) %>%
  bind_rows(
    e4_study_pd %>%
      filter(bin >= 500 & bin <= 6000) %>%
      group_by(sid, condition) %>%
      summarize(pupil_z = mean(pupil_z), .groups = "drop") %>%
      mutate(window = "complete")  # Add the "complete" window label
  ) %>% arrange(sid, condition, window)-> e4_study_pd_win

e4_study_pd %>%
  filter(bin >= 500 & bin <= 6000) %>%
  mutate(window = ifelse(bin <= 3500, 'early', 'late')) %>%
  group_by(sid, word, window) %>%
  summarize(pupil_z = mean(pupil_z)) %>%
  bind_rows(
    e4_study_pd %>%
      filter(bin >= 500 & bin <= 6000) %>%
      group_by(sid, word) %>%
      summarize(pupil_z = mean(pupil_z), .groups = "drop") %>%
      mutate(window = "complete")  # Add the "complete" window label
  ) %>% pivot_wider(names_from='window', values_from='pupil_z') -> e4_study_pd_per_trial

e4_study_pd %>%
  filter(bin >= 500 & bin <= 6000) %>%
  mutate(window = ifelse(bin <= 3500, 'early', 'late')) %>%
  mutate(mem = ifelse(said_yes==0, 'no', 'yes')) %>%
  filter(bin >= 500) %>%
  group_by(sid, condition, mem, window) %>%
  summarize(pupil_z = mean(pupil_z)) %>%
  bind_rows(
    e4_study_pd %>%
      filter(bin >= 500 & bin <= 5000) %>%
      group_by(sid, condition, mem) %>%
      mutate(mem = ifelse(said_yes==0, 'no', 'yes')) %>%
      summarize(pupil_z = mean(pupil_z), .groups = "drop") %>%
      mutate(window = "complete")  # Add the "complete" window label
  ) %>% arrange(sid, condition, mem, window) -> e4_study_pd_win_acc

e4_study_pd_win_acc_filtered <- e4_study_pd_win_acc %>%
  group_by(sid) %>%
  filter(n_distinct(condition, mem, window) == nrow(expand.grid(unique(condition), unique(mem), unique(window)))) %>%
  ungroup()

# Analysis of the windowed pupil
e4_study_pd_win_both = e4_study_pd_win %>% filter(window != 'complete')
e4_anova_win_both = ezANOVA(e4_study_pd_win_both , .(pupil_z), .(sid), .(condition, window))
e4_anova_win_both$Descriptives = ezStats(e4_study_pd_win_both, .(pupil_z), .(sid), .(condition, window))

e4_study_pd_win_comp = e4_study_pd_win %>% filter(window == 'complete')
e4_anova_win = ezANOVA(e4_study_pd_win_comp , .(pupil_z), .(sid), .(condition))
e4_anova_win$Descriptives = ezStats(e4_study_pd_win_comp, .(pupil_z), .(sid), .(condition))

# Analysis of the windowed pupil w/ accuracy
e4_study_pd_win_acc_filtered_both = e4_study_pd_win_acc_filtered %>% filter(window != 'complete', !is.na(mem))
e4_anova_win_acc = ezANOVA(e4_study_pd_win_acc_filtered_both, .(pupil_z), .(sid), .(condition, mem, window))
e4_anova_win_acc$Descriptives = ezStats(e4_study_pd_win_acc_filtered_both, .(pupil_z), .(sid), .(condition, mem, window))

e4_study_pd_win_acc_filtered_comp = e4_study_pd_win_acc_filtered %>% filter(window == 'complete', !is.na(mem))
e4_anova_win_acc = ezANOVA(e4_study_pd_win_acc_filtered_comp, .(pupil_z), .(sid), .(condition, mem))
e4_anova_win_acc$Descriptives = ezStats(e4_study_pd_win_acc_filtered_comp, .(pupil_z), .(sid), .(condition, mem))

#
#
# Functional Data Analysis ------------------------------------------------
#
#

e4_a_m_s = fit_massu_t(e4_study_pd_sum, 'Aloud - Silent', c('aloud', 'silent'))
e4_a_m_c = fit_massu_t(e4_study_pd_sum, 'Aloud - Control', c('aloud', 'control'))
e4_c_m_s = fit_massu_t(e4_study_pd_sum, 'Control - Silent', c('control', 'silent'))

e4_massu_cas = ar1_cluster_correction_t(e4_a_m_s,
                                        alpha          = 0.05,
                                        n_sim          = 3000,     # maybe use 3000 or 5000
                                        two_sided      = TRUE,
                                        df = length(unique(e4_study_pd_sum$sid))-1
)

e4_a_m_s = e4_a_m_s %>%
  mutate(sig_clust = map_lgl(bin, ~ any(e4_massu_cas$observed_clusters$start_tp <= .x & .x <= e4_massu_cas$observed_clusters$end_tp)))


e4_massu_cac = ar1_cluster_correction_t(e4_a_m_c,
                                        alpha          = 0.05,
                                        n_sim          = 3000,     # maybe use 3000 or 5000
                                        two_sided      = TRUE,
                                        df = length(unique(e4_study_pd_sum$sid))-1
)

e4_a_m_c = e4_a_m_c %>%
  mutate(sig_clust = map_lgl(bin, ~ any(e4_massu_cac$observed_clusters$start_tp <= .x & .x <= e4_massu_cac$observed_clusters$end_tp)))


e4_massu_ccs = ar1_cluster_correction_t(e4_c_m_s,
                                        alpha          = 0.05,
                                        n_sim          = 3000,     # maybe use 3000 or 5000
                                        two_sided      = TRUE,
                                        df = length(unique(e4_study_pd_sum$sid))-1
)

e4_c_m_s = e4_c_m_s %>%
  mutate(sig_clust = map_lgl(bin, ~ any(e4_massu_ccs$observed_clusters$start_tp <= .x & .x <= e4_massu_ccs$observed_clusters$end_tp)))

e4_massu = e4_a_m_s %>%
  bind_rows(e4_a_m_c) %>%
  bind_rows(e4_c_m_s)

#
#
# GAMM --------------------------------------------------------------------
#
#

# Analysis ----------------------------------------------------------------

e4_study_pgam_ww = fit_gam_study_w_word(e4_study_pd, 'data/E4/models/e4_study_pgam_ww.rds')
e4_study_mem_pgam = fit_gam_study_mem(e4_study_pd, 'data/E4/models/e4_study_mem_pgam.rds')

#
#
# Predicting Study Phase Accuracy -----------------------------------------
#
#

# Basic Accuracy Model
temp = e4_test_dat %>% 
  filter(condition != 'new') %>%
  left_join(e4_study_pd_per_trial) %>%
  filter(!is.na(early) & !is.na(late)) %>%
  mutate(early = scale(early)[,1], late = scale(late)[,1])

e4_bayes_m1_pupilpred = brm(said_yes ~ condition + condition:early + condition:late-1 + (condition-1 | sid) + (condition-1 | word),
                            family=brms::bernoulli(link='probit'),
                            backend = 'cmdstan',
                            cores = ncores_brm,
                            chains = nchains_brm,
                            prior = e2_acc_pupilpred_priors,
                            iter = niter_brm,
                            sample_prior = 'yes',
                            threads = threading(4),
                            control = list(adapt_delta = .95),
                            file = 'data/E4/models/e4_bayes_m1_pupilpred',
                            file_refit = 'on_change',
                            data=temp)

# Basic Accuracy Model using overall pupil size, not windowed
temp = e4_test_dat %>% 
  filter(condition != 'new') %>%
  left_join(e4_study_pd_per_trial) %>%
  filter(!is.na(early) & !is.na(late)) %>%
  mutate(early = scale(early)[,1], late = scale(late)[,1], complete = scale(complete)[,1])

e4_bayes_m1_pupilpred_comp = brm(said_yes ~ condition + condition:complete-1 + (condition-1 | sid) + (condition-1 | word),
                                 family=brms::bernoulli(link='probit'),
                                 backend = 'cmdstan',
                                 cores = ncores_brm,
                                 chains = nchains_brm,
                                 prior = e2_acc_pupilpred_priors,
                                 iter = niter_brm,
                                 sample_prior = 'yes',
                                 threads = threading(4),
                                 control = list(adapt_delta = .95),
                                 file = 'data/E4/models/e4_bayes_m1_pupilpred_comp',
                                 file_refit = 'on_change',
                                 data=temp)

#
#
# Rolling Window ----------------------------------------------------------
#
#

e4_rolling_window_cluster = fit_model_across_timepoints_cluster(e4_study_pd, seq(50, 5940,10), 'data/E4/models/e4_rolling_clust.rds')

e4_rolling_window_cluster_cs <- ar1_cluster_correction(
  results_df     = e4_rolling_window_cluster,
  effect_z_col   = "silent_z",
  effect_p_col   = "silent_p",
  alpha          = 0.05,
  n_sim          = 3000,     # maybe use 3000 or 5000
  two_sided      = TRUE,
  z_thresh       = 1.96
)

e4_rolling_window_cluster_ca <- ar1_cluster_correction(
  results_df     = e4_rolling_window_cluster,
  effect_z_col   = "aloud_z",
  effect_p_col   = "aloud_p",
  alpha          = 0.05,
  n_sim          = 3000,     # maybe use 3000 or 5000
  two_sided      = TRUE,
  z_thresh       = 1.96
)

e4_rolling_window_cluster_cc <- ar1_cluster_correction(
  results_df     = e4_rolling_window_cluster,
  effect_z_col   = "control_z",
  effect_p_col   = "control_p",
  alpha          = 0.05,
  n_sim          = 3000,     # maybe use 3000 or 5000
  two_sided      = TRUE,
  z_thresh       = 1.96
)

e4_rolling_window_cluster = e4_rolling_window_cluster %>%
  mutate(s_clust = map_lgl(time_point, ~ any(e4_rolling_window_cluster_cs$observed_clusters$start_tp <= .x & .x <= e4_rolling_window_cluster_cs$observed_clusters$end_tp)),
         a_clust = FALSE,
         c_clust = map_lgl(time_point, ~ any(e4_rolling_window_cluster_cc$observed_clusters$start_tp <= .x & .x <= e4_rolling_window_cluster_cc$observed_clusters$end_tp)))
