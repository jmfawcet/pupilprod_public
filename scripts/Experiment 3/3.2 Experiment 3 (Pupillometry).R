
#
# Windowed analysis ----------------------------------------------------------------
# 

# Preprocessing
# Note that 0 - 1500 ms is the period the instruction was on the screen, with
# the word coming on from 1500 - 5000 ms and a fixation from 5000 - 6000 ms.
e3_study_pd %>%
  filter(bin >= 500 & bin <= 5000) %>%
  mutate(window = ifelse(bin <= 1500, 'early', 'late')) %>%
  group_by(sid, condition, window) %>%
  summarize(pupil_z = mean(pupil_z)) %>%
  bind_rows(
    e3_study_pd %>%
      filter(bin >= 500 & bin <= 5000) %>%
      group_by(sid, condition) %>%
      summarize(pupil_z = mean(pupil_z), .groups = "drop") %>%
      mutate(window = "complete")  # Add the "complete" window label
  ) %>% arrange(sid, condition, window)-> e3_study_pd_win

e3_study_pd %>%
  filter(bin >= 500 & bin <= 5000) %>%
  mutate(window = ifelse(bin <= 1500, 'early', 'late')) %>%
  group_by(sid, word, window) %>%
  summarize(pupil_z = mean(pupil_z)) %>%
  bind_rows(
    e3_study_pd %>%
      filter(bin >= 500 & bin <= 5000) %>%
      group_by(sid, word) %>%
      summarize(pupil_z = mean(pupil_z), .groups = "drop") %>%
      mutate(window = "complete")  # Add the "complete" window label
  ) %>% pivot_wider(names_from='window', values_from='pupil_z') -> e3_study_pd_per_trial

e3_study_pd %>%
  filter(bin >= 500 & bin <= 5000) %>%
  mutate(window = ifelse(bin <= 1500, 'early', 'late')) %>%
  mutate(mem = ifelse(said_yes==0, 'no', 'yes')) %>%
  filter(bin >= 500) %>%
  group_by(sid, condition, mem, window) %>%
  summarize(pupil_z = mean(pupil_z)) %>%
  bind_rows(
    e3_study_pd %>%
      filter(bin >= 500 & bin <= 5000) %>%
      group_by(sid, condition, mem) %>%
      mutate(mem = ifelse(said_yes==0, 'no', 'yes')) %>%
      summarize(pupil_z = mean(pupil_z), .groups = "drop") %>%
      mutate(window = "complete")  # Add the "complete" window label
  ) %>% arrange(sid, condition, mem, window)-> e3_study_pd_win_acc

e3_study_pd_win_acc_filtered <- e3_study_pd_win_acc %>%
  group_by(sid) %>%
  filter(n_distinct(condition, mem, window) == nrow(expand.grid(unique(condition), unique(mem), unique(window)))) %>%
  ungroup()

# Analysis of the windowed pupil
e3_study_pd_win_both = e3_study_pd_win %>% filter(window != 'complete')
e3_anova_win_both = ezANOVA(e3_study_pd_win_both , .(pupil_z), .(sid), .(condition, window))
e3_anova_win_both$Descriptives = ezStats(e3_study_pd_win_both, .(pupil_z), .(sid), .(condition, window))

e3_study_pd_win_comp = e3_study_pd_win %>% filter(window == 'complete')
e3_anova_win = ezANOVA(e3_study_pd_win_comp , .(pupil_z), .(sid), .(condition))
e3_anova_win$Descriptives = ezStats(e3_study_pd_win_comp, .(pupil_z), .(sid), .(condition))

# Analysis of the windowed pupil w/ accuracy
e3_study_pd_win_acc_filtered_both = e3_study_pd_win_acc_filtered %>% filter(window != 'complete', !is.na(mem))
e3_anova_win_acc_both = ezANOVA(e3_study_pd_win_acc_filtered_both, .(pupil_z), .(sid), .(condition, mem, window))
e3_anova_win_acc_both$Descriptives = ezStats(e3_study_pd_win_acc_filtered_both, .(pupil_z), .(sid), .(condition, mem, window))

e3_study_pd_win_acc_filtered_comp = e3_study_pd_win_acc_filtered %>% filter(window == 'complete', !is.na(mem))
e3_anova_win_acc = ezANOVA(e3_study_pd_win_acc_filtered_comp, .(pupil_z), .(sid), .(condition, mem))
e3_anova_win_acc$Descriptives = ezStats(e3_study_pd_win_acc_filtered_comp, .(pupil_z), .(sid), .(condition, mem))

#
#
# Functional Data Analysis ------------------------------------------------
#
#

e3_a_m_s = fit_massu_t(e3_study_pd_sum, 'Aloud - Silent', c('aloud', 'silent'))
e3_a_m_c = fit_massu_t(e3_study_pd_sum, 'Aloud - Control', c('aloud', 'control'))
e3_c_m_s = fit_massu_t(e3_study_pd_sum, 'Control - Silent', c('control', 'silent'))

e3_massu_cas = ar1_cluster_correction_t(e3_a_m_s,
                                        alpha          = 0.05,
                                        n_sim          = 3000,     # maybe use 3000 or 5000
                                        two_sided      = TRUE,
                                        df = length(unique(e3_study_pd_sum$sid))-1
)

e3_a_m_s = e3_a_m_s %>%
  mutate(sig_clust = map_lgl(bin, ~ any(e3_massu_cas$observed_clusters$start_tp <= .x & .x <= e3_massu_cas$observed_clusters$end_tp)))


e3_massu_cac = ar1_cluster_correction_t(e3_a_m_c,
                                        alpha          = 0.05,
                                        n_sim          = 3000,     # maybe use 3000 or 5000
                                        two_sided      = TRUE,
                                        df = length(unique(e3_study_pd_sum$sid))-1
)

e3_a_m_c = e3_a_m_c %>%
  mutate(sig_clust = map_lgl(bin, ~ any(e3_massu_cac$observed_clusters$start_tp <= .x & .x <= e3_massu_cac$observed_clusters$end_tp)))


e3_massu_ccs = ar1_cluster_correction_t(e3_c_m_s,
                                        alpha          = 0.05,
                                        n_sim          = 3000,     # maybe use 3000 or 5000
                                        two_sided      = TRUE,
                                        df = length(unique(e3_study_pd_sum$sid))-1
)

e3_c_m_s = e3_c_m_s %>%
  mutate(sig_clust = map_lgl(bin, ~ any(e3_massu_ccs$observed_clusters$start_tp <= .x & .x <= e3_massu_ccs$observed_clusters$end_tp)))

e3_massu = e3_a_m_s %>%
  bind_rows(e3_a_m_c) %>%
  bind_rows(e3_c_m_s)

#
#
# GAMM --------------------------------------------------------------------
#
#

# Analysis ----------------------------------------------------------------

e3_study_pgam_ww = fit_gam_study_w_word(e3_study_pd, 'data/E3/models/e3_study_pgam_ww.rds')
e3_study_mem_pgam = fit_gam_study_mem(e3_study_pd, 'data/E3/models/e3_study_mem_pgam.rds')

#
#
# Predicting Study Phase Accuracy -----------------------------------------
#
#

# Basic Accuracy Model
temp = e3_test_dat %>% 
  filter(condition != 'new') %>%
  left_join(e3_study_pd_per_trial) %>%
  filter(!is.na(early) & !is.na(late)) %>%
  mutate(early = scale(early)[,1], late = scale(late)[,1])

e3_bayes_m1_pupilpred = brm(said_yes ~ condition + condition:early + condition:late-1 + (condition-1 | sid) + (condition-1 | word),
                            family=brms::bernoulli(link='probit'),
                            backend = 'cmdstan',
                            cores = ncores_brm,
                            chains = nchains_brm,
                            prior = e2_acc_pupilpred_priors,
                            iter = niter_brm,
                            sample_prior = 'yes',
                            threads = threading(4),
                            control = list(adapt_delta = .95),
                            file = 'data/E3/models/e3_bayes_m1_pupilpred',
                            file_refit = 'on_change',
                            data=temp)

# Basic Accuracy Model using overall pupil size, not windowed
temp = e3_test_dat %>% 
  filter(condition != 'new') %>%
  left_join(e3_study_pd_per_trial) %>%
  filter(!is.na(early) & !is.na(late)) %>%
  mutate(early = scale(early)[,1], late = scale(late)[,1], complete = scale(complete)[,1])

e3_bayes_m1_pupilpred_comp = brm(said_yes ~ condition + condition:complete-1 + (condition-1 | sid) + (condition-1 | word),
                            family=brms::bernoulli(link='probit'),
                            backend = 'cmdstan',
                            cores = ncores_brm,
                            chains = nchains_brm,
                            prior = e2_acc_pupilpred_priors,
                            iter = niter_brm,
                            sample_prior = 'yes',
                            threads = threading(4),
                            control = list(adapt_delta = .95),
                            file = 'data/E3/models/e3_bayes_m1_pupilpred_comp',
                            file_refit = 'on_change',
                            data=temp)

e3_bayes_m1_pupilpred_early = brm(said_yes ~ condition + condition:early-1 + (condition-1 | sid) + (condition-1 | word),
                                 family=brms::bernoulli(link='probit'),
                                 backend = 'cmdstan',
                                 cores = ncores_brm,
                                 chains = nchains_brm,
                                 prior = e2_acc_pupilpred_priors,
                                 iter = niter_brm,
                                 sample_prior = 'yes',
                                 threads = threading(4),
                                 control = list(adapt_delta = .95),
                                 file = 'data/E3/models/e3_bayes_m1_pupilpred_early',
                                 file_refit = 'on_change',
                                 data=temp)

e3_bayes_m1_pupilpred_late = brm(said_yes ~ condition + condition:late-1 + (condition-1 | sid) + (condition-1 | word),
                                 family=brms::bernoulli(link='probit'),
                                 backend = 'cmdstan',
                                 cores = ncores_brm,
                                 chains = nchains_brm,
                                 prior = e2_acc_pupilpred_priors,
                                 iter = niter_brm,
                                 sample_prior = 'yes',
                                 threads = threading(4),
                                 control = list(adapt_delta = .95),
                                 file = 'data/E3/models/e3_bayes_m1_pupilpred_late',
                                 file_refit = 'on_change',
                                 data=temp)

#
#
# Rolling Window ----------------------------------------------------------
#
#

e3_rolling_window_cluster = fit_model_across_timepoints_cluster(e3_study_pd, seq(50, 5940,10), 'data/E3/models/e3_rolling_clust.rds')

e3_rolling_window_cluster_cs <- ar1_cluster_correction(
  results_df     = e3_rolling_window_cluster,
  effect_z_col   = "silent_z",
  effect_p_col   = "silent_p",
  alpha          = 0.05,
  n_sim          = 3000,     # maybe use 3000 or 5000
  two_sided      = TRUE,
  z_thresh       = 1.96
)
e3_rolling_window_cluster_cs$observed_clusters = e3_rolling_window_cluster_cs$observed_clusters %>% filter(significant)

e3_rolling_window_cluster_ca <- ar1_cluster_correction(
  results_df     = e3_rolling_window_cluster,
  effect_z_col   = "aloud_z",
  effect_p_col   = "aloud_p",
  alpha          = 0.05,
  n_sim          = 3000,     # maybe use 3000 or 5000
  two_sided      = TRUE,
  z_thresh       = 1.96
)

e3_rolling_window_cluster_cc <- ar1_cluster_correction(
  results_df     = e3_rolling_window_cluster,
  effect_z_col   = "control_z",
  effect_p_col   = "control_p",
  alpha          = 0.05,
  n_sim          = 3000,     # maybe use 3000 or 5000
  two_sided      = TRUE,
  z_thresh       = 1.96
)

e3_rolling_window_cluster = e3_rolling_window_cluster %>%
  mutate(s_clust = map_lgl(time_point, ~ any(e3_rolling_window_cluster_cs$observed_clusters$start_tp <= .x & .x <= e3_rolling_window_cluster_cs$observed_clusters$end_tp)),
         a_clust = FALSE,
         c_clust = FALSE)
