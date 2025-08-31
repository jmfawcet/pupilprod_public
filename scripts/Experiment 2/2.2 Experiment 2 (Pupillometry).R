
#
# Windowed analysis ----------------------------------------------------------------
# 

# Preprocessing
e2_study_pd %>%
  filter(bin >= 500 & bin <= 3500) %>%
  group_by(sid, condition) %>%
  summarize(pupil_z = mean(pupil_z)) -> e2_study_pd_win

e2_study_pd %>%
  filter(bin >= 500 & bin <= 3500) %>%
  group_by(sid, word) %>%
  summarize(pupil_z = mean(pupil_z)) -> e2_study_pd_per_trial

e2_study_pd %>%
  mutate(mem = ifelse(said_yes==0, 'no', 'yes')) %>%
  filter(bin >= 500 & bin <= 3500) %>%
  group_by(sid, condition, mem) %>%
  summarize(pupil_z = mean(pupil_z)) -> e2_study_pd_win_acc
  
e2_study_pd_win_acc_filtered <- e2_study_pd_win_acc %>%
  group_by(sid) %>%
  filter(n_distinct(condition, mem) == nrow(expand.grid(unique(condition), unique(mem)))) %>%
  ungroup()

# Analysis of the windowed pupil
e2_anova_win = ezANOVA(e2_study_pd_win, .(pupil_z), .(sid), .(condition))
e2_anova_win$Descriptives = ezStats(e2_study_pd_win, .(pupil_z), .(sid), .(condition))

# Analysis of the windowed pupil w/ accuracy
e2_anova_win_acc = ezANOVA(e2_study_pd_win_acc, .(pupil_z), .(sid), .(condition, mem))
e2_anova_win_acc$Descriptives = ezStats(e2_study_pd_win_acc, .(pupil_z), .(sid), .(condition, mem))

#
#
# Functional Data Analysis ------------------------------------------------
#
#

e2_a_m_s = fit_massu_t(e2_study_pd_sum, 'Aloud - Silent', c('aloud', 'silent'))
e2_a_m_c = fit_massu_t(e2_study_pd_sum, 'Aloud - Control', c('aloud', 'control'))
e2_c_m_s = fit_massu_t(e2_study_pd_sum, 'Control - Silent', c('control', 'silent'))

e2_massu_cas = ar1_cluster_correction_t(e2_a_m_s,
                                        alpha          = 0.05,
                                        n_sim          = 3000,     # maybe use 3000 or 5000
                                        two_sided      = TRUE,
                                        df = length(unique(e2_study_pd_sum$sid))-1
)

e2_a_m_s = e2_a_m_s %>%
  mutate(sig_clust = map_lgl(bin, ~ any(e2_massu_cas$observed_clusters$start_tp <= .x & .x <= e2_massu_cas$observed_clusters$end_tp)))


e2_massu_cac = ar1_cluster_correction_t(e2_a_m_c,
                                       alpha          = 0.05,
                                       n_sim          = 3000,     # maybe use 3000 or 5000
                                       two_sided      = TRUE,
                                       df = length(unique(e2_study_pd_sum$sid))-1
)

e2_a_m_c = e2_a_m_c %>%
  mutate(sig_clust = map_lgl(bin, ~ any(e2_massu_cac$observed_clusters$start_tp <= .x & .x <= e2_massu_cac$observed_clusters$end_tp)))


e2_massu_ccs = ar1_cluster_correction_t(e2_c_m_s,
                                       alpha          = 0.05,
                                       n_sim          = 3000,     # maybe use 3000 or 5000
                                       two_sided      = TRUE,
                                       df = length(unique(e2_study_pd_sum$sid))-1
)

e2_c_m_s = e2_c_m_s %>%
  mutate(sig_clust = map_lgl(bin, ~ any(e2_massu_ccs$observed_clusters$start_tp <= .x & .x <= e2_massu_ccs$observed_clusters$end_tp)))


e2_massu = e2_a_m_s %>%
  bind_rows(e2_a_m_c) %>%
  bind_rows(e2_c_m_s)

#
#
# GAMM --------------------------------------------------------------------
#
#

# Analysis ----------------------------------------------------------------

e2_study_pgam_ww = fit_gam_study_w_word(e2_study_pd, 'data/E2/models/e2_study_pgam_ww.rds')
e2_study_mem_pgam = fit_gam_study_mem(e2_study_pd, 'data/E2/models/e2_study_mem_pgam.rds')

#
#
# Predicting Study Phase Accuracy -----------------------------------------
#
#

# Basic Accuracy Model
temp = e2_test_dat %>% 
  filter(condition != 'new') %>%
  left_join(e2_study_pd_per_trial) %>%
  filter(!is.na(pupil_z)) %>%
  mutate(pupil_z = scale(pupil_z)[,1])

e2_bayes_m1_pupilpred = brm(said_yes ~ condition + condition:pupil_z-1 + (condition-1 | sid) + (condition-1 | word),
                  family=brms::bernoulli(link='probit'),
                  backend = 'cmdstan',
                  cores = ncores_brm,
                  chains = nchains_brm,
                  prior = e2_acc_pupilpred_priors,
                  iter = niter_brm,
                  sample_prior = 'yes',
                  threads = threading(4),
                  control = list(adapt_delta = .99),
                  file = 'data/E2/models/e2_bayes_m1_pupilpred',
                  file_refit = 'on_change',
                  data=temp)



#
#
# Rolling Window ----------------------------------------------------------
#
#

e2_rolling_window_cluster = fit_model_across_timepoints_cluster(e2_study_pd, seq(50, 4440, 10), 'data/E2/models/e2_rolling_clust.rds')


e2_rolling_window_cluster_cs <- ar1_cluster_correction(
  results_df     = e2_rolling_window_cluster,
  effect_z_col   = "silent_z",
  effect_p_col   = "silent_p",
  alpha          = 0.05,
  n_sim          = 3000,     # maybe use 3000 or 5000
  two_sided      = TRUE,
  z_thresh       = 1.96
)

e2_rolling_window_cluster_ca <- ar1_cluster_correction(
  results_df     = e2_rolling_window_cluster,
  effect_z_col   = "aloud_z",
  effect_p_col   = "aloud_p",
  alpha          = 0.05,
  n_sim          = 3000,     # maybe use 3000 or 5000
  two_sided      = TRUE,
  z_thresh       = 1.96
)
e2_rolling_window_cluster_ca$observed_clusters = e2_rolling_window_cluster_ca$observed_clusters %>% filter(significant)

e2_rolling_window_cluster_cc <- ar1_cluster_correction(
  results_df     = e2_rolling_window_cluster,
  effect_z_col   = "control_z",
  effect_p_col   = "control_p",
  alpha          = 0.05,
  n_sim          = 3000,     # maybe use 3000 or 5000
  two_sided      = TRUE,
  z_thresh       = 1.96
)

e2_rolling_window_cluster = e2_rolling_window_cluster %>%
  mutate(s_clust = FALSE,
         a_clust = map_lgl(time_point, ~ any(e2_rolling_window_cluster_ca$observed_clusters$start_tp <= .x & .x <= e2_rolling_window_cluster_ca$observed_clusters$end_tp)),
         c_clust = FALSE)
