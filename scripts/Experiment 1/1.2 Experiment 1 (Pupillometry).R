
#
# Windowed analysis ----------------------------------------------------------------
# 

# Preprocessing
e1_study_pd %>%
  filter(bin >= 500) %>%
  group_by(sid, condition) %>%
  summarize(pupil_z = mean(pupil_z)) -> e1_study_pd_win

e1_study_pd %>%
  filter(bin >= 500) %>%
  group_by(sid, word) %>%
  summarize(pupil_z = mean(pupil_z)) -> e1_study_pd_per_trial

e1_study_pd %>%
  mutate(mem = ifelse(said_yes==0, 'no', 'yes')) %>%
  filter(bin >= 500) %>%
  group_by(sid, condition, mem) %>%
  summarize(pupil_z = mean(pupil_z)) -> e1_study_pd_win_acc

e1_study_pd_win_acc_filtered <- e1_study_pd_win_acc %>%
  group_by(sid) %>%
  filter(n_distinct(condition, mem) == nrow(expand.grid(unique(condition), unique(mem)))) %>%
  ungroup()
  
# Analysis of the windowed pupil
e1_anova_win = ezANOVA(e1_study_pd_win, .(pupil_z), .(sid), .(condition))
e1_anova_win$Descriptives = ezStats(e1_study_pd_win, .(pupil_z), .(sid), .(condition))

# Analysis of the windowed pupil w/ accuracy
e1_anova_win_acc = ezANOVA(e1_study_pd_win_acc_filtered, .(pupil_z), .(sid), .(condition, mem))
e1_anova_win_acc$Descriptives = ezStats(e1_study_pd_win_acc_filtered, .(pupil_z), .(sid), .(condition, mem))

#
#
# Functional Data Analysis ------------------------------------------------
#
#

e1_massu = fit_massu_t(e1_study_pd_sum, 'Aloud - Silent', c('aloud', 'silent'))

e1_massu_as = ar1_cluster_correction_t(e1_massu,
                                        alpha          = 0.05,
                                        n_sim          = 3000,     # maybe use 3000 or 5000
                                        two_sided      = TRUE,
                                        df = length(unique(e1_study_pd_sum$sid))-1
                                        )

e1_massu = e1_massu %>%
  mutate(as_clust = map_lgl(bin, ~ any(e1_massu_as$observed_clusters$start_tp <= .x & .x <= e1_massu_as$observed_clusters$end_tp)))

#
#
# GAMM --------------------------------------------------------------------
#
#

# Analysis ----------------------------------------------------------------

e1_study_pgam_ww = fit_gam_study_w_word(e1_study_pd, 'data/E1/models/e1_study_pgam_ww.rds')
e1_study_mem_pgam = fit_gam_study_mem(e1_study_pd, 'data/E1/models/e1_study_mem_pgam.rds')

#
#
# Predicting Study Phase Accuracy -----------------------------------------
#
#

# Basic Accuracy Model
temp = e1_test_dat %>% 
  filter(condition != 'new') %>%
  left_join(e1_study_pd_per_trial) %>%
  filter(!is.na(pupil_z)) %>%
  mutate(pupil_z = scale(pupil_z)[,1])

e1_bayes_m1_pupilpred = brm(said_yes ~ condition + condition:pupil_z-1 + (condition-1 | sid) + (condition-1 | word),
                            family=brms::bernoulli(link='probit'),
                            backend = 'cmdstan',
                            cores = ncores_brm,
                            chains = nchains_brm,
                            prior = e1_acc_pupilpred_priors,
                            iter = niter_brm,
                            sample_prior = 'yes',
                            threads = threading(4),
                            control = list(adapt_delta = .99),
                            file = 'data/E1/models/e1_bayes_m1_pupilpred',
                            file_refit = 'on_change',
                            data=temp)

#
#
# Rolling Window ----------------------------------------------------------
#
#

e1_rolling_window_cluster = fit_model_across_timepoints_cluster(e1_study_pd, seq(10, 1990,10), 'data/E1/models/e1_rolling_clust.rds')

e1_rolling_window_cluster_cs <- ar1_cluster_correction(
  results_df     = e1_rolling_window_cluster,
  effect_z_col   = "silent_z",
  effect_p_col   = "silent_p",
  alpha          = 0.05,
  n_sim          = 3000,     # maybe use 3000 or 5000
  two_sided      = TRUE,
  z_thresh       = 1.96
)

e1_rolling_window_cluster_ca <- ar1_cluster_correction(
  results_df     = e1_rolling_window_cluster,
  effect_z_col   = "aloud_z",
  effect_p_col   = "aloud_p",
  alpha          = 0.05,
  n_sim          = 3000,     # maybe use 3000 or 5000
  two_sided      = TRUE,
  z_thresh       = 1.96
)

e1_rolling_window_cluster = e1_rolling_window_cluster %>%
  mutate(s_clust = map_lgl(time_point, ~ any(e1_rolling_window_cluster_cs$observed_clusters$start_tp <= .x & .x <= e1_rolling_window_cluster_cs$observed_clusters$end_tp)),
         a_clust = FALSE)


