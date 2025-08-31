e1_study_pd_win %>%
  pivot_wider(id_cols=sid, names_from=condition, values_from=pupil_z) %>%
  mutate(ppe = aloud-silent) -> e1_pupillary_pe_sum

e1_test_sum_total %>%
  filter(condition != 'new') %>%
  select(-n, -yes) %>%
  pivot_wider(id_cols=sid, names_from=condition, values_from=c(d)) %>%
  mutate(dpe = aloud-silent) %>%
  rename(a = aloud, s = silent) %>%
  right_join(e1_pupillary_pe_sum) -> e1_behav_pe_sum

e2_study_pd_win %>%
  pivot_wider(id_cols=sid, names_from=condition, values_from=pupil_z) %>%
  mutate(ppe = aloud-silent, cppe=control-silent)  -> e2_pupillary_pe_sum

e2_test_sum_total %>%
  filter(condition != 'new') %>%
  select(-n, -yes) %>%
  pivot_wider(id_cols=sid, names_from=condition, values_from=c(rho, df, d)) %>%
  mutate(dpe = d_aloud-d_silent, cdpe=d_control-d_silent,
         rpe = rho_aloud - rho_silent, crpe = rho_control - rho_silent,
         fpe = df_aloud - df_silent, cfpe = df_control - df_silent) %>%
  #select_at(vars(contains('pe'))) %>%
  right_join(e2_pupillary_pe_sum) -> e2_behav_pe_sum

e3_study_pd_win %>%
  pivot_wider(id_cols=sid, names_from=c(condition, window), values_from=pupil_z) %>%
  mutate(e_ppe = aloud_early-silent_early, e_cppe=control_early-silent_early,
         l_ppe = aloud_late-silent_late, l_cppe=control_late-silent_late,
         c_ppe = aloud_complete-silent_complete, c_cppe=control_complete-silent_complete)  -> e3_pupillary_pe_sum

e3_test_sum_total %>%
  filter(condition != 'new') %>%
  select(-n, -yes) %>%
  pivot_wider(id_cols=sid, names_from=condition, values_from=c(rho, df, d)) %>%
  mutate(dpe = d_aloud-d_silent, cdpe=d_control-d_silent,
         rpe = rho_aloud - rho_silent, crpe = rho_control - rho_silent,
         fpe = df_aloud - df_silent, cfpe = df_control - df_silent) %>%
  #select_at(vars(contains('pe'))) %>%
  right_join(e3_pupillary_pe_sum) -> e3_behav_pe_sum


e4_study_pd_win %>%
  pivot_wider(id_cols=sid, names_from=c(condition, window), values_from=pupil_z) %>%
  mutate(e_ppe = aloud_early-silent_early, e_cppe=control_early-silent_early,
         l_ppe = aloud_late-silent_late, l_cppe=control_late-silent_late,
         ppe = aloud_complete-silent_complete, cppe=control_complete-silent_complete) -> e4_pupillary_pe_sum

e4_test_sum_total %>%
  filter(condition != 'new') %>%
  select(-n, -yes) %>%
  pivot_wider(id_cols=sid, names_from=condition, values_from=c(rho, df, d)) %>%
  mutate(dpe = d_aloud-d_silent, cdpe=d_control-d_silent,
         rpe = rho_aloud - rho_silent, crpe = rho_control - rho_silent,
         fpe = df_aloud - df_silent, cfpe = df_control - df_silent) %>%
  #select_at(vars(contains('pe'))) %>%
  right_join(e4_pupillary_pe_sum) -> e4_behav_pe_sum

e1_behav_pe_sum %>% 
  mutate(aloud_complete = NA, aloud_early = NA, aloud_late = NA,
         control_complete = NA, control_early = NA, control_late = NA,
         silent_complete = NA, silent_early = NA, silent_late = NA,
         d_control = NA, exp='E1') %>%
  select(sid, exp, d_aloud=a, d_control, d_silent=s,
         aloud_complete=aloud, aloud_early, aloud_late,
         control_complete, control_early, control_late,
         silent_complete=silent, silent_early, silent_late) -> temp1

e2_behav_pe_sum %>% 
  mutate(aloud_complete = aloud, aloud_early = NA, aloud_late = NA,
         control_complete = control, control_early = NA, control_late = NA,
         silent_complete = silent, silent_early = NA, silent_late = NA, exp='E2') %>%
  select(sid, exp, d_aloud, d_control, d_silent,
         aloud_complete, aloud_early, aloud_late,
         control_complete, control_early, control_late,
         silent_complete, silent_early, silent_late)  -> temp2

e3_behav_pe_sum %>% 
  mutate(exp='E3') %>%
  select(sid, exp, d_aloud, d_control, d_silent,
         aloud_complete, aloud_early, aloud_late,
         control_complete, control_early, control_late,
         silent_complete, silent_early, silent_late)  -> temp3

e4_behav_pe_sum %>% 
  mutate(exp='E4') %>%
  select(sid, exp, d_aloud, d_control, d_silent,
         aloud_complete, aloud_early, aloud_late,
         control_complete, control_early, control_late,
         silent_complete, silent_early, silent_late)  -> temp4

temp_comb = bind_rows(temp1, temp2, temp3, temp4)

# Define the pairs of variables you want to plot
pairs_list <- list(
  c("d_aloud", "aloud_complete"),
  c("d_aloud", "aloud_early"),
  c("d_aloud", "aloud_late"),
  
  c("d_control", "control_complete"),
  c("d_control", "control_early"),
  c("d_control", "control_late"),
  
  c("d_silent", "silent_complete"),
  c("d_silent", "silent_early"),
  c("d_silent", "silent_late")
)

cor_results_indv <- compute_correlations_with_outlier_removal(temp_comb[-1], group_var = 'exp', pairs_list)
cor_results_indv_overall <- compute_correlations_with_outlier_removal(temp_comb[-1], pairs_list)
