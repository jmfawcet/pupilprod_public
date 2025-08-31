
#
#
# Frequentist Models ------------------------------------------------------
#
#

# Accuracy Summary Statistics ------------------------------------------------------

# Calculate summaries
e2_test_dat %>%
  group_by(sid, condition) %>%
  summarize(yes=mean(said_yes), n=n()) -> e2_test_sum

# Create a wide data frame
e2_test_sum_wide = e2_test_sum %>%
  pivot_wider(id_cols=sid, names_from=condition, values_from=yes) %>%
  mutate(across(
    matches("^(aloud|silent|control)"),  # Apply to all hit rate columns
    ~ qnorm(d_prime_correct(., 40)) - qnorm(d_prime_correct(new, 80)),  # Compute d'
    .names = "d_{.col}"  # Name new columns with "d_"
  ))

# Create a long d' data frame
e2_test_d_sum = e2_test_sum_wide %>%
  select(sid, starts_with('d_')) %>%
  pivot_longer(cols=c(d_aloud:d_silent), names_to='condition', values_to='d') %>%
  mutate(condition = str_remove(condition, "^d_"))
           
# Preprocess data for DPSDT
e2_dpsdt_dat = preprocess_dpsdt(e2_test_dat)

# Combining into a single test phrase
e2_test_sum %>%
  left_join(e2_dpsdt_dat) %>%
  left_join(e2_test_d_sum) -> e2_test_sum_total
 
# ANOVAs --------------------------------------------------------

# Analysis of the "hits"
e2_anova_acc = ezANOVA(e2_test_sum_total, .(yes), .(sid), .(condition))
e2_anova_acc$Descriptives = ezStats(e2_test_sum_total, .(yes), .(sid), .(condition))  %>%
  mutate(SE = SD / sqrt(N))

# Analysis of the d'
e2_anova_d = ezANOVA(e2_test_sum_total %>% filter(condition != 'new'), .(d), .(sid), .(condition))
e2_anova_d$Descriptives = ezStats(e2_test_sum_total %>% filter(condition != 'new'), .(d), .(sid), .(condition))

# Analysis of the Rho
e2_anova_r = ezANOVA(e2_test_sum_total %>% filter(condition != 'new'), .(rho), .(sid), .(condition))
e2_anova_r$Descriptives = ezStats(e2_test_sum_total %>% filter(condition != 'new'), .(rho), .(sid), .(condition))

# Analysis of the Familiarity
e2_anova_f = ezANOVA(e2_test_sum_total %>% filter(condition != 'new'), .(df), .(sid), .(condition))
e2_anova_f$Descriptives = ezStats(e2_test_sum_total %>% filter(condition != 'new'), .(df), .(sid), .(condition))

#
#
# Bayesian Models ------------------------------------------------------
#
#

# Basic Accuracy Model
temp = e2_test_dat %>% mutate(condition=relevel(factor(condition), 'new'))
e2_bayes_m1 = brm(said_yes ~ condition + (condition | sid) + (condition | word),
                  family=brms::bernoulli(link='probit'),
                  backend = 'cmdstan',
                  cores = ncores_brm,
                  chains = nchains_brm,
                  prior = e1_acc_priors,
                  iter = niter_brm,
                  sample_prior = 'yes',
                  threads = threading(4),
                  control = list(adapt_delta = .95),
                  file = 'data/E2/models/e2_bayes_m1',
                  file_refit = 'on_change',
                  data=temp)

e2_bayes_m1_means = emmeans(e2_bayes_m1, ~condition)
e2_bayes_m1_comps = pairs(e2_bayes_m1_means)

# Analysis of the Rho
temp = e2_test_sum_total %>% filter(condition != 'new')
e2_bayes_rho1 = brm(rho ~ condition-1 + (1 | sid),
                  backend = 'cmdstan',
                  cores = ncores_brm,
                  chains = nchains_brm,
                  prior = e2_rho_priors,
                  iter = niter_brm,
                  sample_prior = 'yes',
                  threads = threading(4),
                  control = list(adapt_delta = .999, max_treedepth=20),
                  file = 'data/E2/models/e2_bayes_rho1',
                  file_refit = 'on_change',
                  data=temp)

e2_bayes_rho1_means = emmeans(e2_bayes_rho1, ~condition)
e2_bayes_rho1_comps = pairs(e2_bayes_rho1_means)

# Analysis of the Familiarity
temp = e2_test_sum_total %>% filter(condition != 'new')
e2_bayes_f1 = brm(df ~ condition-1 + (1 | sid),
                    backend = 'cmdstan',
                    cores = ncores_brm,
                    chains = nchains_brm,
                    prior = e2_f_priors,
                    iter = niter_brm,
                    sample_prior = 'yes',
                    threads = threading(4),
                    control = list(adapt_delta = .999, max_treedepth=20),
                    file = 'data/E2/models/e2_bayes_f1',
                    file_refit = 'on_change',
                    data=temp)

e2_bayes_f1_means = emmeans(e2_bayes_f1, ~condition)
e2_bayes_f1_comps = pairs(e2_bayes_f1_means)

# Output Summaries --------------------------------------------------------

write.xlsx(e2_test_sum_wide, 'data/E2/summaries/e2_test_sum_total.xlsx')
