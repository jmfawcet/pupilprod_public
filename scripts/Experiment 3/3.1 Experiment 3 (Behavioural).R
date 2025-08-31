
#
#
# Frequentist Models ------------------------------------------------------
#
#

# Accuracy Summary Statistics ------------------------------------------------------

# Calculate summaries
e3_test_dat %>%
  group_by(sid, condition) %>%
  summarize(yes=mean(said_yes), n=n()) -> e3_test_sum

# Create a wide data frame
e3_test_sum_wide = e3_test_sum %>%
  pivot_wider(id_cols=sid, names_from=condition, values_from=yes) %>%
  mutate(across(
    matches("^(aloud|silent|control)"),  # Apply to all hit rate columns
    ~ qnorm(d_prime_correct(., 40)) - qnorm(d_prime_correct(new, 80)),  # Compute d'
    .names = "d_{.col}"  # Name new columns with "d_"
  ))

# Create a long d' data frame
e3_test_d_sum = e3_test_sum_wide %>%
  select(sid, starts_with('d_')) %>%
  pivot_longer(cols=c(d_aloud:d_silent), names_to='condition', values_to='d') %>%
  mutate(condition = str_remove(condition, "^d_"))

# Preprocess data for DPSDT
e3_dpsdt_dat = preprocess_dpsdt(e3_test_dat)

# Combining into a single test phrase
e3_test_sum %>%
  left_join(e3_dpsdt_dat) %>%
  left_join(e3_test_d_sum) -> e3_test_sum_total

# ANOVAs --------------------------------------------------------

# Analysis of the "hits"
e3_anova_acc = ezANOVA(e3_test_sum_total, .(yes), .(sid), .(condition))
e3_anova_acc$Descriptives = ezStats(e3_test_sum_total, .(yes), .(sid), .(condition))  %>%
  mutate(SE = SD / sqrt(N))

# Analysis of the d'
e3_anova_d = ezANOVA(e3_test_sum_total %>% filter(condition != 'new'), .(d), .(sid), .(condition))
e3_anova_d$Descriptives = ezStats(e3_test_sum_total %>% filter(condition != 'new'), .(d), .(sid), .(condition))

# Analysis of the Rho
e3_anova_r = ezANOVA(e3_test_sum_total %>% filter(condition != 'new'), .(rho), .(sid), .(condition))
e3_anova_r$Descriptives = ezStats(e3_test_sum_total %>% filter(condition != 'new'), .(rho), .(sid), .(condition))

# Analysis of the Familiarity
e3_anova_f = ezANOVA(e3_test_sum_total %>% filter(condition != 'new'), .(df), .(sid), .(condition))
e3_anova_f$Descriptives = ezStats(e3_test_sum_total %>% filter(condition != 'new'), .(df), .(sid), .(condition))

#
#
# Bayesian Models ------------------------------------------------------
#
#

# Basic Accuracy Model
temp = e3_test_dat %>% mutate(condition=relevel(factor(condition), 'new'))
e3_bayes_m1 = brm(said_yes ~ condition + (condition | sid) + (condition | word),
                  family=brms::bernoulli(link='probit'),
                  backend = 'cmdstan',
                  cores = ncores_brm,
                  chains = nchains_brm,
                  prior = e1_acc_priors,
                  iter = niter_brm,
                  sample_prior = 'yes',
                  threads = threading(4),
                  control = list(adapt_delta = .95),
                  file = 'data/E3/models/e3_bayes_m1',
                  file_refit = 'on_change',
                  data=temp)

e3_bayes_m1_means = emmeans(e3_bayes_m1, ~condition)
e3_bayes_m1_comps = pairs(e3_bayes_m1_means)

# Analysis of the Rho
temp = e3_test_sum_total %>% filter(condition != 'new')
e3_bayes_rho1 = brm(rho ~ condition-1 + (1 | sid),
                    backend = 'cmdstan',
                    cores = ncores_brm,
                    chains = nchains_brm,
                    prior = e2_rho_priors,
                    iter = niter_brm,
                    sample_prior = 'yes',
                    threads = threading(4),
                    control = list(adapt_delta = .999, max_treedepth=20),
                    file = 'data/E3/models/e3_bayes_rho1',
                    file_refit = 'on_change',
                    data=temp)

e3_bayes_rho1_means = emmeans(e3_bayes_rho1, ~condition)
e3_bayes_rho1_comps = pairs(e3_bayes_rho1_means)

# Analysis of the Familiarity
temp = e3_test_sum_total %>% filter(condition != 'new')
e3_bayes_f1 = brm(df ~ condition-1 + (1 | sid),
                  backend = 'cmdstan',
                  cores = ncores_brm,
                  chains = nchains_brm,
                  prior = e2_f_priors,
                  iter = niter_brm,
                  sample_prior = 'yes',
                  threads = threading(4),
                  control = list(adapt_delta = .999, max_treedepth=20),
                  file = 'data/E3/models/e3_bayes_f1',
                  file_refit = 'on_change',
                  data=temp)

e3_bayes_f1_means = emmeans(e3_bayes_f1, ~condition)
e3_bayes_f1_comps = pairs(e3_bayes_f1_means)

# Output Summaries --------------------------------------------------------

write.xlsx(e3_test_sum_wide, 'data/E3/summaries/e3_test_sum_total.xlsx')
