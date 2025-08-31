
#
#
# Frequentist Models ------------------------------------------------------
#
#

# Accuracy Summary Statistics ------------------------------------------------------

# Calculate summaries
e4_test_dat %>%
  group_by(sid, condition) %>%
  summarize(yes=mean(said_yes), n=n()) -> e4_test_sum

# Create a wide data frame
e4_test_sum_wide = e4_test_sum %>%
  pivot_wider(id_cols=sid, names_from=condition, values_from=yes) %>%
  mutate(across(
    matches("^(aloud|silent|control)"),  # Apply to all hit rate columns
    ~ qnorm(d_prime_correct(., 40)) - qnorm(d_prime_correct(new, 80)),  # Compute d'
    .names = "d_{.col}"  # Name new columns with "d_"
  ))

# Create a long d' data frame
e4_test_d_sum = e4_test_sum_wide %>%
  select(sid, starts_with('d_')) %>%
  pivot_longer(cols=c(d_aloud:d_silent), names_to='condition', values_to='d') %>%
  mutate(condition = str_remove(condition, "^d_"))

# Preprocess data for DPSDT
e4_dpsdt_dat = preprocess_dpsdt(e4_test_dat)

# Combining into a single test phrase
e4_test_sum %>%
  left_join(e4_dpsdt_dat) %>%
  left_join(e4_test_d_sum) -> e4_test_sum_total

# ANOVAs --------------------------------------------------------

# Analysis of the "hits"
e4_anova_acc = ezANOVA(e4_test_sum_total, .(yes), .(sid), .(condition))
e4_anova_acc$Descriptives = ezStats(e4_test_sum_total, .(yes), .(sid), .(condition))  %>%
  mutate(SE = SD / sqrt(N))

# Analysis of the d'
e4_anova_d = ezANOVA(e4_test_sum_total %>% filter(condition != 'new'), .(d), .(sid), .(condition))
e4_anova_d$Descriptives = ezStats(e4_test_sum_total %>% filter(condition != 'new'), .(d), .(sid), .(condition))

# Analysis of the Rho
e4_anova_r = ezANOVA(e4_test_sum_total %>% filter(condition != 'new'), .(rho), .(sid), .(condition))
e4_anova_r$Descriptives = ezStats(e4_test_sum_total %>% filter(condition != 'new'), .(rho), .(sid), .(condition))

# Analysis of the Familiarity
e4_anova_f = ezANOVA(e4_test_sum_total %>% filter(condition != 'new'), .(df), .(sid), .(condition))
e4_anova_f$Descriptives = ezStats(e4_test_sum_total %>% filter(condition != 'new'), .(df), .(sid), .(condition))

#
#
# Bayesian Models ------------------------------------------------------
#
#

# Basic Accuracy Model
temp = e4_test_dat %>% mutate(condition=relevel(factor(condition), 'new'))
e4_bayes_m1 = brm(said_yes ~ condition + (condition | sid) + (condition | word),
                  family=brms::bernoulli(link='probit'),
                  backend = 'cmdstan',
                  cores = ncores_brm,
                  chains = nchains_brm,
                  prior = e1_acc_priors,
                  iter = niter_brm,
                  sample_prior = 'yes',
                  threads = threading(4),
                  control = list(adapt_delta = .95),
                  file = 'data/E4/models/e4_bayes_m1',
                  file_refit = 'on_change',
                  data=temp)

e4_bayes_m1_means = emmeans(e4_bayes_m1, ~condition)
e4_bayes_m1_comps = pairs(e4_bayes_m1_means)

# Analysis of the Rho
temp = e4_test_sum_total %>% filter(condition != 'new')
e4_bayes_rho1 = brm(rho ~ condition-1 + (1 | sid),
                    backend = 'cmdstan',
                    cores = ncores_brm,
                    chains = nchains_brm,
                    prior = e2_rho_priors,
                    iter = niter_brm,
                    sample_prior = 'yes',
                    threads = threading(4),
                    control = list(adapt_delta = .999, max_treedepth=20),
                    file = 'data/E4/models/e4_bayes_rho1',
                    file_refit = 'on_change',
                    data=temp)

e4_bayes_rho1_means = emmeans(e4_bayes_rho1, ~condition)
e4_bayes_rho1_comps = pairs(e4_bayes_rho1_means)

# Analysis of the Familiarity
temp = e4_test_sum_total %>% filter(condition != 'new')
e4_bayes_f1 = brm(df ~ condition-1 + (1 | sid),
                  backend = 'cmdstan',
                  cores = ncores_brm,
                  chains = nchains_brm,
                  prior = e2_f_priors,
                  iter = niter_brm,
                  sample_prior = 'yes',
                  threads = threading(4),
                  control = list(adapt_delta = .999, max_treedepth=20),
                  file = 'data/E4/models/e4_bayes_f1',
                  file_refit = 'on_change',
                  data=temp)

e4_bayes_f1_means = emmeans(e4_bayes_f1, ~condition)
e4_bayes_f1_comps = pairs(e4_bayes_f1_means)

# Output Summaries --------------------------------------------------------

write.xlsx(e4_test_sum_wide, 'data/E4/summaries/e4_test_sum_total.xlsx')
