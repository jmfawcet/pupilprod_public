
#
#
# Frequentist Models ------------------------------------------------------
#
#

# Accuracy Summary Statistics ------------------------------------------------------

# Calculate summaries
e1_test_dat %>%
  group_by(sid, condition) %>%
  summarize(yes=mean(said_yes), n=n()) -> e1_test_sum

# Create a wide data frame
e1_test_sum_wide = e1_test_sum %>%
  pivot_wider(id_cols=sid, names_from=condition, values_from=yes) %>%
  mutate(across(
    matches("^(aloud|silent|control)"),  # Apply to all hit rate columns
    ~ qnorm(d_prime_correct(., 40)) - qnorm(d_prime_correct(new, 80)),  # Compute d'
    .names = "d_{.col}"  # Name new columns with "d_"
  ))

# Create a long d' data frame
e1_test_d_sum = e1_test_sum_wide %>%
  select(sid, starts_with('d_')) %>%
  pivot_longer(cols=c(d_aloud:d_silent), names_to='condition', values_to='d') %>%
  mutate(condition = str_remove(condition, "^d_"))
           
# Combining into a single test phrase
e1_test_sum %>%
  left_join(e1_test_d_sum) -> e1_test_sum_total
 
# ANOVAs --------------------------------------------------------

# Analysis of the "hits"
e1_anova_acc = ezANOVA(e1_test_sum_total, .(yes), .(sid), .(condition))
e1_anova_acc$Descriptives = ezStats(e1_test_sum_total, .(yes), .(sid), .(condition)) %>%
  mutate(SE = SD / sqrt(N))

# Analysis of the d'
e1_anova_d = ezANOVA(e1_test_sum_total %>% filter(condition != 'new'), .(d), .(sid), .(condition))
e1_anova_d$Descriptives = ezStats(e1_test_sum_total %>% filter(condition != 'new'), .(d), .(sid), .(condition))

#
#
# Bayesian Models ------------------------------------------------------
#
#

# Basic Accuracy Model
temp = e1_test_dat %>% mutate(condition=relevel(factor(condition), 'new'))
e1_bayes_m1 = brm(said_yes ~ condition + (condition | sid) + (condition | word),
                  family=brms::bernoulli(link='probit'),
                  backend = 'cmdstan',
                  cores = ncores_brm,
                  chains = nchains_brm,
                  prior = e1_acc_priors,
                  iter = niter_brm,
                  sample_prior = 'yes',
                  threads = threading(4),
                  control = list(adapt_delta = .95),
                  file = 'data/E1/models/e1_bayes_m1',
                  file_refit = 'on_change',
                  data=temp)

e1_bayes_m1_means = emmeans(e1_bayes_m1, ~condition)
e1_bayes_m1_comps = pairs(e1_bayes_m1_means)

# Output Summaries --------------------------------------------------------

write.xlsx(e1_test_sum_wide, 'data/E1/summaries/e1_test_sum_total.xlsx')
