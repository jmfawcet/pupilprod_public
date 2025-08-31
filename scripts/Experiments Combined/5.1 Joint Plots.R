
# Behavioural Plot --------------------------------------------------------

e1_behav_plot_dat = prepare_behav_dat_figure_e1(e1_bayes_m1_means, e1_bayes_m1_comps, exp='Experiment 1')
e2_behav_plot_dat = prepare_behav_dat_figure(e2_bayes_m1_means, e2_bayes_m1_comps, 
                                             e2_bayes_rho1_means, e2_bayes_rho1_comps,
                                             e2_bayes_f1_means, e2_bayes_f1_comps, exp='Experiment 2')
e3_behav_plot_dat = prepare_behav_dat_figure(e3_bayes_m1_means, e3_bayes_m1_comps, 
                                             e3_bayes_rho1_means, e3_bayes_rho1_comps,
                                             e3_bayes_f1_means, e3_bayes_f1_comps, exp='Experiment 3')
e4_behav_plot_dat = prepare_behav_dat_figure(e4_bayes_m1_means, e4_bayes_m1_comps, 
                                             e4_bayes_rho1_means, e4_bayes_rho1_comps,
                                             e4_bayes_f1_means, e4_bayes_f1_comps, exp='Experiment 4')
ecomb_behav_plot_dat = e1_behav_plot_dat %>% bind_rows(e2_behav_plot_dat, e3_behav_plot_dat, e4_behav_plot_dat)

notes_frame <- data.frame(
  experiment      = c("Experiment 1", "Experiment 1", 
               "Experiment 1", "Experiment 1"),  # matches facet levels
  type     = c('Rec.', 'Rec. Contrasts', 'Fam.', 'Fam. Contrasts'),
  condition        = c("control", "aloud - control", "control", "aloud - control"),        # where arrow starts horizontally
  y        = c(rep(.5, 4)),        # arrow's start
  label    = c(rep('Not Measured', 4))
) %>%
  mutate(experiment = factor(experiment, levels=c('Experiment 1', 'Experiment 2', 'Experiment 3', 'Experiment 4')),
         condition = factor(condition, levels = c("aloud", "control", "silent", "aloud - silent", "aloud - control", "control - silent")))

ecomb_behav_plot_dat = ecomb_behav_plot_dat %>%
  mutate(
    condition_group = case_when(
      condition %in% c("aloud", "silent", "control") ~ condition,
      condition %in% c("aloud - silent", "aloud - control", "control - silent") ~ "contrast"
    ), condition_group = factor(condition_group, levels=c('aloud', 'control', 'silent', 'contrast'))
  )

g1 = behav_dat_plot(ecomb_behav_plot_dat) + geom_text(
  data = notes_frame,
  aes(x = condition, y = y, label = label),
  inherit.aes = FALSE,
  vjust = -0.5,       # Nudges text above the arrow's start
  hjust = 0.5,      # Center horizontally
  size  = 3,        # Font size
  color = "black"
) + theme(legend.background = element_blank())


ggsave('figures/Figure2_BehaviouralData.pdf', g1, width=10, height=10, bg='transparent')
ggsave('figures/Figure2_BehaviouralData.tiff', g1, width=10, height=10, bg='transparent')

# Pupil Plots -------------------------------------------------------------

e1_study_pgam_plot = plot_smooth(e1_study_pgam_ww, view = "bin", 
                                  main = "Estimated Pupil Response Over Time", 
                                  rug = TRUE, 
                                  plot_all = 'condition')$fv %>%
  select(condition, bin, fit, CI, ll, ul) %>%
  mutate(exp='Experiment 1', model='GAMM')

e1_study_pd_sum_figure = e1_study_pd %>%
  group_by(bin, condition) %>%
  summarize(m = mean(pupil_z), se = sd(pupil_z)/n()**.5, ci=se*1.96) %>%
  mutate(model='Mass Univariate') %>%
  bind_rows(
    e1_study_pgam_plot %>% 
      mutate(se = CI/1.96) %>%
      select(bin, condition, m=fit, se=se, ci=CI, model)
  ) %>%
  mutate(model = factor(model, levels=c('Mass Univariate', 'GAMM')), exp='Experiment 1')

e1_study_pd_sum_figure = e1_study_pd_sum_figure
e2_study_pd_sum_figure = create_study_pd_sum_figure(e2_study_pd, e2_study_pgam_ww, exp='Experiment 2')
e3_study_pd_sum_figure = create_study_pd_sum_figure(e3_study_pd, e3_study_pgam_ww, exp='Experiment 3')
e4_study_pd_sum_figure = create_study_pd_sum_figure(e4_study_pd, e4_study_pgam_ww, exp='Experiment 4')

e1_gamm_plot_sig = plot_diff(e1_study_pgam_ww, view='bin', 
                              comp=list(condition=c('aloud', 'silent')),
                              n.grid=1000, plot=FALSE) %>%
  select(est, bin, CI) %>%
  mutate(comparison = 'Aloud - Silent', sig = est + CI < 0 | est - CI > 0) %>%
  filter(sig) %>%
  mutate(y_offset = case_when(
    comparison == "Aloud - Silent" ~ -0.6,
    comparison == "Aloud - Control" ~ -0.65,
    comparison == "Control - Silent" ~ -0.7
  ), exp='Experiment 1', model='GAMM') %>%
  select(bin, comparison, y_offset, exp, model)

e1_plot_sig = e1_massu %>%
  filter(as_clust) %>%
  mutate(y_offset = case_when(
    comparison == "Aloud - Silent" ~ -0.6,
    comparison == "Aloud - Control" ~ -0.65,
    comparison == "Control - Silent" ~ -0.7
  ), exp='Experiment 1', model='Mass Univariate') %>%
  select(bin, comparison, y_offset, exp, model) %>%
  rownames_to_column("original_names") %>%
  mutate(original_names = row_number()) %>%
  column_to_rownames("original_names") %>%
  bind_rows(
    e1_gamm_plot_sig
  ) %>%
  mutate(model = factor(model, levels=c('Mass Univariate', 'GAMM')))

e1_plot_sig = e1_plot_sig
e2_plot_sig = create_plot_sig_cluster(e2_study_pgam_ww, e2_massu, exp='Experiment 2')
e3_plot_sig = create_plot_sig_cluster(e3_study_pgam_ww, e3_massu, exp='Experiment 3')
e4_plot_sig = create_plot_sig_cluster(e4_study_pgam_ww, e4_massu, exp='Experiment 4')

ecomb_study_pd_sum_figure = e1_study_pd_sum_figure %>%
  bind_rows(e2_study_pd_sum_figure,
            e3_study_pd_sum_figure,
            e4_study_pd_sum_figure) %>%
  mutate(exp = factor(exp, levels=c('Experiment 1', 'Experiment 2', 'Experiment 3', 'Experiment 4')))

ecomb_plot_sig = e1_plot_sig %>%
  bind_rows(e2_plot_sig,
            e3_plot_sig,
            e4_plot_sig) %>%
  mutate(y_offset = case_when(
    y_offset == -.6 ~ -.63,
    y_offset == -.65 ~ -.7,
    y_offset == -.7 ~ -.77
  )) %>%
  mutate(exp = factor(exp, levels=c('Experiment 1', 'Experiment 2', 'Experiment 3', 'Experiment 4')))

arrow_data2 <- data.frame(
  exp      = c("Experiment 1", "Experiment 1", 
               "Experiment 2", "Experiment 2", 
               'Experiment 3', 'Experiment 3', 'Experiment 3', 
               'Experiment 4', 'Experiment 4', 'Experiment 4'),  # matches facet levels
  x        = c(0, 1500, 0, 1500, 0, 1500, 3000, 0, 3500, 5000),        # where arrow starts horizontally
  xend     = c(0, 1500, 0, 1500, 0, 1500, 3000, 0, 3500, 5000),        # the same x, so it's vertical
  y        = c(0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7),        # arrow's start
  yend     = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),         # arrow's end (below start => downward arrow)
  label    = c('Word\nInstruction\nResponse', 'E', 'Word\nInstruction\nResponse', "E", "Instruction", "Word\nResponse", "E", "Word\nInstruction", "Response", "E")
) %>%
  mutate(exp = factor(exp, levels=c('Experiment 1', 'Experiment 2', 'Experiment 3', 'Experiment 4'))) %>%
  filter(label!='E') %>%
  mutate(vjust = -.35 + .125*str_count(label, '\n'))

# Define labels and y-offsets
label_positions <- data.frame(
  comparison = c("aloud - silent", "aloud - control", "control - silent"),
  y_label_offset = c(-0.64, -0.71, -0.78)  # Slightly above segments for clarity
)

g1 = ecomb_study_pd_sum_figure %>%
  mutate() %>% 
  ggplot(aes(x = bin, y = m, fill=condition, color = condition, linetype = condition, ymin=m-ci, ymax=m+ci)) +
  geom_line() +
  geom_ribbon(alpha=.5) +
  geom_segment(data = ecomb_plot_sig,
               aes(x = bin, xend = bin, y = y_offset, yend = y_offset - 0.02),
               inherit.aes = FALSE, size = 1) +
  geom_text(data = label_positions,
            aes(x = min(e1_massu$bin) + 600,  # Align all labels at the same x position
                y = y_label_offset, 
                label = comparison),
            inherit.aes = FALSE, hjust = 1, size = 2.5) +
  labs(x = "Time (ms)", y = "z(Pupil)") +
  theme_classic() +
  geom_hline(
    yintercept = 0,
    colour='grey',
    alpha=.2
  ) + 
  facet_grid(exp~model) +
  theme(
    # Make facet strips black with white text
    strip.background = element_rect(fill = "black", color = NA),
    strip.text.x     = element_text(color = "white"),
    strip.text.y     = element_text(color = "white"),
    
    
    # Transparent panel and plot background
    panel.background = element_rect(fill = NA, color = NA),
    plot.background  = element_rect(fill = NA, color = NA),
    
    # Make axis & legend text visible on transparent background
    axis.text        = element_text(color = "black"),
    axis.title       = element_text(color = "black"),
    #legend.text      = element_text(color = "white"),
    #legend.title     = element_text(color = "white")
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.background = element_blank(),
    legend.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  geom_segment(
    data = arrow_data2,
    aes(x = x, xend = xend, y = y, yend = yend),
    inherit.aes = FALSE,
    # Add arrow at the "last" end of the segment (or "first", "both")
    arrow = arrow(
      length = unit(0.2, "cm"),   # arrow size
      ends = "last",             # which end to place arrow on
      type = "closed"            # closed arrow head
    ),
    color = "black",
    size  = 1
  ) +
  geom_text(
    data = arrow_data2,
    aes(x = x, y = y, label = label, vjust=vjust),
    inherit.aes = FALSE,
    #vjust = -.2,       # Nudges text above the arrow's start
    hjust = 0.5,      # Center horizontally
    size  = 2.5,        # Font size
    color = "black",
    lineheight = 0.9
  ) +
  coord_cartesian(clip = "off", ylim=c(-.8,.8)) +
  scale_color_manual(values = c("aloud" = "#F8766D", "control" = "#7CAE00", "silent" = "#00BFC4"))  # Customize colors

ggsave('figures/Figure3_PupillometricData.pdf', g1, width=12, height=14, bg='transparent')
ggsave('figures/Figure3_PupillometricData.tiff', g1, width=12, height=14, bg='transparent')
