#
#
# Read in Pupillometry Data -----------------------------------------------
#
#

# Get the Pupil file names ------------------------------------------------

files = list.files(
  path = 'data/E1/pupil/'
  , pattern = '.edf'
  , full.names = T
  , recursive = T
)

# Read in the data --------------------------------------------------------

overwrite_rds = FALSE
e1_study_pd = NULL
for(edf_path in files)
{
  sid = str_remove(str_split(edf_path, "/", simplify = TRUE)[5], ".edf")
  
  rds_file = str_replace(edf_path, 'edf', 'rds')
  if(file.exists(rds_file) & !overwrite_rds)
  {
    message(sprintf("File Exists. Reading file %s", rds_file))
    pupil_data = readRDS(str_replace(edf_path, 'edf', 'rds'))
  } else {
    
    # Read in the data
    message(sprintf("Reading file %s", edf_path))
    edf_data <- eyelinkReader::read_edf(edf_path, import_samples=TRUE)
    
    # Store current sampling rate
    sampling_rate = edf_data$recordings$sample_rate[1]
    
    # Store eye measured
    sampled_eye = edf_data$recordings$eye[1]
    message(sprintf("    Measuring %s Eye", sampled_eye))
    
    # Store Events
    events = edf_data$events %>%
      select(trial, time, sttime, entime, sttime_rel, entime_rel, message)
    
    # Create trial conditions
    message("    Reading Trial Meta Data...")
    trial_cond = extract_trial_metadata_e1(events)
    
    # Grab the samples
    message("    Accessing samples...")
    pupil_data <- edf_data$samples %>%
      mutate(paR = ifelse(eye=='RIGHT', paR, paL),
             gxR = ifelse(eye=='RIGHT', gxR, gxL),
             gyR = ifelse(eye=='RIGHT', gyR, gyL)) %>%
      select(trial, time, eye, paR, gxR, gyR) %>%
      rename(rps = paR,
             x = gxR, y = gyR) %>%
      mutate(sid = sid) %>%
      group_by(trial, time, eye, sid) %>%
      summarize(rps=mean(rps, na.rm=TRUE), x=mean(x, na.rm=TRUE), y=mean(y, na.rm=TRUE)) %>%
      ungroup() %>%
      select(trial, time, eye, rps, x, y, sid) %>%
      ungroup() %>%
      arrange(sid, eye, time, trial) %>%  # Ensure time is ordered properly
      distinct(sid, eye, time, .keep_all = TRUE)  # Keep first unique time per eye/sid
    
    # NOTE:
    # By averaging above, we drop from 2000 Hz to 1000 Hz
    sampling_rate = 1000
      
    # Run the basic processing steps
    pupil_data = pupil_data %>%
      { message("    Marking blinks..."); identity(mark_blinks(., edf_data$blinks)) } %>%
      { message("    Detecting dropouts..."); identity(detect_dropouts(., sampling_rate=sampling_rate)) } %>%
      { message("    Detecting artifacts..."); identity(detect_artifacts(., sampling_rate=sampling_rate)) } %>%
      { message("    Marking data missing due to gaze..."); identity(missing_gaze(., lower_x=600, upper_x=1200, lower_y=150, upper_y=1300)) } %>%
      { message("    Marking unreliable data..."); identity(mark_unreliable_data(.)) } %>%
      { message("    Marking bad trials (missing data)..."); identity(remove_excessive_missing(.)) } %>%
      { message("    Interpolating missing data..."); identity(interpolate_pupil(.)) } %>% 
      { message("    Apply a 4 hz low-pass filter..."); identity(lowpass_filter(., sampling_rate=sampling_rate)) } %>%
      { message("    Combining meta data..."); identity(left_join(., trial_cond)) } %>%
      { message("    Epoching data..."); identity(epoch_data(., lower_end=-200, upper_end=2000)) } %>%
      { message("    Baselining data..."); identity(baseline_pupil(.)) } %>%
      { message("    Downsampling to 50Hz..."); identity(downsample_pupil(.)) } %>%
      { message("    Flag bad trials..."); identity(bad_trials_by_missing_data(., missing_threshold_per_bin = 0.8)) } %>%
      { message("    Creating normalized values..."); identity(normalize_pupil_zscore(.)) }
    
    message("    Plotting study phase trials...")
  
    g1 = pupil_data %>% filter(phase=='study') %>% ggplot(aes(x = bin, y = pupil_z)) +
      facet_wrap(~ trial) +
      geom_point(aes(color = pmissing)) +
      scale_color_gradient(low = "green", high = "red") +
      theme_minimal()
    
    ggsave(g1, file=sprintf("data/E1/pupil_plots/%s_study.pdf", sid), width=12, height=12)
    
    saveRDS(pupil_data, str_replace(edf_path, 'edf', 'rds'))
  }
  
  e1_study_pd = bind_rows(e1_study_pd, pupil_data %>% filter(phase=='study'))
}

# Check Bad Trials --------------------------------------------------------

e1_study_pd %>%
  group_by(sid, trial) %>%
  summarize(bad_trials = mean(bad_trial), pmissing = mean(pmissing)) %>%
  group_by(sid) %>%
  summarize(bad_trials = mean(bad_trials), pmissing = mean(pmissing)) -> e1_bad_trial_summary

# Exclude Bad Participants ------------------------------------------------

# Use a gentler exclusion % used here (50%)
# 1 = 100% bad
# 5 = 59% bad
# 13 = 65% bad
# 14 = 68% bad
# 16 = 81% bad
# 20 = 94% bad
# 24 = 58% bad
# 27 = 74% bad
# 32 = 59% bad -> Fixed with 80% change
# 39 = 100% bad
# 43 = 66% bad
# 44 = 93% bad
# 48 = 62% bad
# 49 = 54% bad
e1_too_many_bad_trials = e1_bad_trial_summary %>% filter(bad_trials > .5) %>% pull(sid)

e1_flagged_participants = c(
  # No one flagged other than above
)

e1_bad_part = c(e1_too_many_bad_trials, e1_flagged_participants)

# Exclude Bad Trials ------------------------------------------------------

e1_study_pd = e1_study_pd %>%
  filter(!bad_trial, !(sid %in% e1_bad_part)) %>%
  mutate(condition = factor(condition), sid = factor(sid), word = factor(word))

e1_bad_trial_summary %>% filter(!(sid %in% e1_bad_part)) %>% summarize_if(is.numeric, list(mean=mean, sd=sd))

e1_study_pd %>%
  group_by(sid, trial) %>%
  summarize(bad_trials = mean(bad_trial), pmissing = mean(pmissing)) %>%
  group_by(sid) %>%
  summarize(bad_trials = mean(bad_trials), pmissing = mean(pmissing)) %>%
  summarize(m=mean(pmissing), sd=sd(pmissing))

#
#
# Read in Behavioural Data ------------------------------------------------
#
#

# Get the Behavioural file names ------------------------------------------

e1_test_dat_dir = 'data/E1/behavioural'

# Read-in data
e1_test_dat = NULL
for(l in list.files(e1_test_dat_dir, pattern='RESULTS_FILE', recursive=TRUE, full.names = TRUE))
{
  print(l)
  sdat = read.table(l, header=TRUE) %>%
    filter(TRIAL_INDEX_RECOG == 0) %>%
    mutate(sid=extract_participant_id(l)) %>%
    select(sid, word=TRIAL_WORD_ENCODING, condition=TRIAL_CONDITION_ENCODING) %>%
    mutate(word=str_to_lower(word), condition=str_to_lower(condition))
  
  dat = read.table(l, header=TRUE) %>%
    filter(TRIAL_INDEX_RECOG > 0) %>%
    mutate(sid=extract_participant_id(l), exp='E1', said_yes = as.numeric(correctanswer_1==KEYPRESS_RECOG)) %>%
    mutate(said_yes = ifelse(condition_1=='old', said_yes, 1-said_yes)) %>%
    select(sid, exp, word= wordrecog_1, said_yes) %>%
    mutate(word = str_to_lower(word)) %>%
    left_join(sdat) %>%
    mutate(condition = ifelse(is.na(condition), 'new', condition), key=-1) %>%
    select(sid, exp, condition, word, key, said_yes)

  e1_test_dat = rbind(e1_test_dat, dat)
}

# Exclude the same test phase participants as pupil participants
e1_test_dat %>%
  filter(!(sid %in% e1_bad_part)) -> e1_test_dat

#
#
# Combine Pupil and Behavioural -------------------------------------------
#
#

e1_study_pd = e1_study_pd %>%
  mutate(sid=as.character(sid), condition=as.character(condition), word=as.character(word)) %>%
  left_join(e1_test_dat) %>%
  mutate(mem = ifelse(said_yes==0, 'no', 'yes'), condmem = paste(condition, mem))
  
# Create a summary data frame
e1_study_pd_sum = e1_study_pd %>%
  group_by(sid, bin, condition) %>%
  summarize(pupil_z = mean(pupil_z))
