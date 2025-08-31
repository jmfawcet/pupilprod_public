
# Get the file names ------------------------------------------------------

files = list.files(
  path = 'data/E4/pupil'
  , pattern = '.edf'
  , full.names = T
  , recursive = T
)

# Read in the data --------------------------------------------------------

overwrite_rds = FALSE
e4_study_pd = NULL
for(edf_path in files)
{
  sid = str_split(edf_path, "/", simplify = TRUE)[4]
  
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
    trial_cond = extract_trial_metadata_e4(events)
    
    # Grab the samples
    message("    Accessing samples...")
    pupil_data <- edf_data$samples %>%
      mutate(paR = ifelse(eye=='RIGHT', paR, paL),
             gxR = ifelse(eye=='RIGHT', gxR, gxL),
             gyR = ifelse(eye=='RIGHT', gyR, gyL)) %>%
      select(trial, time, eye, paR, gxR, gyR) %>%
      rename(rps = paR,
             x = gxR, y = gyR) %>%
      mutate(sid = sid)
      
    # Run the basic processing steps
    pupil_data = pupil_data %>%
      { message("    Marking blinks..."); identity(mark_blinks(., edf_data$blinks)) } %>%
      { message("    Detecting dropouts..."); identity(detect_dropouts(., sampling_rate=sampling_rate)) } %>%
      { message("    Detecting artifacts..."); identity(detect_artifacts(., sampling_rate=sampling_rate)) } %>%
      { message("    Marking data missing due to gaze..."); identity(missing_gaze(.)) } %>%
      { message("    Marking unreliable data..."); identity(mark_unreliable_data(.)) } %>%
      { message("    Marking bad trials (missing data)..."); identity(remove_excessive_missing(.)) } %>%
      { message("    Interpolating missing data..."); identity(interpolate_pupil(.)) } %>% 
      { message("    Apply a 4 hz low-pass filter..."); identity(lowpass_filter(., sampling_rate=sampling_rate)) } %>%
      { message("    Combining meta data..."); identity(left_join(., trial_cond)) } %>%
      { message("    Epoching data..."); identity(epoch_data(., lower_end=-200, upper_end=6000)) } %>%
      { message("    Baselining data..."); identity(baseline_pupil(.)) } %>%
      { message("    Downsampling to 50Hz..."); identity(downsample_pupil(.)) } %>%
      { message("    Flag bad trials..."); identity(bad_trials_by_missing_data(., missing_threshold_per_bin=.8)) } %>%
      { message("    Creating normalized values..."); identity(normalize_pupil_zscore(.)) }
    
    message("    Plotting study phase trials...")
  
    g1 = pupil_data %>% filter(phase=='study') %>% ggplot(aes(x = bin, y = pupil_z)) +
      facet_wrap(~ trial) +
      geom_point(aes(color = pmissing)) +
      scale_color_gradient(low = "green", high = "red") +
      theme_minimal()
    
    ggsave(g1, file=sprintf("data/E4/pupil_plots/%s_study.pdf", sid), width=12, height=12)
  
    saveRDS(pupil_data, str_replace(edf_path, 'edf', 'rds'))
  }
  
  e4_study_pd = bind_rows(e4_study_pd, pupil_data %>% filter(phase=='study'))
}

# Check Bad Trials --------------------------------------------------------

e4_study_pd %>%
  group_by(sid, trial) %>%
  summarize(bad_trials = mean(bad_trial), pmissing = mean(pmissing)) %>%
  group_by(sid) %>%
  summarize(bad_trials = mean(bad_trials), pmissing = mean(pmissing)) -> e4_bad_trial_summary


# Exclude Bad Participants ------------------------------------------------

# No one with too many bad trials
e4_too_many_bad_trials = e4_bad_trial_summary %>% filter(bad_trials > bad_trial_threshold_per_part) %>% pull(sid)

# pilot1 was not a real participant, reflecting instead the lead RA running themselves, no EDF was recorded
# PP3003 had a technical (user) error where the study phase was started unexpectedly when trying to calibrate, missing a handful of trials while trying to regroup
# PP3014 was hard to calibrate due to glasses and had to leave before the task completed, no data file
# PP3041 is missing the EDF file - this happened the week before the lockdown and the reason was not recorded, excluded behavioural for consistency
e4_flagged_participants = c(
  'pilot1', 'pp3003', 'pp3014', 'pp3041'
)

e4_bad_part = c(e4_too_many_bad_trials, e4_flagged_participants)

# Exclude Bad Trials ------------------------------------------------------

e4_study_pd = e4_study_pd %>%
  filter(!bad_trial, !(sid %in% e4_bad_part)) %>%
  mutate(condition = factor(condition), sid = factor(sid), word = factor(word))

e4_bad_trial_summary %>% filter(!(sid %in% e4_bad_part)) %>% summarize_if(is.numeric, list(mean=mean, sd=sd))

e4_study_pd %>%
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

e4_test_dat_dir = 'data/E4/behavioural'

# Read-in data
e4_test_dat = NULL
for(l in list.files(e4_test_dat_dir, pattern='Test_Result', recursive=TRUE, full.names = TRUE))
{
  print(l)
  dat = read.table(l, header=TRUE) %>%
    mutate(sid=extract_participant_id(l), exp='E4') %>%
    select(sid, exp, condition=test_condition, word=test_word, key=Key_Pressed) %>%
    mutate(said_yes = as.numeric(key > 3),
           condition = recode(condition, 'prod'='aloud', 'check'='control'),
           word = tolower(word))
  
  e4_test_dat = rbind(e4_test_dat, dat)
}

# Exclude the same test phase participants as pupil participants
e4_test_dat %>%
  filter(!(sid %in% e4_bad_part)) -> e4_test_dat

#
#
# Combine Pupil and Behavioural -------------------------------------------
#
#

e4_study_pd = e4_study_pd %>%
  mutate(sid=as.character(sid), condition=as.character(condition), word=as.character(word)) %>%
  left_join(e4_test_dat) %>%
  mutate(mem = ifelse(said_yes==0, 'no', 'yes'), condmem = paste(condition, mem))


# Create a summary data frame
e4_study_pd_sum = e4_study_pd %>%
  group_by(sid, bin, condition) %>%
  summarize(pupil_z = mean(pupil_z))
