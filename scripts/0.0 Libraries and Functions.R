
# Load libraries ----------------------------------------------------------

# Some of these may be legacy from older approaches, but in finalizing my script
# I was unsure and so left them all in rather than spend time trying to figure
# out which were no longer needed.
library(simts)
library(progressr)
library(pbapply)
library(broom.mixed)
library(forecast)
library(emmeans)
library(bayesCorr)
library(pracma)
library(fda)
library(MASS)
library(signal)
library(PupillometryR)
library(tidyverse)
library(rio)
library(ez)
library(mgcv)
library(openxlsx)
library(lme4)
library(brms)
library(itsadug)
library(eyetrackingR)
library(eyelinker)
library(eyelinkReader)
library(progress)
library(zoo)
library(data.table)
library(stringr)
library(parallel)
library(parallelly)
library(tidyverse)


# Color scheme ------------------------------------------------------------

aloud_col = rgb(149,155,218,max=255)
control_col = rgb(180,150,200,max=255)
silent_col = rgb(216,140,163,max=255)

condition_colors <- c("aloud" = "#F8766D", "control" = "#7CAE00", "silent" = "#00BFC4", "contrast" = "black")
condition_shapes <- c("aloud" = 15, "control" = 17, "silent" = 18, "contrast" = 16)  # Square, Triangle, Diamond, Circle

# Task properties ---------------------------------------------------------

# Define threshold values (adjust based on data distribution)
missing_threshold_per_bin <- 0.5  # Flag bins missing more than 50%
bad_bin_threshold_per_trial <- 0.3  # Flag trials where >30% of bins are bad
bad_trial_threshold_per_part <- 0.3  # Flag participants as bad if >30% bad trials

# For parallel models
# Set the number of threads for mgcv
# Use all but one core to avoid system slowdowns
options(mc.cores = parallelly::availableCores() - 1)

# BRMS parameters; this can be changed as desired
ncores_brm = 4
nchains_brm = 4
niter_brm = 6000

# Priors used for the accuracy model for E1 and E2 (also used to inform later 
# accuracy models)
e1_acc_pupilpred_priors = c(
  prior(normal(0.8, 0.5), class='b', coef='conditionaloud'),
  prior(normal(0.8, 0.5), class='b', coef='conditionsilent'),
  prior(normal(0, 0.5), class='b'),
  prior(normal(0, 1), class='sd'),
  prior(lkj(2), class = 'cor')
)

e1_acc_priors = c(
  prior(normal(-.5, 0.5), class='Intercept'),
  prior(normal(1, 1), class='b'),
  prior(normal(0, 1), class='sd'),
  prior(lkj(2), class = 'cor')
)

e2_acc_pupilpred_priors = c(
  prior(normal(0.8, 0.5), class='b', coef='conditionaloud'),
  prior(normal(0.8, 0.5), class='b', coef='conditioncontrol'),
  prior(normal(0.8, 0.5), class='b', coef='conditionsilent'),
  prior(normal(0, 0.5), class='b'),
  prior(normal(0, 1), class='sd'),
  prior(lkj(2), class = 'cor')
)

# Priors used for the rho model for E2 (used to inform later 
# accuracy models)
e2_rho_priors = c(
  prior(normal(0.3, 0.15), class='b'),
  prior(normal(0, 0.1), class='sd'),
  prior(normal(0, 0.1), class='sigma')
)

# Priors used for the f model for E2 (used to inform later 
# accuracy models)
e2_f_priors = c(
  prior(normal(1, 1), class='b'),
  prior(normal(0, 1), class='sd'),
  prior(normal(0, 1), class='sigma')
)

# Functions ---------------------------------------------------------------

#
#
# D Prime Functions -------------------------------------------------------
#
#

# d' correction function
d_prime_correct = function(x, num_items=40)
{
  if(x == 0)
  {
    x = .5/num_items
  } else if(x==1)
  {
    x = (num_items-.5)/num_items
  }
  
  return(x)
}

#
#
# Dual Process Analyses --------------------------------------------------------
#
#

# Parameter definitions ---------------------------------------------------

pars = c(rho=1, d=1, c1=1, c2=1, c3=1, c4=1, c5=1)

# Function used to fit model  ---------------------------------------------

OptPars = function(parameters, hits, fas, outputall=FALSE)
{
  predhit = parameters[1] + ((1-parameters[1]) * pnorm(parameters[c(3:7)], -parameters[2], 1))
  sq_diff_hit = (hits - predhit) ** 2
  predfa = (1-0) * pnorm(parameters[c(3:7)], 0, 1)
  sq_diff_far = (fas - predfa) ** 2
  
  if(outputall)
  {
    return(data.frame(fas=fas, hit=hits, rho=parameters[1], m=0, d=parameters[2], c=parameters[c(3:7)], predhit=predhit, sqdiffhit=sq_diff_hit, predfa=predfa, sqdifffa=sq_diff_fa))
  }
  else
  {
    return(sum(sq_diff_hit) + sum(sq_diff_far))
  }
}


# Function used to fit the model while checking for Rho  ------------------

fit_dpsdt = function(hits, fas, control=list(maxit=3000))
{
  output = optim(pars, OptPars, hits=hits, fas=fas, control=control, method='L-BFGS-B', lower=c(0, rep(-Inf, 6)))
  temp = output$par
  
  return(c(output$par, converged=output$convergence))
}

gen_sum = function(x)
{
  temp = NULL
  for(i in 1:5)
  {
    temp = c(temp, mean(x > i))
  }
  
  return(temp)
}

preprocess_dpsdt <- function(data) {
  temp_dat <- NULL
  
  for (i in unique(data$sid)) {
    c_dat <- data %>% filter(sid == i)
    fa_curve <- gen_sum(c_dat %>% filter(condition == 'new') %>% pull(key))
    
    for (j in c('silent', 'aloud', 'control')) {
      c_curve <- gen_sum(c_dat %>% filter(condition == j) %>% pull(key))
      c_fit <- fit_dpsdt(c_curve, fa_curve)
      
      temp_dat <- rbind(temp_dat, data.frame(sid = i, condition = j, rho = c_fit[[1]], df = c_fit[[2]]))
    }
  }
  
  return(temp_dat)
}


# Helper to grab data files -----------------------------------------------

extract_participant_id <- function(file_path) {
  file_name <- basename(file_path)  # Get the file name from the path
  participant_id <- sub("_.*", "", file_name)  # Extract the portion before the first underscore
  return(participant_id)
}

# Extract Trial Metadata --------------------------------------------------------

extract_trial_metadata_e1 <- function(events) {
  
  # Get study phase indices
  fixation_onset = events %>% filter(str_detect(message, 'FIXATIONTIME$')) 
  words_onset = events %>% filter(str_detect(message, 'word_encoding')) %>%
    select(trial, word=message) %>%
    mutate(word = str_to_lower(str_remove(word, "^!V TRIAL_VAR word_encoding\\s*")))
  cond_onset = events %>% filter(str_detect(message, 'condition ')) %>%
    select(trial, condition=message) %>%
    mutate(condition = str_to_lower(str_remove(condition, "^!V TRIAL_VAR condition \\s*")))
  
  t_dat = data.frame(trial=fixation_onset$trial, 
                     phase='study',
                     word_onset=fixation_onset$sttime) %>%
    left_join(words_onset) %>%
    left_join(cond_onset) %>%
    select(trial, phase, word, condition, word_onset)
  
  return(t_dat)
}

extract_trial_metadata_e2 <- function(events) {
  
  # Get study phase indices
  events2 = events %>% filter(entime < 10, message != "", !str_detect(message, "TRACKER_TIME|TRIAL_VAR"))
  study_indices = which(events2$message=='StudyWord') 
  study_indices = c(study_indices, study_indices+1, study_indices+2)
  study_trials = events2[study_indices,] %>%
    arrange(trial, sttime) %>%
    select(trial, sttime, message) %>%
    filter(message!='StudyWord') %>%
    mutate(category = ifelse(str_detect(message, 'Plus|Exes|Check'), 'Cond', 'Word'), phase='study') %>%
    select(trial, phase, category, value=message, sttime)
  
  events_processed = study_trials
  
  # Extract word onset time (only from Word category)
  word_onset_df <- events_processed %>%
    filter(category == "Word") %>%
    select(trial, phase, word_onset = sttime)  # Ensure phase is included for uniqueness
  
  # Remove sttime before pivoting
  events_wide <- events_processed %>%
    select(-sttime) %>%
    pivot_wider(
      names_from = category,
      values_from = value
    ) %>%
    rename(word = Word, condition = Cond)  # Rename columns for clarity
  
  # Clean up extracted values and join with word onset times
  events_wide <- events_wide %>%
    mutate(
      word = str_remove(word, "^Word ") %>% tolower(),  # Remove "Word " prefix & convert to lowercase
      condition = str_remove(condition, "^Cond "),  # Remove "Cond " prefix
      condition = recode(condition, "Exes.png" = "silent", "Plus.png" = "aloud", 'Check.png' = 'control', "foil" = "foil")  # Recode condition names
    ) %>%
    left_join(word_onset_df, by = c("trial", "phase")) %>%  # Add word onset time per trial and phase
    select(trial, phase, word, condition, word_onset)  # Keep only necessary columns
  
  return(events_wide)
}

extract_trial_metadata_e3 <- function(events) {
  
  # Get study phase indices
  events2 = events %>% filter(entime < 10, message != "", !str_detect(message, "TRACKER_TIME|TRIAL_VAR"))

  study_indices = which(str_detect(events2$message,'Plus|Exes|Check'))
  study_indices = c(study_indices, study_indices+1, study_indices+2)
  
  study_trials = events2[study_indices,] %>%
    arrange(trial, sttime) %>%
    select(trial, sttime, message) %>%
    filter(message!='StudyWord') %>%
    mutate(category = ifelse(str_detect(message, 'Plus|Exes|Check'), 'Cond', 'Word'), phase='study') %>%
    select(trial, phase, category, value=message, sttime)

  events_processed = study_trials
  
  # Extract cond onset time (only from Cond category): Still called word_onset so everything still works
  word_onset_df <- events_processed %>%
    filter(category == "Cond") %>%
    select(trial, phase, word_onset = sttime)  # Ensure phase is included for uniqueness
  
  # Remove sttime before pivoting
  events_wide <- events_processed %>%
    select(-sttime) %>%
    pivot_wider(
      names_from = category,
      values_from = value
    ) %>%
    rename(word = Word, condition = Cond)  # Rename columns for clarity
  
  # Clean up extracted values and join with word onset times
  events_wide <- events_wide %>%
    mutate(
      word = str_remove(word, "^Word ") %>% tolower(),  # Remove "Word " prefix & convert to lowercase
      condition = str_remove(condition, "^Cond "),  # Remove "Cond " prefix
      condition = recode(condition, "Exes.png" = "silent", "Plus.png" = "aloud", 'Check.png' = 'control', "foil" = "foil")  # Recode condition names
    ) %>%
    left_join(word_onset_df, by = c("trial", "phase")) %>%  # Add word onset time per trial and phase
    select(trial, phase, word, condition, word_onset)  # Keep only necessary columns
  
  return(events_wide)
}

extract_trial_metadata_e4 <- function(events) {
  
  # Get study phase indices
  events2 = events %>% filter(entime < 10, message != "", !str_detect(message, "TRACKER_TIME|TRIAL_VAR"))
  
  study_indices = which(str_detect(events2$message,'StudyWord'))
  study_indices = c(study_indices, study_indices+1, study_indices+2)
  
  study_trials = events2[study_indices,] %>%
    arrange(trial, sttime) %>%
    select(trial, sttime, message) %>%
    filter(message!='StudyWord') %>%
    mutate(category = ifelse(str_detect(message, 'Plus|Exes|Check'), 'Cond', 'Word'), phase='study') %>%
    select(trial, phase, category, value=message, sttime)
  
  events_processed = study_trials
  
  # Extract study word again
  word_onset_df <- events_processed %>%
    filter(category == "Word") %>%
    select(trial, phase, word_onset = sttime)  # Ensure phase is included for uniqueness
  
  # Remove sttime before pivoting
  events_wide <- events_processed %>%
    select(-sttime) %>%
    pivot_wider(
      names_from = category,
      values_from = value
    ) %>%
    rename(word = Word, condition = Cond)  # Rename columns for clarity
  
  # Clean up extracted values and join with word onset times
  events_wide <- events_wide %>%
    mutate(
      word = str_remove(word, "^Word ") %>% tolower(),  # Remove "Word " prefix & convert to lowercase
      condition = str_remove(condition, "^Cond "),  # Remove "Cond " prefix
      condition = recode(condition, "Exes.png" = "silent", "Plus.png" = "aloud", 'Check.png' = 'control', "foil" = "foil")  # Recode condition names
    ) %>%
    left_join(word_onset_df, by = c("trial", "phase")) %>%  # Add word onset time per trial and phase
    select(trial, phase, word, condition, word_onset)  # Keep only necessary columns
  
  return(events_wide)
}

# Blink Detection via Eyelink 1000 Plus -----------------------------------

mark_blinks <- function(pupil_data, blink_data, pre_blink = 100, post_blink = 100) {
  #message("    Marking blinks...")
  
  pupil_data <- pupil_data %>%
    mutate(blink = FALSE)  # Initialize blink column
  
  # Initialize progress bar
  pb <- progress_bar$new(
    format = "Processing blinks [:bar] :percent (:current/:total) ETA: :eta",
    total = nrow(blink_data),
    width = 60
  )
  
  for (i in seq_len(nrow(blink_data))) {
    pb$tick()  # Update progress bar
    
    # Expand the blink window
    blink_start <- blink_data$sttime[i] - pre_blink
    blink_end <- blink_data$entime[i] + post_blink
    
    pupil_data <- pupil_data %>%
      mutate(
        blink = ifelse(
          trial == blink_data$trial[i] & 
            time >= blink_start & 
            time <= blink_end, 
          TRUE, blink
        )
      )
  }
  
  return(pupil_data)
}

# Drop-out Detection ------------------------------------------------------

detect_dropouts <- function(data, sampling_rate = 1000, threshold = 0.05, window_ms = 500) {  
  #message("    Detecting drop-outs...")
  
  binwidth = max(round((window_ms / 1000) * sampling_rate), 2)
  
  data <- data %>%
    group_by(trial) %>%
    mutate(
      dropout = case_when(
        # If all pupil size values are NA for more than 500ms (50 samples @1000Hz)
        zoo::rollapply(is.na(rps), width = binwidth, FUN = all, fill = FALSE, align = "left") ~ TRUE,
        # If pupil size does not change significantly over 500ms, mark as dropout
        zoo::rollapply(rps, width = binwidth, FUN = function(x) sd(x, na.rm = TRUE) < threshold, fill = FALSE, align = "left") ~ TRUE,
        is.na(x) | is.na(y) ~ TRUE, # Missing gaze data
        rps == 0 ~ TRUE, # Signal is 0
        TRUE ~ FALSE
      )
    ) %>%
    ungroup()
  
  return(data)
}

# Extreme Value Detection -------------------------------------------------

detect_artifacts <- function(data, sampling_rate = 1000, cutoff_factor = 5) {
  # 1) Smooth pupil unless blink/dropout
  # 2) Compute derivative as "difference over actual time delta"
  # 3) Use a global MAD-based threshold
  
  data <- data %>%
    group_by(trial) %>%
    mutate(
      # Only smooth “good” pupil data
      pupil_smoothed = ifelse(blink | dropout, NA, rps),
      
      # Simple rolling mean for smoothing, or skip if < 5 samples
      pupil_smoothed = if (n() < 5) rps else rollapply(
        data    = pupil_smoothed, 
        width   = min(5, n()),
        FUN     = mean,
        fill    = NA,
        align   = "center",
        partial = TRUE
      ),
      
      # Interpolate small gaps
      pupil_smoothed = zoo::na.approx(pupil_smoothed, na.rm = FALSE, maxgap = 100),
      
      # Derivative in pupil-units per ms, using forward difference
      #    (pupil[t] - pupil[t-1]) / (time[t] - time[t-1])
      pupil_diff = (pupil_smoothed - lag(pupil_smoothed)) / (time - lag(time))
    ) %>%
    ungroup()
  
  # Compute a global MAD on the derivative, ignoring NA
  valid_diff <- na.omit(data$pupil_diff)
  global_mad <- mad(valid_diff, na.rm = TRUE)
  
  # As a safeguard, also ensure the MAD is never below some small quantile
  global_mad <- max(global_mad, quantile(valid_diff, 0.01, na.rm = TRUE))
  
  data <- data %>%
    mutate(
      # Flag artifact if derivative is > cutoff * global MAD
      artifact = abs(pupil_diff) > (cutoff_factor * global_mad),
      # If derivative is NA (first sample per trial, etc.), do not flag
      artifact = ifelse(is.na(pupil_diff), FALSE, artifact)
    )
  
  return(data)
}

# Mark Extreme Gaze Data --------------------------------------------------

missing_gaze <- function(data, 
                             lower_x = .2*1024, upper_x = .8*1024, 
                             lower_y = .2*768, upper_y = .8*768) {

  data <- data %>%
    mutate(
      missing_gaze = (is.na(x) | (x < lower_x  | x > upper_x)) | 
             (is.na(y) | (y < lower_y | y > upper_y))
      )
  
  return(data)
}

# Mark Unreliable Data ----------------------------------------------------

mark_unreliable_data <- function(data) {
  #message("    Marking unreliable data...")
  
  data <- data %>%
    mutate(
      ps = ifelse(blink | dropout | artifact | missing_gaze, NA, rps),
    )
  return(data)
}

# Check for excessive missing data per trial ------------------------------

remove_excessive_missing <- function(data, exclusion_threshold = 0.30) {
  # Compute % of missing data per trial
  missing_summary <- data %>%
    group_by(trial) %>%
    summarise(
      total_samples = n(),
      missing_samples = sum(is.na(ps)),
      missing_percent = missing_samples / total_samples
    )
  
  # Identify trials that exceed the threshold
  excluded_trials <- missing_summary %>%
    filter(missing_percent > exclusion_threshold) %>%
    pull(trial)
  
  data = data %>% mutate(bad_trial_missing = trial %in% excluded_trials)
  
  message(sprintf("    Flagged %d trials with >%.0f%% missing data", length(excluded_trials), exclusion_threshold * 100))
  
  return(data)
}

# Interpolation -----------------------------------------------------------

interpolate_pupil <- function(data) {
  
  # Identify valid (non-NA) pupil size indices
  valid_idx <- !is.na(data$ps)
  
  if(sum(valid_idx) < 2) {
    # Not enough valid data to interpolate, return as-is
    data <- data %>%
      mutate(
        interpolated = is.na(ps),
        ps_pc = ps,
        ps_cs = ps
      )
    return(data)
  }
  
  # Get valid time and pupil values
  x_valid <- data$time[valid_idx]
  y_valid <- data$ps[valid_idx]
  
  # Clamp time values to the valid range
  min_x <- min(x_valid)
  max_x <- max(x_valid)
  xi <- pmin(pmax(data$time, min_x), max_x)
  
  # Interpolation using pchip
  interpolated_values <- signal::interp1(
    x = x_valid,
    y = y_valid,
    xi = xi,
    method = "pchip"
  )
  
  # Apply interpolated values only to missing data
  data <- data %>%
    mutate(
      interpolated = is.na(ps),
      ps_pc = ifelse(interpolated, interpolated_values, ps),
      ps_cs = if(sum(valid_idx) >= 2) zoo::na.spline(ps, x = time, na.rm = FALSE) else ps
    )
  
  return(data)
}

# Low Pass Filter ---------------------------------------------------------

lowpass_filter <- function(data, sampling_rate = 1000, cutoff = 4) {

  # 4th-order Butterworth low-pass filter
  nyquist <- sampling_rate / 2
  bf <- butter(n = 4, W = cutoff / nyquist, type = "low")  # Create filter
  
  data %>%
    mutate(
      ps_pc = zoo::na.locf(ps_pc, na.rm = FALSE, fromLast = TRUE),  # Fill trailing NAs
      ps_pc = zoo::na.locf(ps_pc, na.rm = FALSE, fromLast = FALSE), # Fill leading NAs
      psf = filtfilt(bf, ps_pc)  # Apply zero-phase filtering
    ) %>%
    ungroup() -> data

  return(data)
}

# Epoch the Data ----------------------------------------------------------

epoch_data <- function(data, lower_end = -200, upper_end = 7000) {
  data = data %>%
    mutate(epoch_time = time - word_onset) %>%  # Calculate epoch time relative to word onset
    filter(epoch_time >= lower_end, epoch_time <= upper_end)  # Keep only relevant time range
  
  return(data)
}

# Baseline Correct Data ---------------------------------------------------

baseline_pupil <- function(data, baseline_window = c(-200, 0)) {
 
  # Perform baselining for each trial
  data <- data %>%
    group_by(trial) %>%
    mutate(
      # Identify baseline period
      is_baseline = epoch_time >= baseline_window[1] & epoch_time <= baseline_window[2],
      
      # Compute baseline mean for each trial
      baseline_mean = mean(psf[is_baseline], na.rm = TRUE),
      
      # Subtract baseline mean to normalize pupil size
      ps_bs = psf - baseline_mean, ps_bd = 100 * psf / baseline_mean - 100
    ) %>%
    ungroup() %>%
    select(-is_baseline)  # Remove helper column
  
  return(data)
}

# Downsample --------------------------------------------------------------

downsample_pupil <- function(data, target_rate = 50) {
  # Compute the downsampling interval (time step in ms)
  downsample_interval <- 1000 / target_rate  # e.g., 1000ms / 50Hz = 20ms
  
  # Create bin edges based on the epoch time range
  bin_edges <- seq(min(data$epoch_time, na.rm = TRUE), max(data$epoch_time, na.rm = TRUE), by = downsample_interval)
  
  # Compute midpoints of each bin
  bin_midpoints <- bin_edges[-length(bin_edges)] + (downsample_interval / 2)
  
  # Assign each epoch_time a bin using midpoints
  data <- data %>%
    mutate(bin = cut(epoch_time, breaks = bin_edges, labels = bin_midpoints, include.lowest = TRUE)) %>%
    mutate(bin = as.numeric(as.character(bin)))  # Convert factor to numeric
  
  # Downsample each trial using the same time grid
  data %>%
    group_by(sid, trial, bin, phase, word, condition) %>%
    reframe(
      x = mean(x),
      y = mean(y),
      pupil = mean(ps_bs),  # Interpolate to align with grid
      pmissing = mean(interpolated)
    ) -> data
  
  return(data)
}

# Flag Bad Trials ---------------------------------------------------------

bad_trials_by_missing_data <- function(data, missing_threshold_per_bin = 0.5, bad_bin_threshold_per_trial = 0.3) {
  # Calculate the proportion of bins with missing data above the threshold for each trial
  bad_trials <- data %>%
    group_by(sid, trial) %>%
    summarise(
      bins_with_high_missing = mean(pmissing > missing_threshold_per_bin, na.rm = TRUE),
      total_bins = n()
    ) %>%
    ungroup() %>%
    mutate(bad_trial = bins_with_high_missing > bad_bin_threshold_per_trial)
  
  # Join bad trial information back to the original data
  data <- data %>%
    left_join(bad_trials %>% select(sid, trial, bad_trial), by = c("sid", "trial"))
  
  return(data)
}

# Normalize the Data ------------------------------------------------------

normalize_pupil_zscore <- function(data) {
  data %>%
    mutate(
      pupil_z = (pupil - mean(pupil, na.rm = TRUE)) / sd(pupil, na.rm = TRUE)
    ) -> data
  
  return(data)
}

#
#
# Analysis Functions ------------------------------------------------------
#
#


# GAMM Functions ----------------------------------------------------------

fit_gam_study_w_word <- function(data, filename) {
  if(file.exists(filename)) {
    message('Loading file...')
    gam_model = readRDS(filename)
  } else {
    
    nc <- parallelly::availableCores() - 1
    
    if (nc > 0) { 
      cl <- parallel::makeCluster(nc) 
    } else {
      cl <- NULL
    }
    
    data %>%
      mutate_at(vars(condition, sid, word), factor) -> data
    
    message('Fitting model...')
    gam_model <- bam(pupil_z ~ condition + 
                       s(bin, by = condition, k = 12) +  # Smooth effect of time (bin), modeled separately for each condition
                       s(bin, sid, by = condition, bs = 'fs', m = 1) +   # Random smooth per subject (to allow individual trajectories)
                       s(bin, word, by = condition, bs = 'fs', m = 1),  # Random smooth per trial (to account for trial-level variability)
                     data = data,                 
                     method = "fREML",
                     family = gaussian(),
                     cluster=cl,
                     AR.start = data$bin == min(data$bin), rho = 0.9)  # Handles autocorrelation
    
    
    if (!is.null(cl)) {
      parallel::stopCluster(cl)
    }
    
    message('Saving file...')
    saveRDS(gam_model, filename)
  }
  
  return(gam_model)
}

fit_gam_study_mem <- function(data, filename) {
  if(file.exists(filename)) {
    message('Loading file...')
    gam_model = readRDS(filename)
  } else {
    
    nc <- parallelly::availableCores() - 1
    
    if (nc > 0) { 
      cl <- parallel::makeCluster(nc) 
    } else {
      cl <- NULL
    }
    
    data %>%
      mutate_at(vars(condmem, sid, word, condition), factor) -> data
    
    message('Fitting model...')
    gam_model <- bam(pupil_z ~ condmem + 
                       s(bin, by = condmem, k = 12) +  # Smooth effect of time (bin), modeled separately for each condition
                       s(bin, sid, by = condmem, bs = 'fs', m = 1) +   # Random smooth per subject (to allow individual trajectories)
                       s(bin, word, by = condition, bs = 'fs', m = 1),  # Random smooth per trial (to account for trial-level variability)
                     data = data,                 
                     method = "fREML",
                     family = gaussian(),
                     cluster=cl,
                     AR.start = data$bin == min(data$bin), rho = 0.9)  # Handles autocorrelation
    
    
    if (!is.null(cl)) {
      parallel::stopCluster(cl)
    }
    
    message('Saving file...')
    saveRDS(gam_model, filename)
  }
  
  return(gam_model)
}

# Mass Univariate ---------------------------------------------------------

fit_massu_t = function(data, comparison, conditions)
{
  data = data %>% 
    filter(condition %in% conditions) %>%
    mutate(condition=as.factor(condition)) %>%
    mutate(condition=droplevels(condition))
  
  results <- data.frame(comparison = character(), bin = numeric(), t_value = numeric(), p_value = numeric())
  
  bins <- unique(data$bin)
  
  for (b in bins) {
    temp_data <- data %>% filter(bin == b)
    
    # Ensure there are exactly two conditions
    conditions <- unique(temp_data$condition)
    if (length(conditions) != 2) {
      next  # Skip if there aren't exactly two conditions
    }
    
    cond1_data <- temp_data %>% filter(condition == conditions[1]) %>% pull(pupil_z)
    cond2_data <- temp_data %>% filter(condition == conditions[2]) %>% pull(pupil_z)
    
    if (length(cond1_data) == length(cond2_data)) {  # Ensure paired sample
      t_test_result <- t.test(cond1_data, cond2_data, paired = TRUE)
      
      results <- rbind(results, data.frame(comparison=comparison, bin = b, t_value = t_test_result$statistic, p_value = t_test_result$p.value))
    }
  }
  
  results = results %>% 
    mutate(adj_p_value = p.adjust(p_value, method = "fdr"))
  
  return(results)
}

# Plotting functions ------------------------------------------------------

create_plot_sig <- function(gamm_model, massu, exp='EXP') {
  # Get significant bins from GAMM comparisons
  gamm_plot_sig <- plot_diff(
    gamm_model,
    view  = 'bin',
    comp  = list(condition = c('silent', 'control')),
    n.grid = 1000,
    plot   = FALSE
  ) %>%
    dplyr::select(est, bin, CI) %>%
    dplyr::mutate(
      comparison = 'Control - Silent',
      sig        = est + CI < 0 | est - CI > 0
    ) %>%
    dplyr::bind_rows(
      plot_diff(
        gamm_model,
        view  = 'bin',
        comp  = list(condition = c('aloud', 'control')),
        n.grid = 1000,
        plot   = FALSE
      ) %>%
        dplyr::select(est, bin, CI) %>%
        dplyr::mutate(
          comparison = 'Aloud - Control',
          sig        = est + CI < 0 | est - CI > 0
        )
    ) %>%
    dplyr::bind_rows(
      plot_diff(
        gamm_model,
        view  = 'bin',
        comp  = list(condition = c('aloud', 'silent')),
        n.grid = 1000,
        plot   = FALSE
      ) %>%
        dplyr::select(est, bin, CI) %>%
        dplyr::mutate(
          comparison = 'Aloud - Silent',
          sig        = est + CI < 0 | est - CI > 0
        )
    ) %>%
    dplyr::filter(sig) %>%
    dplyr::mutate(
      y_offset = dplyr::case_when(
        comparison == "Aloud - Silent"   ~ -0.6,
        comparison == "Aloud - Control"  ~ -0.65,
        comparison == "Control - Silent" ~ -0.7
      ),
      exp   = exp,
      model = 'GAMM'
    ) %>%
    dplyr::select(bin, comparison, y_offset, exp, model)
  
  # Filter the massu data for significant bins and tag them as 'Mass Univariate'
  massu_sig <- massu %>%
    dplyr::filter(adj_p_value < 0.05) %>%
    dplyr::mutate(
      y_offset = dplyr::case_when(
        comparison == "Aloud - Silent"   ~ -0.6,
        comparison == "Aloud - Control"  ~ -0.65,
        comparison == "Control - Silent" ~ -0.7
      ),
      exp   = exp,
      model = 'Mass Univariate'
    ) %>%
    tibble::rownames_to_column("original_names") %>%
    dplyr::mutate(original_names = dplyr::row_number()) %>%
    tibble::column_to_rownames("original_names") %>%
    dplyr::select(bin, comparison, y_offset, exp, model)
  
  # Combine GAMM significance with Mass Univariate significance
  plot_sig <- massu_sig %>%
    dplyr::bind_rows(gamm_plot_sig) %>%
    dplyr::mutate(
      model = factor(model, levels = c('Mass Univariate', 'GAMM'))
    )
  
  return(plot_sig)
}

create_plot_sig_cluster <- function(gamm_model, massu, exp='EXP') {
  # Get significant bins from GAMM comparisons
  gamm_plot_sig <- plot_diff(
    gamm_model,
    view  = 'bin',
    comp  = list(condition = c('silent', 'control')),
    n.grid = 1000,
    plot   = FALSE
  ) %>%
    dplyr::select(est, bin, CI) %>%
    dplyr::mutate(
      comparison = 'Control - Silent',
      sig        = est + CI < 0 | est - CI > 0
    ) %>%
    dplyr::bind_rows(
      plot_diff(
        gamm_model,
        view  = 'bin',
        comp  = list(condition = c('aloud', 'control')),
        n.grid = 1000,
        plot   = FALSE
      ) %>%
        dplyr::select(est, bin, CI) %>%
        dplyr::mutate(
          comparison = 'Aloud - Control',
          sig        = est + CI < 0 | est - CI > 0
        )
    ) %>%
    dplyr::bind_rows(
      plot_diff(
        gamm_model,
        view  = 'bin',
        comp  = list(condition = c('aloud', 'silent')),
        n.grid = 1000,
        plot   = FALSE
      ) %>%
        dplyr::select(est, bin, CI) %>%
        dplyr::mutate(
          comparison = 'Aloud - Silent',
          sig        = est + CI < 0 | est - CI > 0
        )
    ) %>%
    dplyr::filter(sig) %>%
    dplyr::mutate(
      y_offset = dplyr::case_when(
        comparison == "Aloud - Silent"   ~ -0.6,
        comparison == "Aloud - Control"  ~ -0.65,
        comparison == "Control - Silent" ~ -0.7
      ),
      exp   = exp,
      model = 'GAMM'
    ) %>%
    dplyr::select(bin, comparison, y_offset, exp, model)
  
  # Filter the massu data for significant bins and tag them as 'Mass Univariate'
  massu_sig <- massu %>%
    dplyr::filter(sig_clust) %>%
    dplyr::mutate(
      y_offset = dplyr::case_when(
        comparison == "Aloud - Silent"   ~ -0.6,
        comparison == "Aloud - Control"  ~ -0.65,
        comparison == "Control - Silent" ~ -0.7
      ),
      exp   = exp,
      model = 'Mass Univariate'
    ) %>%
    tibble::rownames_to_column("original_names") %>%
    dplyr::mutate(original_names = dplyr::row_number()) %>%
    tibble::column_to_rownames("original_names") %>%
    dplyr::select(bin, comparison, y_offset, exp, model)
  
  # Combine GAMM significance with Mass Univariate significance
  plot_sig <- massu_sig %>%
    dplyr::bind_rows(gamm_plot_sig) %>%
    dplyr::mutate(
      model = factor(model, levels = c('Mass Univariate', 'GAMM'))
    )
  
  return(plot_sig)
}

create_study_pd_sum_figure <- function(pupil_dat, gamm_model, exp='EXP') {
  # From the GAMM model, derive smoothed predictions
  study_pgam_plot <- plot_smooth(
    gamm_model, 
    view      = "bin",
    main      = "Estimated Pupil Response Over Time",
    rug       = TRUE,
    n.grid    = 1000,
    plot_all  = "condition"
  )$fv %>%
    dplyr::select(condition, bin, fit, CI, ll, ul) %>%
    dplyr::mutate(exp = exp, model = 'GAMM')
  
  # Summarize the raw pupil data
  massu_summary <- pupil_dat %>%
    dplyr::group_by(bin, condition) %>%
    dplyr::summarize(
      m  = mean(pupil_z),
      se = sd(pupil_z) / sqrt(dplyr::n()),
      ci = se * 1.96,
      .groups = 'drop'
    ) %>%
    dplyr::mutate(model = 'Mass Univariate')
  
  # Combine the raw data summary with the GAMM predictions
  study_pd_sum_figure <- massu_summary %>%
    dplyr::bind_rows(
      study_pgam_plot %>%
        dplyr::mutate(se = CI / 1.96) %>%
        dplyr::select(
          bin,
          condition,
          m   = fit,
          se,
          ci  = CI,
          model
        )
    ) %>%
    dplyr::mutate(
      model = factor(model, levels = c('Mass Univariate', 'GAMM')),
      exp   = exp
    )
  
  return(study_pd_sum_figure)
}

create_pupil_figure <- function(gamm_model, massu, pupil_dat, exp='EXP') {
  # Use the new helper function to create the significance data (plot_sig)
  plot_sig <- create_plot_sig(
    gamm_model = gamm_model,
    massu      = massu,
    exp        = exp
  )
  
  # Use the new helper function to create the summarized data (study_pd_sum_figure)
  study_pd_sum_figure <- create_study_pd_sum_figure(
    pupil_dat  = pupil_dat,
    gamm_model = gamm_model,
    exp        = exp
  )
  
  # Define labels and offsets for the segment annotations
  label_positions <- data.frame(
    comparison     = c("A-S", "A-C", "C-S"),
    label          = c("A",    "B",   "C"),
    y_label_offset = c(-0.61, -0.66, -0.71)  # Slightly above the segments
  )
  
  # Construct the final ggplot using the newly created data frames
  fig <- study_pd_sum_figure %>%
    ggplot(aes(
      x = bin, 
      y = m,
      fill = condition, 
      color = condition, 
      linetype = condition,
      ymin = m - ci, 
      ymax = m + ci
    )) +
    geom_line() +
    geom_ribbon(alpha = 0.5) +
    geom_segment(
      data = plot_sig,
      aes(x = bin, xend = bin, y = y_offset, yend = y_offset - 0.02),
      inherit.aes = FALSE, 
      size = 1
    ) +
    geom_text(
      data = label_positions,
      aes(
        x = min(massu$bin) - 10,  # Align all labels at a consistent x position
        y = y_label_offset, 
        label = comparison
      ),
      inherit.aes = FALSE, 
      hjust = 0, 
      size = 2.5
    ) +
    labs(x = "Time (ms)", y = "z(Pupil)") +
    theme_classic() +
    geom_hline(yintercept = 0, colour = 'grey', alpha = 0.2) +
    facet_grid(. ~ model) +
    theme(
      # Make facet strips black with white text
      strip.background = element_rect(fill = "black", color = NA),
      strip.text.x     = element_text(color = "white"),
      
      # Transparent panel and plot background
      panel.background = element_rect(fill = NA, color = NA),
      plot.background  = element_rect(fill = NA, color = NA),
      
      # Make axis & legend text visible
      axis.text  = element_text(color = "black"),
      axis.title = element_text(color = "black")
    )
  
  return(fig)
}

bayesian_correlation_test <- function(data, comparisons, filename) {
  if(file.exists(filename)) {
    message("Loading file...")
    results_df = readRDS(filename)
  } else {
    
    results_df = data.frame()
    for(p in comparisons)
    {
      iv1 = p[[1]]
      iv2 = p[[2]]
      
      print(sprintf('%s %s', iv1, iv2))
      
      cor_result <- cor.test(data[[iv1]], data[[iv2]], method='mcmc')
      
      # Extract relevant summary statistics
      estimate <- round(cor_result$estimate, 3)
      ci_lower <- round(cor_result@hdi.l, 3)
      ci_upper <- round(cor_result@hdi.u, 3)
      
      # Mark significant correlations where CI does not include 0
      sig_marker <- ifelse(ci_lower > 0 | ci_upper < 0, "*", "")
      
      results_df = bind_rows(results_df, data.frame(IV1 = iv1, IV2 = iv2, Estimate = estimate, 
                 CI_Lower = ci_lower, CI_Upper = ci_upper, Significant = sig_marker))
    }

    saveRDS(results_df, filename)
  }
  
  return(results_df)
}

prepare_behav_dat_figure_e1 <- function(d_means, d_comps, exp) {
  # Adjust d' values for new-X contrasts (make positive instead of negative)
  d_comps <- d_comps %>%
    data.frame() %>%
    mutate(estimate = ifelse(grepl("new -", contrast), -estimate, estimate),
           lower.HPD = ifelse(grepl("new -", contrast), -lower.HPD, lower.HPD),
           upper.HPD = ifelse(grepl("new -", contrast), -upper.HPD, upper.HPD))
  
  # Create individual data frames for each panel
  df1 <- d_comps %>% filter(grepl("new -", contrast)) %>% mutate(contrast = str_remove(contrast, 'new - '), type = "d'")
  df2 <- d_comps %>% filter(grepl("aloud -|control - silent", contrast)) %>% mutate(type = "d' Contrasts")
  
  # Combine all data frames
  combined_df <- bind_rows(
    df1 %>% rename(value = estimate, condition = contrast),
    df2 %>% rename(value = estimate, condition = contrast),
    
  ) %>% mutate(
    type = factor(type, levels = c("d'", "d' Contrasts", "Rec.", "Rec. Contrasts", "Fam.", "Fam. Contrasts")),
    condition = factor(condition, levels = c("aloud", "control", "silent", "aloud - silent", "aloud - control", "control - silent")),
    experiment = exp
  )
  
  return(combined_df)
}

prepare_behav_dat_figure <- function(d_means, d_comps, rho_means, rho_comps, f_means, f_comps, exp) {
  # Adjust d' values for new-X contrasts (make positive instead of negative)
  d_comps <- d_comps %>%
    data.frame() %>%
    mutate(estimate = ifelse(grepl("new -", contrast), -estimate, estimate),
           lower.HPD = ifelse(grepl("new -", contrast), -lower.HPD, lower.HPD),
           upper.HPD = ifelse(grepl("new -", contrast), -upper.HPD, upper.HPD))
  
  # Create individual data frames for each panel
  df1 <- d_comps %>% filter(grepl("new -", contrast)) %>% mutate(contrast = str_remove(contrast, 'new - '), type = "d'")
  df2 <- d_comps %>% filter(grepl("aloud -|control - silent", contrast)) %>% mutate(type = "d' Contrasts")
  df3 <- rho_means %>% data.frame() %>% mutate(type = "Rec.")
  df4 <- rho_comps %>% data.frame() %>% mutate(type = "Rec. Contrasts")
  df5 <- f_means %>% data.frame() %>% mutate(type = "Fam.")
  df6 <- f_comps %>% data.frame() %>% mutate(type = "Fam. Contrasts")
  
  # Combine all data frames
  combined_df <- bind_rows(
    df1 %>% rename(value = estimate, condition = contrast),
    df2 %>% rename(value = estimate, condition = contrast),
    df3 %>% rename(value = emmean),
    df4 %>% rename(value = estimate, condition = contrast),
    df5 %>% rename(value = emmean),
    df6 %>% rename(value = estimate, condition = contrast)
  ) %>% mutate(
    type = factor(type, levels = c("d'", "d' Contrasts", "Rec.", "Rec. Contrasts", "Fam.", "Fam. Contrasts")),
    condition = factor(condition, levels = c("aloud", "control", "silent", "aloud - silent", "aloud - control", "control - silent")),
    experiment = exp
  )
  
  return(combined_df)
}

behav_dat_plot <- function(combined_df) {
  p <- combined_df %>% ggplot(aes(x = condition, y = value, shape = condition_group, color = condition_group)) +
    geom_point() +
    geom_hline(yintercept = 0, colour = 'grey', alpha = 0.2) +
    geom_errorbar(aes(ymin = lower.HPD, ymax = upper.HPD), width = 0.1) +
    facet_grid(experiment~type, scales='free_x') + #, scales = "free", space = "free") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(y = "Estimate", x = "") +
    theme(
      # Make facet strips black with white text
      strip.background = element_rect(fill = "black", color = NA),
      strip.text.x     = element_text(color = "white"),
      strip.text.y     = element_text(color = "white"),
      
      # Transparent panel and plot background
      panel.background = element_rect(fill = NA, color = NA),
      plot.background  = element_rect(fill = NA, color = NA),
      
      legend.position = "top",  # Move legend under plot
      legend.title = element_blank(),
      
      # Make axis & legend text visible
      axis.text  = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      panel.border = element_rect(color = "black", fill = NA, size = 1)
    ) + scale_color_manual(values = condition_colors) +  # Apply custom colors
    scale_shape_manual(values = condition_shapes) 
  return(p)
}

# Correcting reversed responding ------------------------------------------

invert_responses = function(x, max_val = 6)
{
  x = (max_val+1) - x
}

# Rolling Window w/ Cluster Correction ------------------------------------

fit_model_single_window_cluster <- function(data, time_point) {
  
  # Subset data within the 100 ms window centered at 'time_point'
  window_data <- data[data$bin >= (time_point - 50) & data$bin <= (time_point + 50), ]
  
  # Compute mean pupil size per trial within this window
  window_aggregated <- aggregate(pupil_z ~ sid + word + acc + condition, data = window_data, FUN = mean)
  
  # Fit the logistic mixed model
  model <- tryCatch(
    glmer(acc ~ condition + pupil_z : condition-1 + (condition-1 | sid) + (condition-1 | word),
          data = window_aggregated,
          family = binomial(link = "logit"),
          control = glmerControl(optimizer = "bobyqa")),  # More stable optimizer
    error = function(e) return(NULL)  # Handle model failures
  )
  
  # If model fails to fit, return NA values
  if (is.null(model)) {
    return(data.frame(
      time_point = time_point,
      silent_p = NA,
      silent_OR = NA,
      silent_z = NA,
      aloud_p = NA,
      aloud_OR = NA,
      aloud_z = NA,
      control_p = NA,
      control_OR = NA,
      control_z = NA
    ))
  }
  
  # Extract fixed effects
  results <- broom.mixed::tidy(model, effects = "fixed")
  
  # Extract relevant terms
  silent_effect <- results[results$term == "conditionsilent:pupil_z", ]
  aloud_effect <- results[results$term == "conditionaloud:pupil_z", ]
  control_effect <- results[results$term == "conditioncontrol:pupil_z", ]
  
  aloud_silent = car::linearHypothesis(model, "conditionaloud:pupil_z = conditionsilent:pupil_z") %>% 
    select(`Pr(>Chisq)`) %>% 
    filter(!is.na(`Pr(>Chisq)`)) %>% 
    pull(`Pr(>Chisq)`)
  aloud_control = ifelse(nrow(control_effect)>0, 
                         car::linearHypothesis(model, "conditionaloud:pupil_z = conditioncontrol:pupil_z") %>% 
                           select(`Pr(>Chisq)`) %>% 
                           filter(!is.na(`Pr(>Chisq)`)) %>% 
                           pull(`Pr(>Chisq)`), NA
  ) 
  silent_control = ifelse(nrow(control_effect)>0, 
                          car::linearHypothesis(model, "conditionsilent:pupil_z = conditioncontrol:pupil_z") %>% 
                            select(`Pr(>Chisq)`) %>% 
                            filter(!is.na(`Pr(>Chisq)`)) %>% 
                            pull(`Pr(>Chisq)`), NA
  ) 
  
  # Construct output
  output <- data.frame(
    time_point = time_point,
    silent_p = ifelse(nrow(silent_effect) > 0, silent_effect$p.value, NA),
    silent_OR = ifelse(nrow(silent_effect) > 0, exp(silent_effect$estimate), NA),  # Convert log-odds to OR
    silent_z = ifelse(nrow(silent_effect) > 0, silent_effect$statistic, NA),
    aloud_p = ifelse(nrow(aloud_effect) > 0, aloud_effect$p.value, NA),
    aloud_OR = ifelse(nrow(aloud_effect) > 0, exp(aloud_effect$estimate), NA),
    aloud_z = ifelse(nrow(aloud_effect) > 0, aloud_effect$statistic, NA),
    control_p = ifelse(nrow(control_effect) > 0, control_effect$p.value, NA),
    control_OR = ifelse(nrow(control_effect) > 0, exp(control_effect$estimate), NA),
    control_z = ifelse(nrow(control_effect) > 0, control_effect$statistic, NA),
    aloud_silent = aloud_silent,
    aloud_control = aloud_control,
    silent_control = silent_control
  )
  
  return(output)
}

fit_model_across_timepoints_cluster <- function(full_data, time_points, filename) {
  if(file.exists(filename))
  {
    message('Loading File...')
    results_df = readRDS(filename)
  } else {
    
    full_data = full_data %>%
      mutate(acc = ifelse(mem=='yes', 1, 0), condition = factor(condition, levels=c('silent', 'aloud', 'control')))
    
    # Initialize an empty results list
    results_list <- list()
    
    # Enable progress bar tracking
    handlers("txtprogressbar")  # Ensures progress updates in RStudio
    with_progress({  # Ensures progress updates are forced
      p <- progressor(along = time_points)  # Create progress tracker
      
      # Apply function with progress tracking
      results_list <- lapply(time_points, function(tp) {
        p(sprintf("Processing time point: %d ms", tp))  # Update progress bar
        fit_model_single_window_cluster(full_data, tp)
      })
    })
    # Combine results into a single dataframe
    results_df <- do.call(rbind, results_list)
    saveRDS(results_df, filename)
  }
  
  return(results_df)
}

find_clusters <- function(time_points, stat, z_thresh = 1.96, two_sided = TRUE) {
  
  if (two_sided) {
    above_thresh <- abs(stat) >= z_thresh
  } else {
    above_thresh <- stat >= z_thresh
  }
  
  clusters <- list()
  in_clust <- FALSE
  start_i  <- NA
  
  for (i in seq_along(stat)) {
    if (above_thresh[i] && !in_clust) {
      # Start cluster
      in_clust <- TRUE
      start_i  <- i
    }
    if ((!above_thresh[i] || i == length(stat)) && in_clust) {
      # end cluster
      end_i <- if (above_thresh[i]) i else (i - 1)
      cluster_mass <- sum(abs(stat[start_i:end_i]))
      clusters[[length(clusters)+1]] <- data.frame(
        start_idx = start_i,
        end_idx   = end_i,
        start_tp  = time_points[start_i],
        end_tp    = time_points[end_i],
        size      = (end_i - start_i + 1),
        mass      = cluster_mass
      )
      in_clust <- FALSE
    }
  }
  
  if (length(clusters) == 0) {
    return(NULL)
  } else {
    do.call(rbind, clusters)
  }
}

ar1_cluster_correction <- function(results_df,
                                   effect_z_col    = "silent_z",
                                   effect_p_col    = "silent_p",
                                   alpha           = 0.05,
                                   n_sim           = 2000,
                                   two_sided       = TRUE,
                                   z_thresh        = 1.96) {
  
  # Extract observed data
  z_obs <- results_df[[effect_z_col]]
  p_obs <- results_df[[effect_p_col]]
  tp    <- results_df$time_point
  T_len = length(z_obs)
  
  # Identify observed clusters
  observed_clusters <- find_clusters(time_points = tp, stat = z_obs,
                                     z_thresh = z_thresh, two_sided = two_sided)
  
  # If no cluster crosses threshold, there's nothing to correct
  if (is.null(observed_clusters)) {
    message("No supra-threshold cluster found at z-threshold = ", z_thresh)
    return(list(observed_clusters=NULL, cutoff=NULL, pvals=NULL))
  }
  
  max_cluster_masses <- numeric(n_sim)
  fit_ar <- ar.yw(z_obs, order.max = 20, aic = TRUE)
  phi_hat <- fit_ar$ar
  sigma2_hat <- var(fit_ar$resid, na.rm = TRUE)
  for (i in seq_len(n_sim)) {
    print(i)
    
    # Simulate using this estimated model
    wn <- scale(simts::gen_gts(T_len, simts::AR(phi = phi_hat, sigma2 = sigma2_hat)))

    sim_clusters <- find_clusters(time_points = tp, stat = wn,
                                  z_thresh    = z_thresh,
                                  two_sided   = two_sided)
    
    if (!is.null(sim_clusters)) {
      max_cluster_masses[i] <- max(sim_clusters$mass, na.rm=TRUE)
    } else {
      max_cluster_masses[i] <- 0
    }
  }
  
  # Derive cutoff for FWER control
  cutoff <- quantile(max_cluster_masses, probs=1 - alpha)
  
  # For each observed cluster, compute a "p-value" as proportion 
  # of simulated maxes that exceed that cluster's mass
  obs_cluster_pvals <- sapply(observed_clusters$mass, function(m) {
    mean(max_cluster_masses >= m)
  })
  
  # add them to the data frame
  observed_clusters$cluster_p_value <- obs_cluster_pvals
  observed_clusters$significant     <- (observed_clusters$mass > cutoff)
  
  list(
    observed_clusters      = observed_clusters,
    max_cluster_masses     = max_cluster_masses,
    cutoff                 = cutoff,
    phi_hat = phi_hat
  )
}

# Cluster Correction for T ------------------------------------------------

ar1_cluster_correction_t <- function(results_df,
                                   effect_t_col    = "t_value",
                                   effect_p_col    = "p_value",
                                   alpha           = 0.05,
                                   n_sim           = 2000,
                                   two_sided       = TRUE,
                                   df = 30,
                                   t_thresh        = qt(1 - alpha / 2, df = df),
                                   t_sd = (df / (df-2))**.5) {
  
  # Extract observed data
  t_obs <- results_df[[effect_t_col]]
  p_obs <- results_df[[effect_p_col]]
  tp    <- results_df$bin
  T_len = length(t_obs)
  
  # Identify observed clusters
  observed_clusters <- find_clusters(time_points = tp, stat = t_obs,
                                     z_thresh = t_thresh, two_sided = two_sided)
  
  # If no cluster crosses threshold, there's nothing to correct
  if (is.null(observed_clusters)) {
    message("No supra-threshold cluster found at t-threshold = ", z_thresh)
    return(list(observed_clusters=NULL, cutoff=NULL, pvals=NULL))
  }
  
  max_cluster_masses <- numeric(n_sim)
  fit_ar <- ar.yw(t_obs, order.max = 20, aic = TRUE)
  phi_hat <- fit_ar$ar
  sigma2_hat <- var(fit_ar$resid, na.rm = TRUE)
  for (i in seq_len(n_sim)) {
    print(i)
    
    # Simulate using this estimated model
    wn <- t_sd * scale(simts::gen_gts(T_len, simts::AR(phi = phi_hat, sigma2 = sigma2_hat)))
    
    sim_clusters <- find_clusters(time_points = tp, stat = wn,
                                  z_thresh    = t_thresh,
                                  two_sided   = two_sided)
    
    if (!is.null(sim_clusters)) {
      max_cluster_masses[i] <- max(sim_clusters$mass, na.rm=TRUE)
    } else {
      max_cluster_masses[i] <- 0
    }
  }
  
  # Derive cutoff for FWER control
  cutoff <- quantile(max_cluster_masses, probs=1 - alpha)
  
  # For each observed cluster, compute a "p-value" as proportion 
  # of simulated maxes that exceed that cluster's mass
  obs_cluster_pvals <- sapply(observed_clusters$mass, function(m) {
    mean(max_cluster_masses >= m)
  })
  
  # add them to the data frame
  observed_clusters$cluster_p_value <- obs_cluster_pvals
  observed_clusters$significant     <- (observed_clusters$mass > cutoff)
  
  list(
    observed_clusters      = observed_clusters,
    max_cluster_masses     = max_cluster_masses,
    cutoff                 = cutoff,
    phi_hat = phi_hat
  )
}

# Outlier Detection -------------------------------------------------------

detect_outliers_2d <- function(
    data_2d,
    mcd_fraction = 0.5,
    mcd_alpha = 0.01,
    iso_ntrees = 100,
    iso_sample_size = min(256, nrow(data_2d)),
    iso_threshold = 0.6
) {
  # Ensure data_2d is at least 2 columns and numeric
  if (!is.data.frame(data_2d) && !is.matrix(data_2d)) {
    stop("data_2d must be a data frame or matrix.")
  }
  if (ncol(data_2d) != 2) {
    stop("Function currently expects exactly two columns (2D data).")
  }
  
  # Convert to matrix if it's a data frame (for Routliers/isotree usage)
  mat_data <- as.matrix(data_2d)
  
  # Make sure the Routliers package is loaded
  if (!requireNamespace("Routliers", quietly = TRUE)) {
    stop("Package 'Routliers' is required but not installed.")
  }
  
  res_mcd <- Routliers::outliers_mcd(mat_data, mcd_fraction, mcd_alpha)
  
  mcd_outliers <- rep(FALSE, nrow(mat_data))
  # If there are any outliers flagged:
  if (!is.null(res_mcd$outliers_pos) && length(res_mcd$outliers_pos) > 0) {
    mcd_outliers[res_mcd$outliers_pos] <- TRUE
  }
  
  if (!requireNamespace("isotree", quietly = TRUE)) {
    stop("Package 'isotree' is required but not installed.")
  }
  # Fit isolation forest
  iso_model <- isotree::isolation.forest(
    data = mat_data,
    ntrees = iso_ntrees,
    sample_size = iso_sample_size
  )
  iso_scores <- predict(iso_model, mat_data, type = "score")
  
  # Flag outliers based on iso_threshold
  iso_outliers <- iso_scores > iso_threshold
  
  # Convert input to data frame if needed, then add columns:
  df_out <- as.data.frame(data_2d)  # safe conversion
  df_out$MCD_Outlier <- mcd_outliers
  df_out$IF_Outlier <- iso_outliers
  df_out$both = df_out$MCD_Outlier & df_out$IF_Outlier
  
  return(df_out)
}

compute_correlations_with_outlier_removal <- function(
    df,
    var_pairs,              # list of character pairs
    group_var = NULL,       # optional grouping variable
    mcd_fraction = 0.5,
    mcd_alpha = 0.01,
    iso_ntrees = 100,
    iso_threshold = 0.6,
    iso_sample_size = NULL  # if NULL, will be set per group
) {
  results <- data.frame(
    IV1 = character(),
    IV2 = character(),
    Group = character(),
    n_outliers_removed = integer(),
    corr_gaussian = numeric(),
    low_ci_gaussian = numeric(),
    high_ci_gaussian = numeric(),
    corr_student = numeric(),
    low_ci_student = numeric(),
    high_ci_student = numeric(),
    nu_student = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (pair in var_pairs) {
    var1 <- pair[1]
    var2 <- pair[2]
    
    # Skip pair if either variable is all NA
    if (all(is.na(df[[var1]])) || all(is.na(df[[var2]]))) {
      message(paste("Skipping pair:", var1, "vs", var2, "- one or both variables are all NA"))
      next
    }
    
    # Handle grouping
    group_levels <- if (!is.null(group_var)) unique(df[[group_var]]) else NA
    
    for (group in group_levels) {
      # Subset by group if needed
      sub_df <- if (!is.na(group)) {
        df %>% filter(.data[[group_var]] == group)
      } else {
        df
      }
      
      # Skip if both vars are all NA in this group
      if (all(is.na(sub_df[[var1]])) || all(is.na(sub_df[[var2]]))) {
        message(paste("Skipping:", var1, "vs", var2, "for group", group, "- all values are NA"))
        next
      }
      
      # Subset to just the two vars
      sub_data <- sub_df %>% select(all_of(c(var1, var2)))
      
      # Run outlier detection
      sample_size <- if (is.null(iso_sample_size)) min(256, nrow(sub_data)) else iso_sample_size
      
      outlier_result <- detect_outliers_2d(
        data_2d = sub_data,
        mcd_fraction = mcd_fraction,
        mcd_alpha = mcd_alpha,
        iso_ntrees = iso_ntrees,
        iso_threshold = iso_threshold,
        iso_sample_size = sample_size
      )
      
      # Skip if outlier function returned NULL
      if (is.null(outlier_result)) next
      
      keep_rows <- which(!outlier_result$both)
      n_outliers_removed <- sum(outlier_result$both, na.rm = TRUE)
      
      if (length(keep_rows) < 5) {
        message(paste("Skipping:", var1, "vs", var2, "for group", group, "- too few rows after outlier removal"))
        next
      }
      
      # Scale variables
      model_data <- sub_data[keep_rows, , drop = FALSE] %>%
        mutate(
          !!var1 := (.[[var1]] - mean(.[[var1]], na.rm = TRUE)) / sd(.[[var1]], na.rm = TRUE),
          !!var2 := (.[[var2]] - mean(.[[var2]], na.rm = TRUE)) / sd(.[[var2]], na.rm = TRUE)
        )
      
      # Build formulas
      form <- bf(as.formula(paste0(var1, " ~ ", var2, " - 1")))
      
      # Priors
      student_priors <- c(
        set_prior('normal(0, 0.5)', class = 'b', lb = -1, ub = 1),
        set_prior('constant(1)', class = 'sigma'),
        set_prior('gamma(2, 0.1)', class = 'nu')
      )
      gaussian_priors <- c(
        set_prior('normal(0, 0.5)', class = 'b', lb = -1, ub = 1),
        set_prior('constant(1)', class = 'sigma')
      )
      
      # Fit models
      fit_student <- brm(
        formula = form,
        data = model_data,
        family = student(),
        prior = student_priors,
        cores = 4, iter = 4000, chains = 4, backend = 'cmdstanr', silent = TRUE, refresh = 0
      )
      fit_gaussian <- brm(
        formula = form,
        data = model_data,
        family = gaussian(),
        prior = gaussian_priors,
        cores = 4, iter = 4000, chains = 4, backend = 'cmdstanr', silent = TRUE, refresh = 0
      )
      
      # Extract summaries
      sum_student <- posterior_summary(fit_student, pars = c("b_", "nu"))
      sum_gauss <- posterior_summary(fit_gaussian, pars = "b_")
      
      slope_name <- grep("^b_", rownames(sum_student), value = TRUE)[1]
      
      results <- rbind(results, data.frame(
        IV1 = var1,
        IV2 = var2,
        Group = ifelse(is.na(group), "All", as.character(group)),
        n_outliers_removed = n_outliers_removed,
        corr_gaussian = sum_gauss[slope_name, "Estimate"],
        low_ci_gaussian = sum_gauss[slope_name, "Q2.5"],
        high_ci_gaussian = sum_gauss[slope_name, "Q97.5"],
        corr_student = sum_student[slope_name, "Estimate"],
        low_ci_student = sum_student[slope_name, "Q2.5"],
        high_ci_student = sum_student[slope_name, "Q97.5"],
        nu_student = sum_student["nu", "Estimate"],
        stringsAsFactors = FALSE
      ))
    }
  }
  
  return(results)
}
