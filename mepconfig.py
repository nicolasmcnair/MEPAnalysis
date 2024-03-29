# The polarity of the MEP deflection. Set to -1 for Negative-Positive MEP, or 1 for Positive-Negative MEP
polarity = -1
# Whether to apply a detrend to the data. Allowed values: 'linear_baseline', 'linear', 'constant', 'baseline', or None (N.B. baseline will use the background_boundary definition)
detrend = 'linear_baseline'
# Output file time window (in milliseconds). End-inclusive
# To start the time window from the beginning of the data, set the first value to None
# To end the time window at the end of the data, set the second value to None
output_time_window = (-50, 20)#(-100, 100)
# Type of baseline movement detection: 'SD', 'Voltage', or 'Both'
background_detection = 'both'
# Beginning and end time (in milliseconds) to look for baseline movement. End-inclusive
# To start the time window from the beginning of the data, set the first value to None
# To end the time window at the end of the data, set the second value to None
background_boundary = (-100,-5)
# Threshold for RMS movement detection (in sds from mean)
background_rms_outlier_sds = 3#6
# Threshold for background EMG voltage detection (in microVolts)
background_voltage_threshold = 50#100
# Peak analysis: 'Height' or 'Prom'
peak_detection = 'height'
# Beginning and end time (in milliseconds) to look for peaks. End-inclusive
# To start the time window from the beginning of the data, set the first value to None
# To end the time window at the end of the data, set the second value to None
peak_time_window = (15,50)
# Maximum time allowed between P1 and P2; if no suitable P2 is identified it will re-look without this limit; set to float('Inf') to ignore
peak_max_ptp_interval = 40
# Minimum peak prominence (in microVolts) to qualify for retention (will iterate down if no peaks/valleys are found); set to zero to ignore
# N.B. This is also used for Height analysis
peak_prominence = 50
# Time window analysis: 'Auc' (Area under curve), 'Average', or 'RMS'
time_window_detection = 'auc'
# Beginning and end time (in milliseconds) to look for peaks. End-inclusive
# To start the time window from the beginning of the data, set the first value to None
# To end the time window at the end of the data, set the second value to None
time_window_boundary = (5,None)
# Minimum prominence for peak-detection will reduce by this factor each time detection fails                
prominence_correction_factor = 0.2
# Type of bad mep detection: 'SD', 'Voltage', or 'Both'
bad_mep_detection_method = 'SD'
# Threshold for SD mep classification (in sds from mean)
bad_mep_outlier_sds = 3
# Threshold for minimum mep amplitude (in microVolts)
bad_mep_min_voltage = 50
# Remove MEP values for rejected trials
remove_rejected_mep = True
# Adibin header parameters
adibin_file_header_fmt_string = '<4sld5l2d4l'
adibin_file_header_byte_length = 68
adibin_channel_header_byte_length = 96