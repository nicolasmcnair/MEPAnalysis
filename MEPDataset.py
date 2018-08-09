from __future__ import division
import mepconfig
import numpy as np
from collections import deque
from csv import reader, writer
from datetime import datetime
from os.path import splitext
from sys import version_info
from scipy.signal import detrend
from struct import unpack_from
from plotmep import plot_data

class MEPDataset(object):

    PEAK_DETECTION_METHODS = {'HEIGHT', 'PROM'}
    TIME_WINDOW_DETECTION_METHODS = {'AUC', 'AVERAGE', 'RMS'}
    BACKGROUND_MOVEMENT_DETECTION_METHODS = {'PTP', 'RMS', 'BOTH'}

    @classmethod
    def extract_peaks_and_valleys(cls,
                                  data,
                                  boundary=(0,None),
                                  prominence=mepconfig.peak_prominence,
                                  prominence_correction_factor=mepconfig.prominence_correction_factor):
        """Return all peaks and valleys that exceed specified prominence, sorted by prominence"""

        def getProminences(data,peak_locs):
            """Get prominence all of peaks"""
            peak_prominences = np.zeros(len(peak_locs))
            for idx,peak in enumerate(peak_locs):
                # Find all higher peaks to the left of current peak
                l_peaks = np.where(data[peak_locs[:idx]] >= data[peak])[0]
                # Find the lowest point between the current peak and the closest higher peak, if any, or the left edge
                l_base = np.min(data[peak_locs[l_peaks[-1]] if len(l_peaks) else 0:peak])
                # Find all higher peaks to the right of current peak (correcting for truncated search)
                r_peaks = np.where(data[peak_locs[idx + 1:]] >= data[peak])[0] + (idx + 1)
                # Find the lowest point between the current peak and the closest higher peak, if any, or the right edge
                r_base = np.min(data[peak + 1:peak_locs[r_peaks[0]] if len(r_peaks) else None])
                # Calculate prominence for current peak
                peak_prominences[idx] = data[peak] - max(l_base,r_base)
            return peak_prominences

        # Set up return values
        peaks = valleys = None
        # Isolate time window
        sample_data = data[slice(*boundary)]
        # Calculate first derivative
        derivative = np.diff(sample_data)
        # Replace all zero points in diff (N.B. left point used as anchor when flat)
        zeros, = np.where(derivative == 0)
        while len(zeros):
            derivative[zeros] = np.hstack([derivative[1:], 0.])[zeros]
            zeros, = np.where(derivative == 0)
        # Detect peak locations from the derivative
        peak_locs = np.where((np.hstack([derivative, 0.]) < 0.) & (np.hstack([0., derivative]) > 0.))[0]

        # Get prominence of peaks (if any, otherwise bail now)
        if len(peak_locs):
            prominence_criterion = prominence
            peak_prominences = getProminences(sample_data,peak_locs)
            prominent_peaks = None
            while not prominent_peaks:
                prominent_peaks = [(peak_locs[prom_peak],prom_peak) for prom_peak in np.where(peak_prominences > prominence_criterion)[0]]
                prominence_criterion -= prominence * prominence_correction_factor
            peaks = [(peak_loc + boundary[0],sample_data[peak_loc],prom) for peak_loc,prom in prominent_peaks]
            # Detect locations from the derivative
            valley_locs = np.where((np.hstack([derivative, 0.]) > 0.) & (np.hstack([0., derivative]) < 0.))[0]
            valley_locs = valley_locs[np.where(valley_locs > peak_locs[0])[0]]
            # Get prominence of valleys (if any, otherwise clear peak_locs and bail now)
            if len(valley_locs):
                prominence_criterion = prominence
                valley_prominences = getProminences(sample_data * -1,valley_locs)
                prominent_valleys = None
                while not prominent_valleys:
                    prominent_valleys = [(valley_locs[prom_valley],prom_valley) for prom_valley in np.where(valley_prominences > prominence_criterion)[0]]
                    prominence_criterion -= prominence * prominence_correction_factor
                valleys = [(valley_loc + boundary[0],sample_data[valley_loc],prom) for valley_loc,prom in prominent_valleys]
            else:
                peaks = None
        return [peaks,valleys]

    @classmethod
    def extract_time_window(cls,
                            rect_data,
                            ptp_data,
                            boundary=(0,None),
                            method=mepconfig.time_window_detection,
                            search=True):
        #Set up return values
        window_amplitude = onset = offset = None
        # Get subsample of data
        sample_data = rect_data[slice(*boundary)]
        # Are we searching for onset/offset times?
        if search:
            # Calculate cumulative signal
            cumulative_signal = np.cumsum(sample_data)
            # Calculate reference line
            reference_line = np.interp(range(len(cumulative_signal)),[0,len(cumulative_signal)-1],[cumulative_signal[0],cumulative_signal[-1]])
            # Get integrated profile
            integrated_profile = reference_line - cumulative_signal
            # Estimating time windows is too problematic if the integrated profile does not have both positive and negative numbers
            if len(np.where(integrated_profile < 0)) and len(np.where(integrated_profile > 0)):
                # Max difference is MEP onset
                onset = np.argmax(integrated_profile)
                # Min difference (after onset) is MEP offset
                offset = np.argmin(integrated_profile[onset:]) + onset
        else:
            onset,offset = boundary
        if onset is not None:
            if method.upper() == 'RMS':
                # Get rms amplitude
                window_amplitude = np.sqrt(np.mean(np.square(sample_data[onset:offset + 1])))
            elif method.upper() == 'AVERAGE':
                window_amplitude = np.mean(sample_data[onset:offset + 1])
            elif method.upper() == 'AUC':
                window_amplitude = np.trapz(sample_data[onset:offset + 1],dx=1)
                    
            # Correct samples for starting point when returning
            onset += boundary[0]
            offset += boundary[0]

            # If movement onset/offset don't bound identified peaks, can't trust them
            if search and (not ptp_data[0] or (onset >= ptp_data[1] and offset <= ptp_data[2])):
                window_amplitude = onset = offset = None

        return [window_amplitude, onset, offset]

    def __init__(self,
                 mep_file,
                 n_trials,
                 polarity=mepconfig.polarity,
                 detrend_data=mepconfig.detrend):
        self.header = {key:None for key in ['n_channels','sample_rate','session_time','pretrigger_time','posttrigger_time','file']}
        self.header['n_trials'] = n_trials
        self.header['file'] = mep_file
        self.channels = []
        self._load_file(splitext(mep_file)[1], polarity, detrend_data)

    def _load_file(self,
                   file_type,
                   polarity=mepconfig.polarity,
                   detrend_data=mepconfig.detrend):
        if file_type == '.adibin':
            # Read entire adibin file (N.B. The impact on RAM will be more or less identical to the size of the file)
            with open(self.header['file'], 'rb') as mep_file:
                adibin_contents = mep_file.read()

            # Unpack and then parse the file header
            adibin_header = deque(unpack_from(mepconfig.adibin_file_header_fmt_string, adibin_contents, 0))
            # On valid LabChart binary files the first four bytes spell "CFWB"
            adibin_version = adibin_header.popleft()
            if isinstance(adibin_version, bytes):
                adibin_version = adibin_version.decode('ascii')
                convert_from_bytes = True
            else:
                convert_from_bytes = False
            if adibin_version != 'CFWB':
                raise ValueError('Incorrect file format.')
            adibin_header.popleft() # Dump version info
            self.header['sample_rate'] = adibin_header.popleft() * 1000
            self.header['session_time'] = str(datetime(*([adibin_header.popleft() for _ in range(5)] + [int(adibin_header.popleft()),])))
            self.header['pretrigger_time'] = max(adibin_header.popleft(),0) * 1000
            self.header['n_channels'] = adibin_header.popleft()
            n_samples = adibin_header.popleft()
            time_channel = bool(adibin_header.popleft())
            data_format = ('d', 'f', 'h')[adibin_header.popleft() - 1]
            del adibin_header

            # Validate number of trials
            self.header['samples_per_trial'],remainder = divmod(n_samples, self.header['n_trials'])
            if remainder:
                raise ValueError('Number of trials do not match number of samples.')
            self.header['posttrigger_time'] = self.sample_to_time(self.header['samples_per_trial'])

            # Read in channel headers
            offset = mepconfig.adibin_file_header_byte_length
            self.channels = [None] * self.header['n_channels']
            for channel in range(self.header['n_channels']):
                channel_header = deque(unpack_from("<64c4d", adibin_contents, offset))
                self.channels[channel] = {'header':{}}
                if convert_from_bytes:
                    self.channels[channel]['header']['title'] = ''.join([channel_header.popleft().decode('ascii') for _ in range(32)]).replace('\x00', '')
                    self.channels[channel]['header']['units'] = ''.join([channel_header.popleft().decode('ascii') for _ in range(32)]).replace('\x00', '')
                else:
                    self.channels[channel]['header']['title'] = ''.join([channel_header.popleft() for _ in range(32)]).replace('\x00', '')
                    self.channels[channel]['header']['units'] = ''.join([channel_header.popleft() for _ in range(32)]).replace('\x00', '')
                # Adjust scale to incoporate voltage correction (so output is in microvolts)
                self.channels[channel]['header']['scale'] = channel_header.popleft() * [1000000,1000,1][(['V', 'mV', '?V'].index(self.channels[channel]['header']['units']))]
                self.channels[channel]['header']['offset'] = channel_header.popleft()
                self.channels[channel]['header']['ptp'] = self.channels[channel]['header']['time_window'] = self.channels[channel]['header']['rejected'] = False
                self.channels[channel]['header']['time_window_method'] = None
                # N.B. The last two elements in channel_header (rangeHigh & rangeLow) are unused
                offset += mepconfig.adibin_channel_header_byte_length
                del channel_header

            # Read in all data
            adibin_data = unpack_from("<" + str((self.header['n_channels'] + time_channel) * n_samples) + data_format,adibin_contents,offset)
            del adibin_contents

            # Extract interleaved channel data (time channel is first, if present)
            for channel in range(self.header['n_channels']):
                # Create data array
                self.channels[channel]['data'] = np.array([(self.channels[channel]['header']['scale'] * (x + self.channels[channel]['header']['offset'])) * polarity for x in adibin_data[channel + time_channel::self.header['n_channels'] + time_channel]], dtype=np.float).reshape(self.header['n_trials'], self.header['samples_per_trial'])
                if detrend_data:
                    self.channels[channel]['data'] = detrend(self.channels[channel]['data'],type=detrend_data)
                self.channels[channel]['rect'] = np.fabs(self.channels[channel]['data'])
                self.channels[channel]['ptp'] = self.channels[channel]['time_window'] = [None,None,None]
                self.channels[channel]['rejected'] = None
            del adibin_data

        # Read ASCII file
        elif file_type == '.csv':

            # Helper function
            def is_number(s):
                try:
                    float(s)
                    return True
                except ValueError:
                    return False

            # Read in header line
            with open(self.header['file'], 'rb') as mep_file:
                csv_reader = reader(mep_file)
                header_line = next(csv_reader)
            # Get start time (this assumes first digit entry in headerline is onset time, which it should be if it's a csv file written by this module)
            time_index,self.header['pretrigger_time'] = next((x,float(y) * -1) for x,y in enumerate(header_line) if is_number(y))
            # Work out number of samples
            self.header['samples_per_trial'] = len(header_line) - time_index
            # Work out sampling rate
            self.header['sample_rate'] = (float(header_line[-1]) - float(header_line[time_index])) / (self.header['samples_per_trial'] - 1)
            # Work out posttrigger time
            self.header['posttrigger_time'] = self.sample_to_time(self.header['samples_per_trial'])
            # Get PTP index, if present
            ptp_index = header_line.index('PTP') if 'PTP' in header_line else None
            # Get Time Window index, if present
            window_index = header_line.index('Time Window') if 'Time Window' in header_line else None
            # Get Background movement rejection index, if present
            rejected_index = header_line.index('Rejected') if 'Rejected' in header_line else None
            # Get channel indices
            channel_column = np.loadtxt(self.header['file'], dtype=np.str, delimiter=',',skiprows=1, usecols=[0])
            self.header['n_channels'] = len(set(channel_column))
            channel_indices = sorted([np.where(channel_column==item)[0][0] for item in set(channel_column)])
            self.header['n_trials'] = len(channel_column) - channel_indices[-1]
            channel_indices = zip(channel_indices,channel_indices[1:] + [len(channel_column)])
            # Get data
            ptp_data = window_data = rejected_data = None
            data = np.genfromtxt(self.header['file'], dtype=np.float, delimiter=',', skip_header=1, usecols=range(time_index, len(header_line)))
            if ptp_index:
                ptp_data = np.genfromtxt(self.header['file'], dtype=np.float, delimiter=',', skip_header=1, usecols=range(ptp_index, ptp_index + 3), filling_values=None)
            if window_index:
                window_data = np.genfromtxt(self.header['file'], dtype=np.float, delimiter=',', skip_header=1, usecols=range(window_index, window_index + 3), filling_values=None)
            if rejected_index:
                rejected_data = np.genfromtxt(self.header['file'], dtype=np.float, delimiter=',', skip_header=1, usecols=[rejected_index], filling_values=False)
            # Sort through channels
            for channel in range(self.header['n_channels']):
                self.channels.append({'header':{}})
                self.channels[channel]['header']['title'] = channel_column[channel_indices[channel][0]]
                self.channels[channel]['header']['ptp'] = self.channels[channel]['header']['time_window'] = self.channels[channel]['header']['rejected'] = False
                self.channels[channel]['header']['time_window_method'] = None
                self.channels[channel]['ptp'] = self.channels[channel]['time_window'] = [None,None,None]
                self.channels[channel]['rect'] = self.channels[channel]['rejected'] = None
                # Read in channel data
                self.channels[channel]['data'] = data[slice(*channel_indices[channel])]
                self.channels[channel]['rect'] = np.fabs(self.channels[channel]['data'])
                if ptp_index:
                    self.channels[channel]['header']['ptp'] = True
                    self.channels[channel]['ptp'] = ptp_data[slice(*channel_indices[channel])]
                if window_index:
                    self.channels[channel]['header']['time_window'] = True
                    self.channels[channel]['time_window'] = window_data[slice(*channel_indices[channel])]
                if rejected_index:
                    self.channels[channel]['header']['rejected'] = True
                    self.channels[channel]['rejected'] = rejected_data[slice(*channel_indices[channel])]

    def _parse_boundary(self,boundary):
        return (0 if boundary[0] is None else self.time_to_sample(min((boundary[0],-self.header['pretrigger_time']),key=abs)),
                None if boundary[1] is None else self.time_to_sample(min((boundary[1] + self.header['sample_rate'],self.header['posttrigger_time']),key=abs)))

    def time_to_sample(self, time_point):
        return int((time_point + self.header['pretrigger_time']) // self.header['sample_rate'])

    def sample_to_time(self, sample_point):
        return (sample_point * self.header['sample_rate']) - self.header['pretrigger_time']

    def detect_background_movement(self,
                                   method=mepconfig.background_detection,
                                   boundary = mepconfig.background_boundary,
                                   peak_prominence=mepconfig.background_peak_prominence,
                                   ptp_threshold = mepconfig.background_ptp_threshold,
                                   max_ptp_interval = mepconfig.peak_max_ptp_interval,
                                   rms_threshold = mepconfig.background_rms_threshold):
        # Detection must be 'PTP', 'RMS', or 'BOTH'
        if method.upper() not in MEPDataset.BACKGROUND_MOVEMENT_DETECTION_METHODS:
            raise ValueError('Invalid background movement detection value.')
        # Convert boundary for background movement detection and max ptp interval to samples
        boundary = self._parse_boundary(boundary)
        max_ptp_interval = int(max_ptp_interval / self.header['sample_rate'])
        # Detect background movement
        for channel in self.channels:
            channel['header']['rejected'] = True
            channel['rejected'] = [False] * self.header['n_trials']
            for trial in range(self.header['n_trials']):
                # Handle PTP detection
                if method.upper() in {'PTP','BOTH'}:
                    peaks, valleys = MEPDataset.extract_peaks_and_valleys(channel['data'][trial], boundary, prominence=peak_prominence)
                    # Marry peaks to valleys and get peak-to-peak values
                    for peak in peaks:
                        peak_boundary = max_ptp_interval + peak[0]
                        # Get valleys
                        valley = [x[0] for x in valleys if  peak[0] < x[0] <= peak_boundary]
                        # If we got a valley, check that the peak-to-peak exceeds the threshold
                        if valley and ptp_threshold < (channel['data'][trial][peak[0]] - channel['data'][trial][valley[0]]):
                            channel['rejected'][trial] = True
                            break  
                # Handle RMS detection
                if method.upper() in {'RMS','BOTH'} and not channel['rejected'][trial]:
                    channel['rejected'][trial] = rms_threshold < MEPDataset.extract_time_window(channel['rect'][trial], channel['data'][trial], boundary, 'RMS', search=False)[0]

    def analyse_time_window(self,
                            method=mepconfig.time_window_detection,
                            boundary = mepconfig.time_window_boundary):
        # Detection must be 'AUC', 'AVERAGE', or 'RMS'
        if method.upper() not in MEPDataset.TIME_WINDOW_DETECTION_METHODS:
            raise ValueError('Invalid time window detection value. Must be: AUC, Average, or RMS.')
        # Do we have peak-to-peak information?
        if self.header['n_channels'] and not all([channel['header']['ptp'] for channel in self.channels]):
            self.analyse_peak_to_peak(secondary_to_time_window=True)
        # Convert boundary for peak detection to samples
        boundary = self._parse_boundary(boundary)
        # Create rectified data and then extract time window info for each trial (but only if we detected movement - i.e., have a ptp value)
        for channel in self.channels:
            channel['header']['time_window'] = True
            channel['header']['time_window_method'] = method.upper()
            channel['time_window'] = [MEPDataset.extract_time_window(channel['rect'][trial], channel['ptp'][trial], boundary, method) if channel['ptp'][trial][0] else [None,None,None] for trial in range(self.header['n_trials'])]

    def analyse_peak_to_peak(self,
                             method=mepconfig.peak_detection,
                             boundary = mepconfig.peak_time_window,
                             max_ptp_interval=mepconfig.peak_max_ptp_interval,
                             secondary_to_time_window=False):
        # Detection must be 'HEIGHT' or 'PROM'
        if method.upper() == 'HEIGHT':
            def sort_key_function(tup): return tup[1]
            peak_sort_reverse = True
            valley_sort_reverse = False
        elif method.upper() == 'PROM':
            def sort_key_function(tup): return tup[2]
            peak_sort_reverse = valley_sort_reverse = True
        else:
            raise ValueError('Invalid peak detection value. Must be: Height or Prom.')
        # Convert boundary for peak detection to samples
        boundary = self._parse_boundary(boundary)
        # Extract PTP data
        max_ptp_interval = max_ptp_interval // self.header['sample_rate']
        for channel in self.channels:
            channel['header']['ptp'] = False if secondary_to_time_window else True
            channel['ptp'] = [None] * self.header['n_trials']
            for trial in range(self.header['n_trials']):
                peaks,valleys = MEPDataset.extract_peaks_and_valleys(channel['data'][trial], boundary)
                p1 = p2 = None
                # Are there any peaks?
                if peaks:
                    # Sort by relevant attribute
                    peaks.sort(key=sort_key_function, reverse=peak_sort_reverse)
                    p1 = peaks[0][0]
                # Did we find any valleys? Only bother if we also found p1
                if p1 and valleys:
                    # Make sure our valley occurs after p1
                    peak_boundary = max_ptp_interval + p1
                    valleys = [x for x in valleys if p1 < x[0] <= peak_boundary] or [x for x in valleys if p1 < x[0]]
                    if valleys:
                        valleys.sort(key=sort_key_function, reverse=valley_sort_reverse)
                        p2 = valleys[0][0]
                    else:
                        p1 = None
                # Calculate peak-to-peak value
                channel['ptp'][trial] = [channel['data'][trial][p1] - channel['data'][trial][p2],p1,p2] if p1 and p2 else [None,None,None]

    def write_to_file(self,
                      boundary=mepconfig.output_time_window):
        # Make sure we have channel data
        if self.header['n_channels']:
            # Get MEP header line
            header_string = ['Channel','Trial']
            if any([channel['header']['ptp'] for channel in self.channels]):
                header_string += ['PTP','P1','P2']
            if any([channel['header']['time_window'] for channel in self.channels]):
                index = [channel['header']['time_window'] for channel in self.channels].index(True)
                header_string += [self.channels[index]['header']['time_window_method'].upper(),'Onset','Offset']
            if any([channel['header']['rejected'] for channel in self.channels]):
                header_string += ['Rejected']

            # Get output times
            start_time = boundary[0]
            end_time = boundary[1] or self.sample_to_time(self.header['samples_per_trial'])
            header_string += [str(start_time + (x * self.header['sample_rate'])) for x in range(int((end_time - start_time) / self.header['sample_rate']) + 1)]
            boundary = self._parse_boundary(boundary)

            # Write to file
            with (open(splitext(self.header['file'])[0] + '.csv','w',newline='') if version_info[0] > 2 else open(splitext(self.header['file'])[0] + '.csv','wb')) as output_file:
                csv_writer = writer(output_file)
                # Write header 
                csv_writer.writerow(header_string)    
                # Write MEP data
                for channel in self.channels:
                    for trial in range(self.header['n_trials']):
                        trial_string = [channel['header']['title'], trial + 1]
                        if channel['header']['ptp']:
                            if channel['ptp'][trial][0] and (not channel['header']['rejected'] or not channel['rejected'][trial]):
                                trial_string += [channel['ptp'][trial][0]] + [self.sample_to_time(x) for x in channel['ptp'][trial][1:]]
                            else:
                                trial_string += ['','','']
                        if channel['header']['time_window']:
                            if channel['time_window'][trial][0] and (not channel['header']['rejected'] or not channel['rejected'][trial]):
                                trial_string += [channel['time_window'][trial][0]] + [self.sample_to_time(x) for x in channel['time_window'][trial][1:]]
                            else:
                                trial_string += ['','','']
                        if channel['header']['rejected']:
                            trial_string += ['Yes'] if channel['rejected'][trial] else ['-']
                        csv_writer.writerow(trial_string + [channel['data'][trial][x] for x in range(*boundary)])
        else:
            raise IndexError('No channel information found.')
    
    def query_data(self,query_type=None,reference_line=False):
        # If no query type provided, then query whatever exists
        if query_type is None:
            query_type = ['ptp'] if any([channel['header']['ptp'] for channel in self.channels]) else []
            query_type += ['time_window'] if any([channel['header']['time_window'] for channel in self.channels]) else []
        # Pass to plot_data
        plot_data(self,query_type,reference_line)
        # If we queried the ptp data, then recalculate the ptp values
        if 'ptp' in query_type:
            for channel in self.channels:
                if channel['header']['ptp']:
                    channel['ptp'] = [[channel['data'][trial][p1]-channel['data'][trial][p2],p1,p2] if None not in (p1,p2) else [None,None,None] for trial,(_,p1,p2) in enumerate(channel['ptp'])]
        # If we queried the time window data, then recalculate the time window values
        if 'time_window' in query_type:
            for channel in self.channels:
                if channel['header']['time_window']:
                    if channel['header']['time_window_method'].upper() == 'RMS':
                        channel['time_window'] = [[np.sqrt(np.mean(np.square(channel['rect'][trial][p1:p2 + 1]))),p1,p2] if None not in (p1,p2)  else [None,None,None] for trial,(_,p1,p2) in enumerate(channel['time_window'])]
                    if channel['header']['time_window_method'].upper() == 'AVERAGE':
                        channel['time_window'] = [[np.mean(channel['rect'][trial][p1:p2 + 1]),p1,p2] if None not in (p1,p2)  else [None,None,None] for trial,(_,p1,p2) in enumerate(channel['time_window'])]
                    if channel['header']['time_window_method'].upper() == 'AUC':
                        channel['time_window'] = [[np.trapz(channel['rect'][trial][p1:p2 + 1],dx=1),p1,p2] if None not in (p1,p2)  else [None,None,None] for trial,(_,p1,p2) in enumerate(channel['time_window'])]