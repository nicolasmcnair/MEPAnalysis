from __future__ import division
from Tkinter import Tk
from tkFileDialog import askopenfilename
from os.path import splitext
from sys import exit
from cPickle import dump

# User-Defined Variables
processing = {}
processing['totalNTrials'] = 180                     # Total number of trials; this must evenly divide with the number of samples
processing['skipNTrials'] = 0                        # Skip this many trials at the beginning of the file
processing['polarity'] = -1                          # 1 = Postive-Negative MEP deflection; -1 = Negative-Positive MEP deflection
processing['detrend'] = True                         # Detrend LabChart EMG data
processing['writeCSV'] = True                        # Write output to .csv file as well
processing['timeWindow'] = (-50,100)                 # Beginning and end time (in milliseconds) for output .csv file (set to None to include entire EMG)
processing['peakAnalysis'] = {}
processing['peakAnalysis']['analyse'] = 1            # Peak analysis: 0 = No, 1 = Max Height, 2 = Max Prominence, 3 = First Peak + most prominent of two following two valleys (will only occur if 'auto' or 'query' set to True)
processing['peakAnalysis']['auto'] = True            # Try automated peak analysis (will only occur if ['analyse'] set to something other than 0)
processing['peakAnalysis']['query'] = True           # Manually inspect peaks (will only occur if ['analyse'] set to something other than 0)
processing['peakAnalysis']['timeWindow'] = (15,50)   # Beginning and end time (in milliseconds) to look for MEP (set to None to look at entire EMG after trigger onset)
processing['peakAnalysis']['prominence'] = 100       # Minimum peak prominence (in microVolts) to qualify for retention (will iterate down 10% if no peaks/valleys are found); N.B. This is used for Max value analysis as well
processing['windowAnalysis'] = {}
processing['windowAnalysis']['analyse'] = 1          # Time window analysis: 0 = No, 1 = Area under curve, 2 = Average amplitude (will only occur if 'auto' or 'query' set to True)
processing['windowAnalysis']['auto'] = True          # Try automated time window analysis (will only occur if ['analyse'] set to something other than 0)
processing['windowAnalysis']['query'] = True         # Manually inspect time windows (will only occur if ['analyse'] set to something other than 0)
processing['windowAnalysis']['baseline'] = (-200,-5) # Beginning and end time (in milliseconds) to calculate baseline amplitude (set to None to include everything up to trigger onset)
processing['windowAnalysis']['startTime'] = 15       # Beginning time (in milliseconds) to look for MEP onset
processing['windowAnalysis']['sdThreshold'] = 8      # Standard deviation threshold for detecting MEP onset
processing['windowAnalysis']['nThreshold'] = 10      # Number of trials that signal must be above/below threshold for determining MEP onset/offset

########################################################################################################################################
#
#                                        BE CAREFUL CHANGING ANYTHING BELOW THIS LINE                                                  #
#
########################################################################################################################################

# Get file
Tk().withdraw()
filename = askopenfilename(filetypes = [('MEP analysis file','*.mep'),('LabChart binary files', '*.adibin')])
if not filename:
    exit()

# Check what kind of file we have
skipAnalysis = False
fileType = splitext(filename)[1]
if fileType == '.adibin':

    # Grab file contents
    with open(filename,'rb') as labchartFile:
        labchartContents = labchartFile.read()

    # Read in the file header (using a deque as I'm lazy and it means I don't need to work out offsets as much)
    from collections import deque
    from struct import unpack,unpack_from
    labchartHeader = deque(unpack_from("<4cld5l2d4l",labchartContents,0))

    # On valid LabChart binary files the first four bytes spell "CFWB"
    if ''.join([labchartHeader.popleft() for _ in range(4)]) != 'CFWB':
        print('Incorrect file format.')
        exit()

    # N.B. all times are converted to milliseconds
    from datetime import datetime
    fileHeader = {}
    fileHeader['version'] = labchartHeader.popleft()
    fileHeader['samplingTickrate'] = labchartHeader.popleft() * 1000
    fileHeader['sessionTime'] = str(datetime(*([labchartHeader.popleft() for _ in range(5)] + [int(labchartHeader.popleft()),])))
    fileHeader['pretriggerTime'] = max(labchartHeader.popleft(),0) * 1000
    fileHeader['nChannels'] = labchartHeader.popleft()
    fileHeader['nSamples'] = labchartHeader.popleft()
    fileHeader['timeChannel'] = bool(labchartHeader.popleft())
    fileHeader['dataFormat'] = ('d','f','h')[labchartHeader.popleft() - 1]
    del labchartHeader

    # Check number of trials is appropriate
    if fileHeader['nSamples'] % processing['totalNTrials']:
        print('Incorrect number of trials specified.')
        exit()
    else:
        samplesPerTrial = int(fileHeader['nSamples'] / processing['totalNTrials'])

    # Read in the channel headers (offset is byte 68; each channel header is 96 bytes)
    channels = {}
    offset = 68
    for i in range(fileHeader['nChannels']):
        channelHeader = deque(unpack_from("<64c4d",labchartContents,offset))
        channels[i] = {'header':{},'data':{}}
        channels[i]['header']['title'] = ''.join([channelHeader.popleft() for _ in range(32)]).replace('\x00','')
        channels[i]['header']['units'] = ''.join([channelHeader.popleft() for _ in range(32)]).replace('\x00','')
        # Adjust scale to incoporate voltage correction (so output is in microvolts)
        channels[i]['header']['scale'] = channelHeader.popleft() * [1000000,1000,1][(['V','mV','?V'].index(channels[i]['header']['units']))]
        channels[i]['header']['offset'] = channelHeader.popleft()
        channels[i]['header']['range_high'] = channelHeader.popleft()
        channels[i]['header']['range_low'] = channelHeader.popleft()
        offset += 96
    del channelHeader

    # Read in all data (accounting for time channel if present)
    labchartData = unpack_from("<" + str((fileHeader['nChannels'] + fileHeader['timeChannel']) * fileHeader['nSamples']) + fileHeader['dataFormat'],labchartContents,offset)
    del labchartContents

    # Extract channel data
    if processing['detrend']:
        from scipy.signal import detrend
    for i in range(fileHeader['nChannels']):
        channelData = deque([(channels[i]['header']['scale'] * (x + channels[i]['header']['offset'])) * processing['polarity'] for x in labchartData[i + fileHeader['timeChannel']::fileHeader['nChannels'] + fileHeader['timeChannel']]])
        offset = 0
        channels[i]['trialData'] = []
        # Cycle through each trial    
        for j in range(processing['skipNTrials'],processing['totalNTrials']):
            #Set up channel trial data
            channels[i]['trialData'].append({})
            # Detrend data, if requested
            if processing['detrend']:
                channels[i]['trialData'][j]['data'] = list(detrend([channelData.popleft() for _ in range(samplesPerTrial)]))
            else:
                channels[i]['trialData'][j]['data'] = [channelData.popleft() for _ in range(samplesPerTrial)]
            channels[i]['trialData'][j]['ptp'] = None
            channels[i]['trialData'][j]['peaks'] = set()
            channels[i]['trialData'][j]['valleys'] = set()
            channels[i]['trialData'][j]['window'] = None

    del labchartData, channelData
elif fileType == '.mep':
    from cPickle import load
    from tkMessageBox import askyesno
    # Grab file contents
    with open(filename,'rb') as mepFile:
        mepFileContents = load(mepFile)
    temp,fileHeader,channels = mepFileContents[:]
    temp['peakAnalysis']['auto'] = temp['peakAnalysis']['auto'] = True
    if askyesno('Reanalyse?','Apply current analysis parameters before continuing?'):
        temp['peakAnalysis'] = processing['peakAnalysis']
        temp['windowAnalysis'] = processing['windowAnalysis']
    else:
        skipAnalysis = True
    processing = temp
    del temp
else:
    print('Unrecognised file type.')
    exit()

# Load peak detection function(s) if required
if processing['peakAnalysis']['auto'] and processing['peakAnalysis']['analyse'] and not skipAnalysis:

    # Convert time variables from milliseconds to samples
    if processing['peakAnalysis']['timeWindow']:
        peakTimeIndices = [int((x + fileHeader['pretriggerTime']) / fileHeader['samplingTickrate']) for x in processing['peakAnalysis']['timeWindow']]
    else:
        peakTimeIndices = [int(fileHeader['pretriggerTime'] / fileHeader['samplingTickrate']), samplesPerTrial]

    def getPeakProminence(data,peak,p,v):
        """Get prominence of current peak"""
        from operator import itemgetter
        # Find all higher peaks to the left of current peak (default to left edge if none higher)
        higherPeaks = [x for x in p[:p.index(peak)][::-1] if data[x] > data[peak]] or [0]
        # If there are higher peaks, find the lowest valley between the current peak and the closest higher peak, otherwise default to the left edge
        leftBase = min([(x,data[x]) for x in v[::-1] if higherPeaks[0] < x < peak] or [(0,data[0])], key=itemgetter(1))[0]

        # Find all higher peaks to the right of current peak
        higherPeaks = [x for x in p[p.index(peak) + 1:] if data[x] > data[peak]] or [len(data) - 1]
        # If there are higher peaks, find the lowest valley between the current peak and the closest higher peak, otherwise default to the right edge
        rightBase = min([(x,data[x]) for x in v if peak < x < higherPeaks[0]] or [(len(data) - 1,data[-1])], key=itemgetter(1))[0]
        # get prominence, based on difference between highest base and peak
        return data[peak] - max(data[leftBase],data[rightBase])

    def extractPeaksValleys(data, samplingBoundary, prominence = processing['peakAnalysis']['prominence']):
        """Return all peaks and valleys that exceed specified prominence, sorted by prominence"""
        # Set up return values
        outPeaks = []
        outValleys = []
        # Isolate time window
        data = data[samplingBoundary[0]:samplingBoundary[1] + 1]    
        # Calculate first derivative
        dx = [data[x] - data[x-1] for x in range(1,len(data))]
        # Detect peak dx locations, N.B. add one to found locations to account for differentiatial
        peaks = []
        x = 0
        # Cycle through dx
        while x < len(dx[:-1]):
            # Only positive values of dx are candidates for peaks
            if dx[x] > 0:
                # Cycle through rest of dx
                for y in range(x + 1,len(dx)):
                    # If next value is flat, step over to next value
                    if dx[y] == 0:
                        x = y
                        continue
                    # Otherwise it's a change in elevation; if negative then we've found a peak
                    if dx[y] < 0:
                        peaks.append(x + 1)
                        x = y + 1
                    else:
                        x = y
                    # Go back to looking for peaks
                    break
            # If negative value step over to next value
            else:
                x += 1

        # Detect valley dx locations, N.B. add one to found locations to account for differentiatial
        valleys = []
        x = 0
        # Cycle through dx
        while x < len(dx[:-1]):
            # Only negative values of dx are candidates for valleys
            if dx[x] < 0:
                # Cycle through rest of dx
                for y in range(x + 1,len(dx)):
                    # If next value is flat, step over to next value
                    if dx[y] == 0:
                        x = y
                        continue
                    # Otherwise it's a change in elevation; if positive then we've found a valley
                    if dx[y] > 0:
                        valleys.append(x + 1)
                        x = y + 1
                    else:
                        x = y
                    # Go back to looking for peaks
                    break
            # If positive value step over to next value
            else:
                x += 1

        # Get prominence of peaks
        prominenceCriterion = prominence
        while not outPeaks and (prominenceCriterion >= 0.0):
            for peak in peaks:
                # get prominence
                peakProminence = getPeakProminence(data,peak,peaks,valleys)
                # return peak value, bases, and prominence (if above threshold) - adjusting return locations for restricted sampling
                if peakProminence >= prominenceCriterion:
                    outPeaks.append((peak + samplingBoundary[0],data[peak],peakProminence))
            prominenceCriterion -= (prominence * 0.1)
        # valleys only have relevance if there is an initial peak
        if outPeaks:
            # Get prominence of valleys that occur after first (retained) peak
            prominenceCriterion = prominence
            while not outValleys and (prominenceCriterion >= 0.0):
                for valley in [x for x in valleys if (x + samplingBoundary[0]) > outPeaks[0][0]]:
                    # get prominence (inverting data)
                    peakProminence = getPeakProminence([x * -1 for x in data],valley,valleys,peaks)
                    # return valley value, bases, and prominence (if above threshold) - adjusting return locations for restricted sampling
                    if (peakProminence >= prominenceCriterion):
                        outValleys.append((valley + samplingBoundary[0],data[valley],peakProminence))
                prominenceCriterion -= (prominence * 0.1)
        return [outPeaks,outValleys]

# Load time window detection function(s) if required
if processing['windowAnalysis']['auto'] and processing['windowAnalysis']['analyse'] and not skipAnalysis:
    from numpy import mean,std,trapz

    # Convert time variables from milliseconds to samples
    if processing['windowAnalysis']['baseline']:
        baselineTimeIndices = [int((x + fileHeader['pretriggerTime']) / fileHeader['samplingTickrate']) for x in processing['windowAnalysis']['baseline']]
    else:
        baselineTimeIndices = [0, int(fileHeader['pretriggerTime'] / fileHeader['samplingTickrate'])]
    startTimeIndex = int((processing['windowAnalysis']['startTime'] + fileHeader['pretriggerTime']) / fileHeader['samplingTickrate'])    

    def extractWindow(data,baselineBoundary,startingPoint,sdThreshold,nThreshold,windowType):
        # Rectify signal
        data = map(abs,data)
        # Calculate baseline threshold
        baselineThreshold = std(data[baselineBoundary[0]:baselineBoundary[1] + 1]) * sdThreshold
        # Start looking for values above threshold
        sample = startingPoint
        exceeded = 0
        onset = offset = windowAmplitude = None
        while sample < len(data):
            # If we find a value above threshold, add to current count above threshold
            if data[sample] > baselineThreshold:
                exceeded += 1
                # If we reach the threshold for number of samples over baseline, set onset and break out of while loop
                if exceeded == nThreshold:
                    onset = (sample - nThreshold) + 1
                    break
            # If we find a value below threshold, reset current count to zero
            else:
                exceeded = 0
            sample += 1
        #  If we found a valid onset, start looking for values below threshold (default to end of data if not found)
        if onset:
            offset = len(data)
            exceeded = 0
            sample += 1
            while sample < len(data):
                # If we find a value above threshold, add to current count above threshold
                if data[sample] < baselineThreshold:
                    exceeded += 1
                    # If we reach the threshold for number of samples over baseline, set onset and break out of while loop
                    if exceeded == nThreshold:
                        offset = (sample - nThreshold) + 1
                        break
                # If we find a value below threshold, reset current count to zero
                else:
                    exceeded = 0
                sample += 1
        # If onset found, get window amplitude
        if onset:
            if windowType == 1:
                windowAmplitude = trapz(data[onset:offset],dx=1)
            elif windowType == 2:
                windowAmplitude = mean(data[onset:offset])

        return [windowAmplitude,onset,offset]

# Do automated peak and/or time window analysis
if not skipAnalysis:
    for i in range(fileHeader['nChannels']):
        # Cycle through each trial    
        for j in range(len(channels[i]['trialData'])):
            # Are we extracting peak information?
            if processing['peakAnalysis']['auto'] and processing['peakAnalysis']['analyse']:
                p1 = p2 = None
                prominenceCriterion = processing['peakAnalysis']['prominence']
                # Extract peak and valley information
                channels[i]['trialData'][j]['peaks'], channels[i]['trialData'][j]['valleys'] = extractPeaksValleys(channels[i]['trialData'][j]['data'],peakTimeIndices,prominenceCriterion)
                # Did we find any peaks?
                if channels[i]['trialData'][j]['peaks']:
                    # Sort by relevant attribute                
                    if processing['peakAnalysis']['analyse'] == 1:
                        channels[i]['trialData'][j]['peaks'].sort(key=lambda tup: tup[1],reverse=True)
                    elif processing['peakAnalysis']['analyse'] == 2:
                        channels[i]['trialData'][j]['peaks'].sort(key=lambda tup: tup[2],reverse=True)
                    elif processing['peakAnalysis']['analyse'] == 3:
                        channels[i]['trialData'][j]['peaks'].sort(key=lambda tup: tup[0])
                    # Set P1 to max peak and then strip out data and prominence values from peaks and convert to set
                    p1 = channels[i]['trialData'][j]['peaks'][0][0]
                    channels[i]['trialData'][j]['peaks'] = set([x[0] for x in channels[i]['trialData'][j]['peaks']])
                # Did we find any valleys?
                if channels[i]['trialData'][j]['valleys']:
                    # Make sure our valley occurs after p1
                    p1Valleys = [x for x in channels[i]['trialData'][j]['valleys'] if x[0] > p1]
                    if p1Valleys:
                        # Sort by relevant attribute                
                        if processing['peakAnalysis']['analyse'] == 1:
                            p1Valleys.sort(key=lambda tup: tup[1])
                        elif processing['peakAnalysis']['analyse'] == 2:
                            p1Valleys.sort(key=lambda tup: tup[2],reverse=True)
                        elif processing['peakAnalysis']['analyse'] == 3:
                            p1Valleys.sort(key=lambda tup: tup[0])
                            p1Valleys = p1Valleys[:2]
                            p1Valleys.sort(key=lambda tup: tup[2],reverse=True)
                        # Set P2 to max valley and then strip out data and prominence values from valleys and convert to set
                        p2 = p1Valleys[0][0]
                    channels[i]['trialData'][j]['valleys'] = set([x[0] for x in channels[i]['trialData'][j]['valleys']])
                # Calculate peak-to-peak value
                if p1 and p2:
                    channels[i]['trialData'][j]['ptp'] = [channels[i]['trialData'][j]['data'][p1] - channels[i]['trialData'][j]['data'][p2],p1,p2]

            # Are we extracting time window information?
            if processing['windowAnalysis']['auto'] and processing['windowAnalysis']['analyse']:
                # Extract time window
                channels[i]['trialData'][j]['window'] = extractWindow(channels[i]['trialData'][j]['data'],baselineTimeIndices,startTimeIndex,processing['windowAnalysis']['sdThreshold'],processing['windowAnalysis']['nThreshold'],processing['windowAnalysis']['analyse'])

# Query peaks, if requested
if processing['peakAnalysis']['query'] and processing['peakAnalysis']['analyse']:
    from QueryMEP import queryPeaks
    for i in range(fileHeader['nChannels']):        
        queryPeaks(channels[i]['header']['title'],channels[i]['trialData'],fileHeader['samplingTickrate'],fileHeader['pretriggerTime'])

# Query time windows, if requested
if processing['windowAnalysis']['query'] and processing['windowAnalysis']['analyse']:
    from QueryMEP import queryWindow
    for i in range(fileHeader['nChannels']):        
        queryWindow(channels[i]['header']['title'],channels[i]['trialData'],fileHeader['samplingTickrate'],fileHeader['pretriggerTime'],processing['windowAnalysis']['analyse'])

# Save .mep file
with open(splitext(filename)[0] + '.mep','wb') as output_file:
    dump([processing,fileHeader,channels],output_file)

# Write output file, if requested
if processing['writeCSV']:
    from csv import writer
    with open(splitext(filename)[0] + '.csv','wb') as output_file:

        if processing['timeWindow']:
            timeIndices = [int((x + fileHeader['pretriggerTime']) / fileHeader['samplingTickrate']) for x in processing['timeWindow']]
        else:
            timeIndices = [0, samplesPerTrial]

        csvWriter = writer(output_file)
        headerString = ['Trial',]
        if processing['peakAnalysis']['analyse'] or processing['peakAnalysis']['query']:
            headerString += ['PTP','P1','P2']
        if processing['windowAnalysis']['analyse'] or processing['windowAnalysis']['query']:
            if processing['windowAnalysis']['analyse'] == 1:
                headerString += ['AUC',]
            elif processing['windowAnalysis']['analyse'] == 2:
                headerString += ['Average',]
            headerString += ['Onset','Offset']
        # Write header 
        csvWriter.writerow(headerString + [str(((x * fileHeader['samplingTickrate']) - fileHeader['pretriggerTime'])) for x in range(*timeIndices)])    

        for i in range(fileHeader['nChannels']):
            for j in range(len(channels[i]['trialData'])):
                trialString = [j + 1,]
                if processing['peakAnalysis']['analyse'] and (processing['peakAnalysis']['query'] or processing['peakAnalysis']['auto']):
                    if channels[i]['trialData'][j]['ptp']:
                        trialString += [channels[i]['trialData'][j]['ptp'][0],] + [(x * fileHeader['samplingTickrate']) - fileHeader['pretriggerTime'] for x in channels[i]['trialData'][j]['ptp'][1:]]
                    else:
                        trialString += ['-','-','-']
                if processing['windowAnalysis']['analyse'] and (processing['windowAnalysis']['query'] or processing['windowAnalysis']['auto']):
                    if channels[i]['trialData'][j]['window'] and channels[i]['trialData'][j]['window'][0]:
                        trialString += [channels[i]['trialData'][j]['window'][0],] + [(x * fileHeader['samplingTickrate']) - fileHeader['pretriggerTime'] for x in channels[i]['trialData'][j]['window'][1:]]
                    else:
                        trialString += ['-','-','-']
                csvWriter.writerow(trialString + [channels[i]['trialData'][j]['data'][x] for x in range(*timeIndices)])
