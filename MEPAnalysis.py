from __future__ import division
from Tkinter import Tk
from tkFileDialog import askopenfilename
from os.path import splitext
from sys import exit
from numpy import trapz, interp, mean
from operator import itemgetter
from csv import writer

# User-Defined Variables
processing = {}
processing['totalNTrials'] = 198                     # Total number of trials; this must evenly divide with the number of samples
processing['skipNTrials'] = 0                        # Skip this many trials at the beginning of the file
processing['polarity'] = -1                          # 1 = Postive-Negative MEP deflection; -1 = Negative-Positive MEP deflection
processing['detrend'] = 'linear'                     # Detrend EMG data; allowable values are 'linear' (least-squares regression), 'constant' (mean), or None
processing['timeWindow'] = (-50,100)                 # Beginning and end time (in milliseconds) for output .csv file (set to None to include entire EMG trace)
processing['matchToBehavioural'] = False             # Set to True to if you want to merge the MEP data with a tab-delimited behavioural data file; False otherwise (N.B. Empty line is used to signify beginning of trial data)
processing['peakAnalysis'] = {}
processing['peakAnalysis']['analyse'] = 1            # Peak analysis: 0 = No, 1 = Max Height, 2 = Max Prominence, 3 = First identifiable peak, plus most prominent of subequent two valleys (will only occur if 'auto' or 'query' set to True)
processing['peakAnalysis']['timeWindow'] = (15,50)   # Beginning and end time (in milliseconds) to look for peaks (set to None to look at entire EMG after trigger onset)
processing['peakAnalysis']['maxPeakInterval'] = 10   # Maximum time allowed between P1 and P2; if no suitable P2 is identified it will re-look without set limit (set to None to ignore)
processing['peakAnalysis']['prominence'] = 50        # Minimum peak prominence (in microVolts) to qualify for retention (will iterate down 20% if no peaks/valleys are found); N.B. This is used for Max Height analysis as well (set to zero to ignore)
# Movement onset and offset are determined based on an 'Integrated Profile' of the rectified EMG data (see: Allison G T 2003 Trunk muscle onset detection technique for EMG signals with ECG artefact J. Electromyogr. Kinesiol. 13 209–16)
processing['windowAnalysis'] = {}
processing['windowAnalysis']['analyse'] = 1          # Time window analysis: 0 = No, 1 = Area under curve, 2 = Average amplitude (will only occur if 'auto' or 'query' set to True)
processing['windowAnalysis']['startTime'] = 10       # Beginning time (in milliseconds) to look for MEP onset

########################################################################################################################################
#                                                                                                                                      #
#                                        BE CAREFUL CHANGING ANYTHING BELOW THIS LINE                                                  #
#                                                                                                                                      #
########################################################################################################################################

SAVECSV = True
PROMCORRECTION = 0.2
COMPUTEAVERAGE = True

# Get file
filename = behaviouralFile = fileHeader = channels = subjectInfo = trialInfo = None

Tk().withdraw()
filename = askopenfilename(filetypes = [('LabChart binary files', '*.adibin'),('Processed MEP file','*.csv'),('Processed Subject Datafile','*.txt')])
if not filename:
    exit()

# Some helper functions
def isNumber(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def convertTimeToSample(timePoint):
    return int((timePoint + fileHeader['pretriggerTime']) / fileHeader['samplingTickrate'])

def convertSampleToTime(samplePoint):
    return (samplePoint * fileHeader['samplingTickrate']) - fileHeader['pretriggerTime']

# Parse Labchart Binary file
def readAdibinFile(adibinFile):
    # Grab file contents
    with open(adibinFile,'rb') as labchartFile:
        labchartContents = labchartFile.read()

    # Read in the file header (using a deque as I'm lazy and it means I don't need to work out offsets as much)
    from collections import deque
    from struct import unpack,unpack_from
    labchartHeader = deque(unpack_from("<4sld5l2d4l",labchartContents,0))

    # On valid LabChart binary files the first four bytes spell "CFWB"
    if labchartHeader.popleft() != 'CFWB':
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
        fileHeader['samplesPerTrial'] = int(fileHeader['nSamples'] / processing['totalNTrials'])

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
    fileHeader['nTrials'] = (processing['totalNTrials'] - processing['skipNTrials'])
    if processing['detrend']:
        from scipy.signal import detrend
    for i in range(fileHeader['nChannels']):
        channelData = deque([(channels[i]['header']['scale'] * (x + channels[i]['header']['offset'])) * processing['polarity'] for x in labchartData[i + fileHeader['timeChannel']::fileHeader['nChannels'] + fileHeader['timeChannel']]])
        offset = 0
        #Set up channel trial data
        channels[i]['data'] = [None] * fileHeader['nTrials']
        channels[i]['ptp'] = [None] * fileHeader['nTrials']
        channels[i]['rect'] = [None] * fileHeader['nTrials']
        channels[i]['window'] = [None] * fileHeader['nTrials']
        # Cycle through each trial 
        for j in range(processing['skipNTrials'],processing['totalNTrials']):
            # Detrend data, if requested
            if processing['detrend']:
                channels[i]['data'][j] = list(detrend([channelData.popleft() for _ in range(fileHeader['samplesPerTrial'])],type=processing['detrend']))
            else:
                channels[i]['data'][j] = [channelData.popleft() for _ in range(fileHeader['samplesPerTrial'])]
            channels[i]['rect'][j] = map(abs,channels[i]['data'][j])
            channels[i]['ptp'][j] = [None,None,None]
            channels[i]['window'][j] = [None,None,None]

        if COMPUTEAVERAGE:
            temp = mean([x for x in channels[i]['data']],0)
            if processing['detrend']:
                channels[i]['data'].append(list(detrend(temp)))
            else:
                channels[i]['data'].append(list(temp))
            channels[i]['rect'].append(map(abs,channels[i]['data'][-1]))
            channels[i]['ptp'].append(None)
            channels[i]['window'].append(None)

    del labchartData, channelData

    return (fileHeader,channels)

# Parse csv mep file
def readMEPFile(mepFile):
    from csv import reader
    # Grab csv header
    with open(filename,'rb') as csvFile:
        csvReader = reader(csvFile)
        headerLine = csvReader.next()
        # Get start time (this assumes first number in headerline is onset time, which it should be if it's a csv file written by Analyser)
        fileHeader = {}
        timeIndex,fileHeader['pretriggerTime'] = next((x,float(y) * -1) for x,y in enumerate(headerLine) if isNumber(y))
        # Work out number of samples
        fileHeader['samplesPerTrial'] = len(headerLine) - timeIndex
        # Work out sampling rate
        fileHeader['samplingTickrate'] = (float(headerLine[-1]) - float(headerLine[timeIndex])) / (fileHeader['samplesPerTrial'] - 1)
        # Get PTP index, if present
        processing['peakAnalysis'] = {}
        ptpIndex = headerLine.index('PTP') if 'PTP' in headerLine else None
        processing['peakAnalysis']['analyse'] = 1 if ptpIndex else 0
        # Get Time Window index, if present
        processing['windowAnalysis'] = {}
        if [x for x in headerLine if x in {'AUC','Mean'}]:
            windowIndex = next(x for x,y in enumerate(headerLine) if y in {'AUC','Mean'})
            processing['windowAnalysis']['analyse'] = 1 if headerLine[windowIndex] == 'AUC' else 2
        else:
            windowIndex = None
            processing['windowAnalysis']['analyse'] = 0
        # Start reading in data
        fileHeader['nChannels'] = 0
        channels = {}
        channelIndex = -1
        currentChannel = ''
        for row in csvReader:
            if len(row) == len(headerLine):
                if currentChannel != row[0]:
                    currentChannel = row[0]
                    channelIndex += 1
                    fileHeader['nChannels'] += 1
                    channels[channelIndex] = {'header':{},'data':{}}
                    channels[channelIndex]['header']['title'] = currentChannel
                    channels[channelIndex]['data'] = []
                    channels[channelIndex]['ptp'] = []
                    channels[channelIndex]['rect'] = []
                    channels[channelIndex]['window'] = []
            
                if ptpIndex:
                    if isNumber(row[ptpIndex]):
                        channels[channelIndex]['ptp'].append([float(row[ptpIndex])] + [int((float(x) + fileHeader['pretriggerTime']) / fileHeader['samplingTickrate']) for x in row[ptpIndex + 1:ptpIndex + 3]])
                    else:
                        channels[channelIndex]['ptp'].append([None,None,None])
                if windowIndex:
                    if isNumber(row[windowIndex]):
                        channels[channelIndex]['window'].append([float(row[windowIndex])] + [int((float(x) + fileHeader['pretriggerTime']) / fileHeader['samplingTickrate']) for x in row[windowIndex + 1:windowIndex + 3]])
                    else:
                        channels[channelIndex]['window'].append([None,None,None])
                channels[channelIndex]['data'].append([float(x) for x in row[timeIndex:]])
                channels[channelIndex]['rect'].append(map(abs,channels[channelIndex]['data'][-1]))
    return (fileHeader,channels)

# Parse tab-delimited behavioural file
def readDataFile(dataFile,mepData = False):
    with open(dataFile,'r') as textFile:
        dataFileContents = textFile.readlines()

    dataFileContents = [x.rstrip() for x in dataFileContents]
    trialInfoIndex = dataFileContents.index('')
    subjectInfo = dataFileContents[:trialInfoIndex]
    trialInfo = dataFileContents[trialInfoIndex + 1:]
    del dataFileContents
    # Check if textfile contains MEP data, if so strip it out
    if mepData:
        headerLine = trialInfo[0].split('\t')
        headerLength = len(headerLine)
        fileHeader = {}
        # Get PTP index, if present        
        processing['peakAnalysis'] = {}
        ptpIndex = headerLine.index('PTP') if 'PTP' in headerLine else None
        processing['peakAnalysis']['analyse'] = 1 if ptpIndex else 0
        # Get Time Window index, if present
        processing['windowAnalysis'] = {}
        if [x for x in headerLine if x in {'AUC','Mean'}]:
            windowIndex = next(x for x,y in enumerate(headerLine) if y in {'AUC','Mean'})
            processing['windowAnalysis']['analyse'] = 1 if headerLine[windowIndex] == 'AUC' else 2
        else:
            windowIndex = None
            processing['windowAnalysis']['analyse'] = 0
        # Get start time (this is done after getting PTP and Time Window indices, since a behavioural datafile may have other numbers in the header line)
        timeIndex,fileHeader['pretriggerTime'] = next((x,float(y) * -1) for x,y in enumerate(headerLine) if isNumber(y) and x >= (windowIndex or ptpIndex or 0))
        # Work out number of samples
        fileHeader['samplesPerTrial'] = headerLength - timeIndex
        # Work out sampling rate
        fileHeader['samplingTickrate'] = (float(headerLine[-1]) - float(headerLine[timeIndex])) / (fileHeader['samplesPerTrial'] - 1)
        # Start reading in data
        channelIndex = headerLine.index('Channel') if 'Channel' in headerLine else None
        channels = {}
        currentChannel = ''
        # Do we have more than 1 channel?
        if channelIndex:
            fileHeader['nChannels'] = 0
            channelIndex = -1
        else:
            fileHeader['nChannels'] = 1
            channelIndex = 0
            channels[0] = {'header':{},'data':{}}
            channels[0]['header']['title'] = ''
            channels[0]['data'] = []
            channels[0]['ptp'] = []
            channels[0]['rect'] = []
            channels[0]['window'] = []
        del headerLine[(ptpIndex or windowIndex):]
        trialInfo[0] = '\t'.join(headerLine)
        # Cycle through each trial
        for i in range(1,len(trialInfo)):
            # Pass MEP data to channels
            tmp = trialInfo[i].split('\t')
            if len(tmp) == headerLength:
                if channelIndex and currentChannel != trialInfo[i][channelIndex]:
                    currentChannel = tmp[channelIndex]
                    channelIndex += 1
                    fileHeader['nChannels'] += 1
                    channels[channelIndex] = {'header':{},'data':{}}
                    channels[channelIndex]['header']['title'] = currentChannel
                    channels[channelIndex]['data'] = []
                    channels[channelIndex]['ptp'] = []
                    channels[channelIndex]['rect'] = []
                    channels[channelIndex]['window'] = []
            
                if ptpIndex:
                    if isNumber(tmp[ptpIndex]):
                        channels[channelIndex]['ptp'].append([float(tmp[ptpIndex])] + [int((float(x) + fileHeader['pretriggerTime']) / fileHeader['samplingTickrate']) for x in tmp[ptpIndex + 1:ptpIndex + 3]])
                    else:
                        channels[channelIndex]['ptp'].append([None,None,None])
                if windowIndex:
                    if isNumber(tmp[windowIndex]):
                        channels[channelIndex]['window'].append([float(tmp[windowIndex])] + [int((float(x) + fileHeader['pretriggerTime']) / fileHeader['samplingTickrate']) for x in tmp[windowIndex + 1:windowIndex + 3]])
                    else:
                        channels[channelIndex]['window'].append([None,None,None])
                channels[channelIndex]['data'].append([float(x) for x in tmp[timeIndex:]])
                channels[channelIndex]['rect'].append(map(abs,channels[channelIndex]['data'][-1]))
            # Strip MEP data from trial data
            del tmp[(ptpIndex or windowIndex):]
            trialInfo[i] = '\t'.join(tmp)

        return (subjectInfo,trialInfo,fileHeader,channels)
    else:
        return (subjectInfo,trialInfo)

# Check what kind of file we have
doAnalysis = False
fileType = splitext(filename)[1]
if fileType == '.adibin':
    doAnalysis = True
    # Are we matching this with a behavioural datafile?
    if processing['matchToBehavioural']:
        behaviouralFile = askopenfilename(filetypes = [('Behavioural Datafile','*.txt')])
        if not behaviouralFile:
            exit()
        subjectInfo, trialInfo = readDataFile(behaviouralFile)

    # Read adibin file
    fileHeader, channels = readAdibinFile(filename)
    if processing['matchToBehavioural'] and (fileHeader['nTrials'] != (len(trialInfo) - 1)):
        print('Trial numbers do not match between Labchart file (' + str(fileHeader['nTrials']) + ') and Behavioural data file (' + str(len(trialInfo) - 1) + ').')
        exit()

elif fileType == '.csv':
    # Read csv file
    fileHeader, channels = readMEPFile(filename)

elif fileType == '.txt':
    # Read txt file
    subjectInfo, trialInfo, fileHeader, channels = readDataFile(filename,True)

else:
    print('Unrecognised file type.')
    exit()

# Do automated peak and/or time window analysis
if doAnalysis:

    # Convert time variables from milliseconds to samples
    if processing['peakAnalysis']['timeWindow']:
        peakTimeIndices = [convertTimeToSample(x) for x in processing['peakAnalysis']['timeWindow']]
    else:
        peakTimeIndices = [convertTimeToSample(0), fileHeader['samplesPerTrial']]

    def getPeakProminence(data,peak,p,v):
        """Get prominence of current peak"""
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
        sampData = data[samplingBoundary[0]:samplingBoundary[1] + 1]    
        # Calculate first derivative
        dx = [sampData[x] - sampData[x-1] for x in range(1,len(sampData))]
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
        while not outPeaks and (prominenceCriterion > 0.00001):
            for peak in peaks:
                # get prominence
                peakProminence = getPeakProminence(sampData,peak,peaks,valleys)
                # return peak value, bases, and prominence (if above threshold) - adjusting return locations for restricted sampling
                if peakProminence >= prominenceCriterion:
                    outPeaks.append((peak + samplingBoundary[0],sampData[peak],peakProminence))
            prominenceCriterion -= (prominence * 0.1)
        # valleys only have relevance if there is an initial peak
        if outPeaks:
            # Get prominence of valleys that occur after first (retained) peak and prior to interval boundary
            prominenceCriterion = prominence
            while not outValleys and (prominenceCriterion >= 0.0):
                valleys = [x for x in valleys if (x + samplingBoundary[0]) > outPeaks[0][0]]
                for valley in valleys:
                    # get prominence (inverting data)
                    peakProminence = getPeakProminence([x * -1 for x in sampData],valley,valleys,peaks)
                    # return valley value, bases, and prominence (if above threshold) - adjusting return locations for restricted sampling
                    if (peakProminence >= prominenceCriterion):
                        outValleys.append((valley + samplingBoundary[0],sampData[valley],peakProminence))
                prominenceCriterion -= (prominence * PROMCORRECTION)
        return [outPeaks,outValleys]

    # Convert start time from milliseconds to samples
    startTimeIndex = convertTimeToSample(processing['windowAnalysis']['startTime']) 

    def extractWindow(data,startingPoint,peakInfo,windowType):
        #Set up return values
        windowAmplitude = onset = offset = None
        # Get subsample of data
        sampData = data[startingPoint:]
        # Calculate cumulative signal
        cumSignal = [sum(sampData[:x + 1]) for x in range(len(sampData))]
        # Calculate reference line
        refSignal = interp(range(len(cumSignal)),[0,len(cumSignal)-1],[cumSignal[0],cumSignal[-1]])
        # Get integrated profile
        integProf = [refSignal[x] - cumSignal[x] for x in range(len(cumSignal))]
        # Estimating time windows problematic if integrated profile does not have both positive and negative numbers
        if any(1 for _ in integProf if _ < 0) and any(1 for _ in integProf if _ > 0):
            # Max difference is MEP onset
            onset = max(enumerate(integProf), key=itemgetter(1))[0]
            # Min difference (after onset) is MEP offset
            offset = min(enumerate(integProf[onset:]), key=itemgetter(1))[0] + onset
            # Get window amplitude
            windowAmplitude = trapz(sampData[onset:offset + 1],dx=1)
            # Correct samples for starting point when returning
            onset += startingPoint
            offset += startingPoint
            # Do movement onset/offset bound identified peaks?
            if peakInfo[0] and onset < peakInfo[1] and offset > peakInfo[2]:
                # Average amplitude is based on area
                if windowType == 2:
                    windowAmplitude /= (offset - onset) + 1
            # If not, can't trust onset/offset times
            else:
                onset = offset = None

        return [windowAmplitude,onset,offset]

    for i in range(fileHeader['nChannels']):
        # Cycle through each trial    
        for j in range(len(channels[i]['data'])):
            # We always extract peak information, since it's used to bound time windows
            PEAKS = VALLEYS = p1 = p2 = None
            prominenceCriterion = processing['peakAnalysis']['prominence']
            # Extract peak and valley information
            PEAKS, VALLEYS = extractPeaksValleys(channels[i]['data'][j],peakTimeIndices,prominenceCriterion)
            # Did we find any peaks?
            if PEAKS:
                # Sort by relevant attribute                
                if processing['peakAnalysis']['analyse'] == 1:
                    PEAKS.sort(key=lambda tup: tup[1],reverse=True)
                elif processing['peakAnalysis']['analyse'] == 2:
                    PEAKS.sort(key=lambda tup: tup[2],reverse=True)
                elif processing['peakAnalysis']['analyse'] == 3:
                    PEAKS.sort(key=lambda tup: tup[0])
                # Set P1 to max peak and then strip out data and prominence values from peaks and convert to set
                p1 = PEAKS[0][0]
            # Did we find any valleys?
            if VALLEYS:
                # Make sure our valley occurs after p1
                if processing['peakAnalysis']['maxPeakInterval']:
                    peakBoundary = (processing['peakAnalysis']['maxPeakInterval'] / fileHeader['samplingTickrate']) + p1
                else:
                    peakBoundary = fileHeader['samplesPerTrial']
                VALLEYS = [x for x in VALLEYS if peakBoundary >= x[0] > p1] or [x for x in VALLEYS if x[0] > p1]
                if VALLEYS:
                    # Sort by relevant attribute                
                    if processing['peakAnalysis']['analyse'] == 1:
                        VALLEYS.sort(key=lambda tup: tup[1])
                    elif processing['peakAnalysis']['analyse'] == 2:
                        VALLEYS.sort(key=lambda tup: tup[2],reverse=True)
                    elif processing['peakAnalysis']['analyse'] == 3:
                        VALLEYS.sort(key=lambda tup: tup[0])
                        VALLEYS = VALLEYS[:2]
                        VALLEYS.sort(key=lambda tup: tup[2],reverse=True)
                    # Set P2 to max valley and then strip out data and prominence values from valleys and convert to set
                    p2 = VALLEYS[0][0]
            # Calculate peak-to-peak value
            if p1 and p2:
                channels[i]['ptp'][j] = [channels[i]['data'][j][p1] - channels[i]['data'][j][p2],p1,p2]
            del PEAKS, VALLEYS

            # Are we extracting time window information?
            if processing['windowAnalysis']['analyse']:
                # Only extract time window info if we found peaks (i.e., we detected movement)
                if channels[i]['ptp'][j][0]:
                    channels[i]['window'][j] = extractWindow(channels[i]['rect'][j],startTimeIndex,channels[i]['ptp'][j],processing['windowAnalysis']['analyse'])

# Query peaks, if analysing
if processing['peakAnalysis']['analyse']:
    from QueryMEP import queryData
    for i in range(fileHeader['nChannels']):        
        queryData(filename,channels[i]['header']['title'],channels[i]['data'],channels[i]['ptp'],fileHeader['samplingTickrate'],fileHeader['pretriggerTime'])

# Query time windows, if analysing
if processing['windowAnalysis']['analyse']:
    from QueryMEP import queryData
    for i in range(fileHeader['nChannels']):       
        queryData(filename,channels[i]['header']['title'],channels[i]['rect'],channels[i]['window'],fileHeader['samplingTickrate'],fileHeader['pretriggerTime'],processing['windowAnalysis']['analyse'])

# Do we at least have any channel data loaded?
if channels:
    # Write output file, if requested (only ever False for debugging)
    if SAVECSV:
        # Get samples for output times
        if processing['timeWindow']:
            timeIndices = [convertTimeToSample(x) for x in processing['timeWindow']]
        else:
            timeIndices = [0, fileHeader['samplesPerTrial']]

        # Get MEP header line
        headerString = []
        if not processing['matchToBehavioural'] or (fileHeader['nChannels'] > 1):
            headerString += ['Channel',]
        if not processing['matchToBehavioural']:
            headerString += ['Trial',]
        if processing['peakAnalysis']['analyse']:
            headerString += ['PTP','P1','P2']
        if processing['windowAnalysis']['analyse']:
            if processing['windowAnalysis']['analyse'] == 1:
                headerString += ['AUC',]
            elif processing['windowAnalysis']['analyse'] == 2:
                headerString += ['Average',]
            headerString += ['Onset','Offset']
        headerString += [str(convertSampleToTime(x)) for x in range(*timeIndices)]

        # Check if writing behavioural file or MEP file
        if processing['matchToBehavioural']:
            with open(splitext(filename)[0] + '_MEP.txt','w') as outputFile:
                # Write subject info
                for row in subjectInfo:
                    outputFile.write(row + '\n')
                outputFile.write('\n')
                # Write header line
                outputFile.write(trialInfo[0] + '\t' + '\t'.join(headerString) + '\n')
                # Write MEP data
                for i in range(fileHeader['nChannels']):
                    for j,row in enumerate(trialInfo[1:]):
                        if fileHeader['nChannels'] > 1:
                            trialString = [channels[i]['header']['title'],]
                        else:
                            trialString = []
                        if processing['peakAnalysis']['analyse']:
                            if channels[i]['ptp'][j][0]:
                                trialString += [str(channels[i]['ptp'][j][0])] + [str(convertSampleToTime(x)) for x in channels[i]['ptp'][j][1:]]
                            else:
                                trialString += ['','','']
                        if processing['windowAnalysis']['analyse']:
                            if channels[i]['window'][j][0]:
                                trialString += [str(channels[i]['window'][j][0])] + [str(convertSampleToTime(x)) for x in channels[i]['window'][j][1:]]
                            else:
                                trialString += ['','','']
                        outputFile.write(row + '\t' + '\t'.join(trialString[:] + [str(channels[i]['data'][j][x]) for x in range(*timeIndices)]) + '\n')

        else:
            with open(splitext(filename)[0] + '_MEP.csv','wb') as outputFile:
                csvWriter = writer(outputFile)
                # Write header 
                csvWriter.writerow(headerString)    
                # Write MEP data
                for i in range(fileHeader['nChannels']):
                    for j in range(len(channels[i]['data'])):
                        trialString = [channels[i]['header']['title'],j + 1,]
                        if processing['peakAnalysis']['analyse']:
                            if channels[i]['ptp'][j][0]:
                                trialString += [channels[i]['ptp'][j][0],] + [convertSampleToTime(x) for x in channels[i]['ptp'][j][1:]]
                            else:
                                trialString += ['','','']
                        if processing['windowAnalysis']['analyse']:
                            if channels[i]['window'][j][0]:
                                trialString += [channels[i]['window'][j][0],] + [convertSampleToTime(x) for x in channels[i]['window'][j][1:]]
                            else:
                                trialString += ['','','']
                        csvWriter.writerow(trialString + [channels[i]['data'][j][x] for x in range(*timeIndices)])
