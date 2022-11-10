import MEPDataset, mepconfig
import numpy as np
from plotmep import queryManager
from scipy.signal import detrend

class mepQueryManager(queryManager):

    def __init__(self,*args,**kwargs):
        super(mepQueryManager,self).__init__(*args,**kwargs)
        # This is a guard, because when the figure plot is closed __init__ is called again (not sure why) which would reset 'accept' back to True
        if not hasattr(self, 'accept'): 
            self.accept = True

    def key_release(self,event):
        if event.key in {'y','n'}:
            if event.key == 'y':
                self.accept = True
            elif event.key == 'n':
                self.accept = False
            self.close()
        elif event.key == 'escape':
            super(mepQueryManager,self).key_release(event)

class MEP(MEPDataset.MEPDataset):

    def __init__(self,
                 sampling_rate,
                 n_samples,
                 zero_sample,
                 units='?V',
                 processing_steps=(MEPDataset.MEPDataset.detect_background_movement,
                                   MEPDataset.MEPDataset.analyse_peak_to_peak,
                                   MEPDataset.MEPDataset.analyse_time_window)):

        # Create pseudo MEPDataset that has one channel and one trial
        self.processing_steps = processing_steps
        self.header = {key:None for key in ['n_channels','sample_rate','session_time','pretrigger_time','posttrigger_time','file']}
        self.header['n_trials'] = 1
        self.header['n_channels'] = 1
        self.header['sample_rate'] = 1000 / sampling_rate
        self.header['samples_per_trial'] = n_samples
        self.header['session_time'] = None
        ms_per_sample = 1000 / sampling_rate
        pretrigger_samples = len(list(range(n_samples))[:zero_sample])
        posttrigger_samples = len(list(range(n_samples))[zero_sample:])
        self.header['pretrigger_time'] = pretrigger_samples * ms_per_sample
        self.header['posttrigger_time'] = posttrigger_samples * ms_per_sample
        self.channels = [{'header':{'title' : 'MEP',
                                    'units' : units,
                                    'scale' : [1000000,1000,1,1][(['V', 'mV', 'µV', '?V'].index(units))],
                                    'offset' : 0.0,
                                    'rejected' : {'background':False,'mep':False,'other':False},
                                    'ptp' : False,
                                    'time_window' : False},
                          'data' : [],
                          'rect' : [],
                          'ptp' : [None,None,None],
                          'time_window' : [None,None,None],
                          'rejected' : None}]

    def evaluate(self,new_data,show_mep=True,evaluate_using='ptp'):
        # Load new data
        self.channels[0]['data'] = [np.array(new_data)]
        # Detrend if necessary
        mepconfig.background_boundary = [-50,-5]
        if mepconfig.detrend:
            if mepconfig.detrend == 'baseline':
                self.channels[0]['data'][0] -= np.mean(self.channels[0]['data'][0][slice(*[self.time_to_sample(x) for x in mepconfig.background_boundary])])
            else:
                self.channels[0]['data'][0] = detrend(self.channels[0]['data'][0],type=mepconfig.detrend)
        # Calculate rectified data
        self.channels[0]['rect'] = [np.fabs(self.channels[0]['data'][0])]
        # Reset rejection information
        self.channels[0]['ptp'] = self.channels[0]['time_window'] = [None,None,None]
        self.channels[0]['rejected'] = None
        self.channels[0]['header']['rejected'] = {'background':False,'mep':False,'other':False}
        self.channels[0]['header']['ptp'] = False
        self.channels[0]['header']['time_window'] = False
        # Run processing steps on new data
        for process in self.processing_steps:
            process(self)
        # Evaluate using mepQueryManager if 
        if show_mep:
            accepted = mepQueryManager(self.channels[0],self.header,evaluate_using).accept
        else:
            accepted = any(list(self.channels[0]['header']['rejected'].values()))
        if not evaluate_using == 'ptp':
            evaluate_using = 'time_window'
        mep_data = self.channels[0][evaluate_using]
        return (accepted,mep_data)