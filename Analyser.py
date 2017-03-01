from MEPDataset import MEPDataset
from Tkinter import Tk
from tkFileDialog import askopenfilename
from tkSimpleDialog import askinteger

specified_n_trials = 308

Tk().withdraw()
file_name = askopenfilename(filetypes = [('LabChart binary files', '*.adibin'),('Processed MEP file','*.csv')])
if not file_name:
    exit()
n_trials = specified_n_trials or askinteger('Number of trials','Enter the number of trials',minvalue=1)
if not n_trials:
    exit()
mep_dataset = MEPDataset(file_name,n_trials)
mep_dataset.detect_background_movement()
mep_dataset.analyse_peak_to_peak()
mep_dataset.analyse_time_window()
mep_dataset.query_data()
mep_dataset.write_to_file()