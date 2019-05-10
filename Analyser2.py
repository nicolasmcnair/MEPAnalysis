from MEPDataset import MEPDataset
try:
    from Tkinter import Tk
    from tkFileDialog import askopenfilename
    from tkSimpleDialog import askinteger 
except ModuleNotFoundError:
    from tkinter import Tk
    from tkinter.filedialog import askopenfilename
    from tkinter.simpledialog import askinteger 


# Either set the number of trials here, or set to None to have a dialog query the number when you load a file
specified_n_trials = None

#############################           DON'T EDIT BETWEEN THESE LINES         ###################################
Tk().withdraw()                                                                                                  #
file_name = askopenfilename(filetypes = [('LabChart binary files', '*.adibin'),('Processed MEP file','*.csv')])  #
if not file_name:                                                                                                #
    exit()                                                                                                       #
n_trials = specified_n_trials or askinteger('Number of trials','Enter the number of trials',minvalue=1)          #
if not n_trials:                                                                                                 #
    exit()                                                                                                       #
mep_dataset = MEPDataset(file_name,n_trials)                                                                     #
#############################           DON'T EDIT BETWEEN THESE LINES         ###################################

# Comment in/out the various analyses you'd like to conduct here
# These will use the parameters outlined in mepconfig.py by default
# However, you can also change them by passing them as arguments in the functions below
# See the MEPDataset.py files and plotmep.py files for details on the various arguments
mep_dataset.detect_background_movement()
mep_dataset.analyse_peak_to_peak()
mep_dataset.detect_bad_meps()
#mep_dataset.analyse_time_window()
mep_dataset.query_data()
mep_dataset.write_to_file()
