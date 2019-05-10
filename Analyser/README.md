# MEPAnalysis
Python script for analysing LabChart binary files containing MEP data. These must be in .adibin format, and you must know the number of trials in the file (you can check the original .adicht file in LabReader to determine this). It requires the numpy and scipy package libraries to run, and may have issues if not using Qt. The easiest way to ensure you have these is to use a Python distribution such as Anaconda. Using Anaconda, it has been tested as working on Windows and MacOSX with both Python 2 and 3.


Please check the mepconfig.py file and tailor to your particular paradigm before use. The overall analysis pipeline can be altered by commenting/uncommenting lines 30 through 34 in Analyser.py. 


Figure Controls:

&nbsp;&nbsp;&nbsp;&nbsp;Left click: Set P1/Onset

&nbsp;&nbsp;&nbsp;&nbsp;Right click: Set P2/Offset

&nbsp;&nbsp;&nbsp;&nbsp;Enter: Accept all trials (N.B. closing figure window will reject all changes, but continue with original points)

&nbsp;&nbsp;&nbsp;&nbsp;Delete: Remove peaks from current trial

&nbsp;&nbsp;&nbsp;&nbsp;Backspace: Toggle flag for background movement detection

&nbsp;&nbsp;&nbsp;&nbsp;Left/Right arrow: Cycle backward/forward through trials

&nbsp;&nbsp;&nbsp;&nbsp;Scroll Up/Down: Cycle backward/forward through trials

&nbsp;&nbsp;&nbsp;&nbsp;Escape: Quit entire MEP data analysis (this will mean no file is written)


Note: Depending on the Python version/distro you are running, you may need to alter line 8 in plotmep.py to a graphics backend that is supported in your installation. Recent installs can likely use Qt5Agg, whereas older ones might need to try Qt4Agg.


Further Note: Analyser won't be able to let you supervise peak selection with Spyder when using Qt5Agg in Python 3. It will open a window displaying MEPs and, depending on how long the kernel stays running, will allow you to interact with it for a short while before freezing. The reasons for this are to do with Spyder/iPython and there is no satisfactory way for me to work around it.
