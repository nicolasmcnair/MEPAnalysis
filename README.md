# MEPAnalysis
Python script for analysing LabChart binary files containing MEP data.

Figure Controls:
Left click: Set P1/Onset
Right click: Set P2/Offset
Enter: Accept all trials (N.B. closing figure window will reject all changes, but continue with original points)
Delete: Remove peaks from current trial
Backspace: Toggle flag for background movement detection
Left/Right arrow: Cycle backward/forward through trials
Scroll Up/Down: Cycle backward/forward through trials
Escape: Quit entire MEP data analysis (this will mean no file is written)

Note: Depending on the Python version/distro you are running, you may need to alter line 8 in plotmep.py to a graphics backend that is supported in your installation. Recent installs can likely use Qt5Agg, whereas older ones might need to try Qt4Agg.

Further Note: Analyser won't be able to let you supervise peak selection with Spyder when using Qt5Agg in Python 3. It will open a window displaying MEPs and, depending on how long the kernel stays running, will allow you to interact with it for a short while before freezing. The reasons for this are to do with Spyder/iPython and there is no satisfactory way for me to work around it.
