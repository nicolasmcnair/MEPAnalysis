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
