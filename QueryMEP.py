from tkMessageBox import askokcancel
import warnings
warnings.simplefilter(action = "ignore", category = FutureWarning)
import matplotlib
matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy.spatial as spatial

class dataCursor(object):
    """Display the x,y location of the nearest data point."""
    def __init__(self, ax, x, y, rescale, shift, offsets=(-20, 20)):
        x = np.asarray(x, dtype='float')
        y = np.asarray(y, dtype='float')
        self._points = np.column_stack((x, y))
        self.rescale = rescale
        self.shift = shift
        self.tree = spatial.cKDTree(self._points)
        self.ax = ax
        self.fig = ax.figure
        self.dot = ax.scatter([x.min()], [y.min()], s=130, color='green', alpha=0.7)
        self.annotation = self.ax.annotate('', xy=(0, 0), ha = 'right', xytext = offsets, textcoords = 'offset points', va = 'bottom', bbox = dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.75), arrowprops = dict(arrowstyle='->', connectionstyle='arc3,rad=0'))
        plt.connect('motion_notify_event', self)

    def __call__(self, event):
        if event.inaxes:
            x, y = self._points[self.tree.query((event.xdata, event.ydata), k=1, p=2)[1]]
            self.annotation.xy = x, y
            self.annotation.set_text('x: {x:0.2f}\ny: {y:0.2f}'.format(x=(x * self.rescale) + self.shift, y=y))
            self.dot.set_offsets((x, y))
            event.canvas.draw()

def queryWindow(channelName,channelData,tickrate,pretrigger_time,windowType):
    global queryManager
    queryManager = {'currentTrial':None,'rectData':None,'onset':None,'onsetMarker':None,'offset':None,'offsetMarker':None,'continueQuery':True,'cursor':None,'fill':None}

    redArrow = {'facecolor':'red','shrink':0.05,'frac':1.0}
    leftOffset = (-30,15)
    rightOffset = (30,15)
    highlightColor = 'cyan'

    def keyRelease(event):
        global queryManager

        if event.key == 'ctrl+enter':
            if askokcancel('Accept All?','Accept all remaining trials?'):
                queryManager['continueQuery'] = False
                plt.close()
        elif event.key == 'enter':
            plt.close()
        elif event.key in {'delete','backspace'}:
            queryManager['onset'] = queryManager['offset'] = None
            plt.close()
        elif event.key == 'escape':
            if askokcancel('Quit?','Quit analysis?'):
                exit()

    def mousePress(event):
        global queryManager

        if (event.button == 1) and event.xdata:
            queryManager['onset'] = int(round(queryManager['cursor'].annotation.xy[0]))
            xyLoc = (queryManager['onset'],queryManager['rectData'][queryManager['onset']])
            if queryManager['onsetMarker']:
                queryManager['onsetMarker'].xy = xyLoc
                queryManager['onsetMarker'].xyann = leftOffset
            else:
                queryManager['onsetMarker'] = event.inaxes.annotate('Onset',xy=xyLoc,xytext=leftOffset,textcoords="offset points",arrowprops=redArrow, horizontalalignment='center')
            if queryManager['offset']:
                if queryManager['fill']:
                    queryManager['fill'].remove()
                queryManager['fill'] = event.inaxes.fill_between(range(queryManager['onset'],queryManager['offset']),queryManager['rectData'][queryManager['onset']:queryManager['offset']],facecolor=highlightColor)
            event.canvas.draw()
        elif (event.button == 3) and event.xdata:
            queryManager['offset'] = int(round(queryManager['cursor'].annotation.xy[0]))
            xyLoc = (queryManager['offset'],queryManager['rectData'][queryManager['offset']])
            if queryManager['offsetMarker']:
                queryManager['offsetMarker'].xy = xyLoc
                queryManager['offsetMarker'].xyann = rightOffset
            else:
                queryManager['offsetMarker'] = event.inaxes.annotate('Offset',xy=xyLoc,xytext=rightOffset,textcoords="offset points",arrowprops=redArrow, horizontalalignment='center')
            if queryManager['onset']:
                if queryManager['fill']:
                    queryManager['fill'].remove()
                queryManager['fill'] = event.inaxes.fill_between(range(queryManager['onset'],queryManager['offset']),queryManager['rectData'][queryManager['onset']:queryManager['offset']],facecolor=highlightColor)
            event.canvas.draw()

    for j in range(len(channelData)): 
        if queryManager['continueQuery']:
            queryManager['currentTrial'] = channelData[j]
            if queryManager['currentTrial']['window'] and queryManager['currentTrial']['window'][0]:
                queryManager['onset'] = queryManager['currentTrial']['window'][1]
                queryManager['offset'] = queryManager['currentTrial']['window'][2]
            else:
                queryManager['onset'] = None
                queryManager['offset'] = None
            fig = plt.figure()
            fig.suptitle(channelName + ': Trial ' + str(j + 1))
            ax = fig.add_subplot(1,1,1)
            queryManager['rectData'] = map(abs,queryManager['currentTrial']['data'])
            plt.plot(queryManager['rectData'],color='black')
            plt.ylim([0,ax.get_ylim()[1]])
            plt.ylabel(r'$\mu$V')
            plt.xlim([0, len(queryManager['currentTrial']['data'])])
            plt.xlabel('time (ms)')
            plt.xticks(range(0,len(queryManager['currentTrial']['data']) + 1,100),[(x * tickrate) - pretrigger_time for x in range(0,len(queryManager['currentTrial']['data']) + 1,100)])
            plt.tight_layout()
            if queryManager['onset']:
                xyLoc = (queryManager['onset'],queryManager['rectData'][queryManager['onset']])
                queryManager['onsetMarker'] = ax.annotate('Onset',xy=xyLoc,xytext=leftOffset,textcoords="offset points",arrowprops=redArrow, horizontalalignment='center')
            elif queryManager['onsetMarker']:
                queryManager['onsetMarker'].remove()
                queryManager['onsetMarker'] = None
            if queryManager['offset']:
                xyLoc = (queryManager['offset'],queryManager['rectData'][queryManager['offset']])
                queryManager['offsetMarker'] = ax.annotate('Offset',xy=xyLoc,xytext=rightOffset,textcoords="offset points",arrowprops=redArrow, horizontalalignment='center')
            elif queryManager['offsetMarker']:
                queryManager['offsetMarker'].remove()
                queryManager['offsetMarker'] = None
            if queryManager['onset']:
                queryManager['fill'] = ax.fill_between(range(queryManager['onset'],queryManager['offset']),queryManager['rectData'][queryManager['onset']:queryManager['offset']],facecolor=highlightColor)
            plt.get_current_fig_manager().window.showMaximized()
            ax.format_coord = lambda x, y: ''
            queryManager['cursor'] = dataCursor(ax, range(len(queryManager['currentTrial']['data'])), queryManager['rectData'],tickrate, -pretrigger_time)
            fig.canvas.mpl_connect('button_press_event', mousePress)
            fig.canvas.mpl_connect('key_release_event', keyRelease)
            plt.show()
            if queryManager['onset'] and queryManager['offset']:
                if windowType == 1:
                    queryManager['currentTrial']['window'] = (np.trapz(queryManager['rectData'][queryManager['onset']:queryManager['offset']],dx=1),queryManager['onset'],queryManager['offset'])
                elif windowType == 2:
                    queryManager['currentTrial']['window'] = (np.mean(queryManager['rectData'][queryManager['onset']:queryManager['offset']]),queryManager['onset'],queryManager['offset'])               
            else:
                queryManager['currentTrial']['window'] = None
            del fig

def queryPeaks(channelName,channelData,tickrate,pretrigger_time):
    global queryManager
    queryManager = {'currentTrial':None,'p1':None,'p1Marker':None,'p2':None,'p2Marker':None,'continueQuery':True,'cursor':None}

    greenArrow = {'facecolor':'green','shrink':0.1,'frac':1.0}
    redArrow = greenArrow.copy()
    redArrow['facecolor'] = 'red'
    upOffset = (0,15)
    downOffset = (0,-20)

    def keyRelease(event):
        global queryManager

        if event.key == 'ctrl+enter':
            if askokcancel('Accept All?','Accept all remaining trials?'):
                queryManager['continueQuery'] = False
                plt.close()
        elif event.key == 'enter':
            plt.close()
        elif event.key in {'delete','backspace'}:
            queryManager['p1'] = queryManager['p2'] = None
            plt.close()
        elif event.key == 'escape':
            if askokcancel('Quit?','Quit analysis?'):
                exit()
        
    def mousePress(event):
        global queryManager

        if (event.button == 1) and event.xdata:
            queryManager['p1'] = int(round(queryManager['cursor'].annotation.xy[0]))
            xyLoc = (queryManager['p1'],queryManager['currentTrial']['data'][queryManager['p1']])
            if queryManager['p1Marker']:
                queryManager['p1Marker'].xy = xyLoc
                queryManager['p1Marker'].xyann = upOffset
            else:
                queryManager['p1Marker'] = event.inaxes.annotate('P1',xy=xyLoc,xytext=upOffset,textcoords="offset points",arrowprops=redArrow, horizontalalignment='center')
            event.canvas.draw()
        elif (event.button == 3) and event.xdata:
            queryManager['p2'] = int(round(queryManager['cursor'].annotation.xy[0]))
            xyLoc = (queryManager['p2'],queryManager['currentTrial']['data'][queryManager['p2']])
            if queryManager['p2Marker']:
                queryManager['p2Marker'].xy = xyLoc
                queryManager['p2Marker'].xyann = downOffset
            else:
                queryManager['p2Marker'] = event.inaxes.annotate('P2',xy=xyLoc,xytext=downOffset,textcoords="offset points",arrowprops=redArrow, horizontalalignment='center')
            event.canvas.draw()
    
    for j in range(len(channelData)): 
        if queryManager['continueQuery']:
            queryManager['currentTrial'] = channelData[j]
            if queryManager['currentTrial']['ptp']:
                queryManager['p1'] = queryManager['currentTrial']['ptp'][1]
                queryManager['p2'] = queryManager['currentTrial']['ptp'][2]
            else:
                queryManager['p1'] = None
                queryManager['p2'] = None
            fig = plt.figure()
            fig.suptitle(channelName + ': Trial ' + str(j + 1))
            ax = fig.add_subplot(1,1,1)
            plt.plot(queryManager['currentTrial']['data'],color='black')
            plt.ylabel(r'$\mu$V')
            plt.xlim([0, len(queryManager['currentTrial']['data'])])
            plt.xlabel('time (ms)')
            plt.xticks(range(0,len(queryManager['currentTrial']['data']) + 1,100),[(x * tickrate) - pretrigger_time for x in range(0,len(queryManager['currentTrial']['data']) + 1,100)])
            plt.tight_layout()
            for p in queryManager['currentTrial']['peaks']:
                xyLoc = (p,queryManager['currentTrial']['data'][p])
                ax.annotate('X',color='white',xy=xyLoc,xytext=upOffset,textcoords="offset points",arrowprops=greenArrow)
            if queryManager['p1']:
                xyLoc = (queryManager['p1'],queryManager['currentTrial']['data'][queryManager['p1']])
                queryManager['p1Marker'] = ax.annotate('P1',xy=xyLoc,xytext=upOffset,textcoords="offset points",arrowprops=redArrow, horizontalalignment='center')
            elif queryManager['p1Marker']:
                queryManager['p1Marker'].remove()
                queryManager['p1Marker'] = None
            for v in queryManager['currentTrial']['valleys']:
                xyLoc = (v,queryManager['currentTrial']['data'][v])
                ax.annotate('X',color='white',xy=xyLoc,xytext=downOffset,textcoords="offset points",arrowprops=greenArrow)
            if queryManager['p2']:
                xyLoc = (queryManager['p2'],queryManager['currentTrial']['data'][queryManager['p2']])
                queryManager['p2Marker'] = ax.annotate('P2',xy=xyLoc,xytext=downOffset,textcoords="offset points",arrowprops=redArrow, horizontalalignment='center')
            elif queryManager['p2Marker']:
                queryManager['p2Marker'].remove()
                queryManager['p2Marker']= None
            plt.get_current_fig_manager().window.showMaximized()
            ax.format_coord = lambda x, y: ''
            queryManager['cursor'] = dataCursor(ax, range(len(queryManager['currentTrial']['data'])), queryManager['currentTrial']['data'],tickrate, -pretrigger_time)
            fig.canvas.mpl_connect('button_press_event', mousePress)
            fig.canvas.mpl_connect('key_release_event', keyRelease)
            plt.show()
            if queryManager['p1'] and queryManager['p2']:
                queryManager['currentTrial']['ptp'] = (queryManager['currentTrial']['data'][queryManager['p1']] - queryManager['currentTrial']['data'][queryManager['p2']],queryManager['p1'],queryManager['p2'])            
                queryManager['currentTrial']['peaks'].add(queryManager['p1'])
                queryManager['currentTrial']['valleys'].add(queryManager['p2'])
            else:
                queryManager['currentTrial']['ptp'] = None
            del fig