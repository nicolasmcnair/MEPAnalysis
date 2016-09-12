from tkMessageBox import askokcancel
import warnings
warnings.simplefilter(action = "ignore", category = FutureWarning)
import matplotlib
matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy.spatial as spatial
from copy import copy,deepcopy

class dataCursor(object):
    """Display the x,y location of the nearest data point."""
    def __init__(self, ax, x, y, rescale, shift, offsets=(-20, 20)):
        self.x = np.asarray(x, dtype='float')
        y = np.asarray(y, dtype='float')
        self._points = np.column_stack((x, y))
        self.rescale = rescale
        self.shift = shift
        self.tree = spatial.cKDTree(self._points)
        self.ax = ax
        self.fig = ax.figure
        self.dot = ax.scatter([0], [y.min()], s=130, color='green', alpha=0.7,animated=True)
        self.annotation = self.ax.annotate('', xy=(0, 0), ha = 'right', xytext = offsets, textcoords = 'offset points', va = 'bottom', bbox = dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.75,animated=True), arrowprops = dict(arrowstyle='->', connectionstyle='arc3,rad=0',animated=True),animated=True)
        self.background = None
        plt.connect('motion_notify_event', self)

    def __call__(self, event):
        if event.inaxes:
            if not self.background:
                self.background = self.fig.canvas.copy_from_bbox(self.ax.bbox)
            x, y = self._points[self.tree.query((event.xdata, event.ydata), k=1, p=2)[1]]
            event.canvas.restore_region(self.background)
            self.dot.visible=True
            self.annotation.visible=True
            self.annotation.xy = x, y
            self.xytext = (20,20)
            self.annotation.set_text('x: {x:0.2f}\ny: {y:0.2f}'.format(x=(x * self.rescale) + self.shift, y=y))
            self.dot.set_offsets((x, y))
            self.ax.draw_artist(self.annotation)
            self.ax.draw_artist(self.dot)
            event.canvas.blit(self.ax.bbox)

    def update(self,y):
        y = np.asarray(y, dtype='float')
        self._points = np.column_stack((self.x, y))
        self.tree = spatial.cKDTree(self._points)

class queryManager(object):

    text_offsets = {'left':(-30,15),'right':(30,15),'top':(0,18),'bottom':(0,-27)}
    arrow_lengths = {'left':15,'right':15,'top':15,'bottom':15}
    colours = {'user':'red','original':'green','fill':'cyan','plot':'black'}

    def __init__(self,filename,channelName,channelData,mepInfo,tickrate,pretriggerTime,windowType):
        self.channelName = channelName
        self.windowType = windowType
        self.tickrate = tickrate
        self.pretriggerTime = pretriggerTime    
        self.currentTrial = 0
        self.currentPlot = None
        self.fig = plt.figure()
        self.fig.canvas.set_window_title(filename)      
        self.axes = self.fig.add_subplot(1,1,1)
        self.fig.canvas.draw()
        self.points = {}
        self.pMarkers = {}
        if windowType:
            self.p1Offset = queryManager.text_offsets['left']
            self.p2Offset = queryManager.text_offsets['right']
            self.p1ArrowLength = queryManager.arrow_lengths['left']
            self.p2ArrowLength = queryManager.arrow_lengths['right']
            self.p1Label = 'Onset'
            self.p2Label = 'Offset'
            self.fill = None
        else:
            self.p1Offset = queryManager.text_offsets['top']
            self.p2Offset = queryManager.text_offsets['bottom']
            self.p1ArrowLength = queryManager.arrow_lengths['top']
            self.p2ArrowLength = queryManager.arrow_lengths['bottom']
            self.p1Label = 'P1'
            self.p2Label = 'P2'
            self.fill = None
        self.data = [copy(channelData[i]) for i in range(len(channelData))]
        self.points['original'] = [(mepInfo[i][1],mepInfo[i][2]) for i in range(len(mepInfo))]
        self.points['user'] = [[mepInfo[i][1],mepInfo[i][2]] for i in range(len(mepInfo))]
        self.axes.set_ylabel(r'$\mu$V')
        self.axes.set_xlabel('time (ms)')
        self.axes.set_xticks(range(0,len(self.data[0]) + 1,100))
        self.axes.set_xticklabels([str((x * tickrate) - pretriggerTime) for x in range(0,len(self.data[0]) + 1,100)])
        self.axes.set_xlim([0, len(self.data[0])])
        self.axes.format_coord = lambda x, y: ''
        self.fig.tight_layout()
        self.fig.canvas.manager.window.showMaximized()
        self.background = self.fig.canvas.copy_from_bbox(self.axes.bbox)
        self.fig.canvas.mpl_connect('key_release_event', self.keyRelease)
        self.fig.canvas.mpl_connect('button_press_event', self.mousePress)
        self.fig.canvas.mpl_connect('scroll_event', self.mouseScroll)
        self.currentPlot, = self.axes.plot(self.data[0],color=queryManager.colours['plot'])
        self.originalMarkers = {}
        self.originalMarkers['p1'] = self.axes.annotate(self.p1Label,color='white',xy=(0,0),xytext=self.p1Offset,textcoords="offset points",arrowprops={'facecolor':queryManager.colours['original'],'shrink':0.05,'headlength':self.p1ArrowLength}, horizontalalignment='center',visible=False)
        self.originalMarkers['p2'] = self.axes.annotate(self.p2Label,color='white',xy=(0,0),xytext=self.p2Offset,textcoords="offset points",arrowprops={'facecolor':queryManager.colours['original'],'shrink':0.05,'headlength':self.p2ArrowLength}, horizontalalignment='center',visible=False)
        self.userMarkers = {}
        self.userMarkers['p1'] = self.axes.annotate(self.p1Label,xy=(0,0),xytext=self.p1Offset,textcoords="offset points",arrowprops={'facecolor':queryManager.colours['user'],'shrink':0.05,'headlength':self.p1ArrowLength}, horizontalalignment='center',visible=False)
        self.userMarkers['p2'] = self.axes.annotate(self.p2Label,xy=(0,0),xytext=self.p2Offset,textcoords="offset points",arrowprops={'facecolor':queryManager.colours['user'],'shrink':0.05,'headlength':self.p2ArrowLength}, horizontalalignment='center',visible=False)
        self.cursor = dataCursor(self.axes, range(len(self.data[0])), self.data[0],tickrate, -pretriggerTime)
        self.updateDisplay()
        plt.show()

    def updateDisplay(self):
        self.fig.suptitle(self.channelName + ': Trial ' + str(self.currentTrial + 1))
        self.currentPlot.set_ydata(self.data[self.currentTrial])
        self.cursor.dot.visible=False
        self.cursor.annotation.visible=False
        ylim = max(250,((max(map(abs,self.data[self.currentTrial])) // 100) + 2) * 100)
        if self.windowType:
            self.axes.set_ylim([0,ylim])
        else:
            self.axes.set_ylim([-ylim,ylim])
        if self.fill:
            self.fill.remove()
            del self.fill
            self.fill = None
        if self.points['original'][self.currentTrial][0]:
            self.originalMarkers['p1'].xy = (self.points['original'][self.currentTrial][0],self.data[self.currentTrial][self.points['original'][self.currentTrial][0]])
            self.originalMarkers['p2'].xy = (self.points['original'][self.currentTrial][1],self.data[self.currentTrial][self.points['original'][self.currentTrial][1]])
            self.originalMarkers['p1'].set_visible(True)
            self.originalMarkers['p2'].set_visible(True)
        else:
            self.originalMarkers['p1'].set_visible(False)
            self.originalMarkers['p2'].set_visible(False)
        if self.points['user'][self.currentTrial][0]:
            self.userMarkers['p1'].xy = (self.points['user'][self.currentTrial][0],self.data[self.currentTrial][self.points['user'][self.currentTrial][0]])
            self.userMarkers['p1'].set_visible(True)
        else:
            self.userMarkers['p1'].set_visible(False)
        if self.points['user'][self.currentTrial][1]:
            self.userMarkers['p2'].xy = (self.points['user'][self.currentTrial][1],self.data[self.currentTrial][self.points['user'][self.currentTrial][1]])
            self.userMarkers['p2'].set_visible(True)
            if self.windowType and self.points['user'][self.currentTrial][0]:
                self.fill = self.axes.fill_between(range(*self.points['user'][self.currentTrial]),self.data[self.currentTrial][self.points['user'][self.currentTrial][0]:self.points['user'][self.currentTrial][1]],facecolor=queryManager.colours['fill'])
                self.axes.draw_artist(self.fill)
        else:
            self.userMarkers['p2'].set_visible(False)
        self.axes.draw_artist(self.originalMarkers['p1'])
        self.axes.draw_artist(self.originalMarkers['p2'])
        self.axes.draw_artist(self.userMarkers['p1'])
        self.axes.draw_artist(self.userMarkers['p2'])
        self.axes.draw_artist(self.currentPlot)

    def mouseScroll(self,event):
            newTrial = max(0,self.currentTrial - 1) if event.step > 0 else min(len(self.data) - 1,self.currentTrial + 1)
            if newTrial != self.currentTrial:
                self.fig.canvas.restore_region(self.background)
                self.currentTrial = newTrial
                self.updateDisplay()
                self.fig.canvas.draw()
                self.cursor.background = self.fig.canvas.copy_from_bbox(self.axes.bbox)
                self.cursor.update(self.data[newTrial])

    def keyRelease(self,event):
        if event.key == 'enter':
            if askokcancel('Accept?','Accept all trials?'):
                plt.close()
        elif event.key in {'delete','backspace'}:
            self.userMarkers['p1'].set_visible(False)
            self.userMarkers['p2'].set_visible(False)
            self.axes.draw_artist(self.userMarkers['p1'])
            self.axes.draw_artist(self.userMarkers['p2'])
            self.points['user'][self.currentTrial] = [None,None]
            if self.windowType and self.fill:
                self.fill.remove()
                del self.fill
                self.fill = None
            self.fig.canvas.draw()
            self.cursor.background = self.fig.canvas.copy_from_bbox(self.axes.bbox)
        elif event.key == 'escape':
            if askokcancel('Quit?','Quit analysis?'):
                exit()
        elif event.key in {'left','right'}:
            newTrial = max(0,self.currentTrial - 1) if event.key == 'left' else min(len(self.data) - 1,self.currentTrial + 1)
            if newTrial != self.currentTrial:
                self.currentTrial = newTrial
                self.updateDisplay()
                self.fig.canvas.draw()
                self.cursor.background = self.fig.canvas.copy_from_bbox(self.axes.bbox)
                self.cursor.update(self.data[newTrial])

    def mousePress(self,event):
        if (event.button == 1) and event.xdata:
            self.points['user'][self.currentTrial][0] = int(round(self.cursor.annotation.xy[0]))
            self.userMarkers['p1'].xy = (self.points['user'][self.currentTrial][0],self.data[self.currentTrial][self.points['user'][self.currentTrial][0]])
            self.userMarkers['p1'].set_visible(True)
            self.axes.draw_artist(self.userMarkers['p1'])
            if self.points['user'][self.currentTrial][1] and self.windowType:
                if self.fill:
                    self.fill.remove()
                    del self.fill
                    self.fill = None
                self.fill = self.axes.fill_between(range(*self.points['user'][self.currentTrial]),self.data[self.currentTrial][self.points['user'][self.currentTrial][0]:self.points['user'][self.currentTrial][1]],facecolor=queryManager.colours['fill'])
            self.fig.canvas.draw()
            self.cursor.background = self.fig.canvas.copy_from_bbox(self.axes.bbox)
        elif (event.button == 3) and event.xdata:
            self.points['user'][self.currentTrial][1] = int(round(self.cursor.annotation.xy[0]))
            self.userMarkers['p2'].xy = (self.points['user'][self.currentTrial][1],self.data[self.currentTrial][self.points['user'][self.currentTrial][1]])
            self.userMarkers['p2'].set_visible(True)
            self.axes.draw_artist(self.userMarkers['p2'])
            if self.points['user'][self.currentTrial][0] and self.windowType:
                if self.fill:
                    self.fill.remove()
                    del self.fill
                    self.fill = None
                self.fill = self.axes.fill_between(range(*self.points['user'][self.currentTrial]),self.data[self.currentTrial][self.points['user'][self.currentTrial][0]:self.points['user'][self.currentTrial][1]],facecolor=queryManager.colours['fill'])
            self.fig.canvas.draw()
            self.cursor.background = self.fig.canvas.copy_from_bbox(self.axes.bbox)

def queryData(filename,channelName,channelData,mepInfo,tickrate,pretriggerTime,windowType = None):
    # Pass all arguments to queryManager
    myQuery = queryManager(**locals())
    if windowType:
        for j in range(len(channelData)):
            onset,offset = myQuery.points['user'][j]
            if onset and offset:
                windowAmplitude = np.trapz(myQuery.data[j][onset:offset + 1],dx=1)
                if windowType == 2:
                    windowAmplitude /= (offset - onset) + 1
            else:
                windowAmplitude = onset = offset = None
            mepInfo[j] = (windowAmplitude,onset,offset)
    else:
        for j in range(len(channelData)):
            p1,p2 = myQuery.points['user'][j]
            if p1 and p2:
                ptp = myQuery.data[j][p1] - myQuery.data[j][p2]
            else:
                ptp = p1 = p2 = None
            mepInfo[j] = (ptp,p1,p2)