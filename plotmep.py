from tkMessageBox import askokcancel
import warnings
warnings.simplefilter(action = "ignore", category = FutureWarning)
import matplotlib
matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy.spatial as spatial
from copy import deepcopy

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
    colours = {'user':'red','original':'green','fill':'cyan','good_plot':'black','bad_plot':'red'}

    def __init__(self,channel, header, query_type):
        self.channel_name = channel['header']['title']
        self.query_type = query_type
        self.current_trial = 0
        self.current_plot = None
        self.fig = plt.figure()
        self.fig.canvas.set_window_title(header['file'])      
        self.axes = self.fig.add_subplot(1,1,1)
        self.fig.canvas.draw()
        self.points = {}
        if query_type == 'time_window':
            self.peak1_offset = queryManager.text_offsets['left']
            self.peak2_offset = queryManager.text_offsets['right']
            self.peak1_arrowlength = queryManager.arrow_lengths['left']
            self.peak2_arrowlength = queryManager.arrow_lengths['right']
            self.peak1_label = 'Onset'
            self.peak2_label = 'Offset'
            self.fill = None
            self.data = channel['rect']
        elif query_type == 'ptp':
            self.peak1_offset = queryManager.text_offsets['top']
            self.peak2_offset = queryManager.text_offsets['bottom']
            self.peak1_arrowlength = queryManager.arrow_lengths['top']
            self.peak2_arrowlength = queryManager.arrow_lengths['bottom']
            self.peak1_label = 'P1'
            self.peak2_label = 'P2'
            self.fill = None
            self.data = channel['data']
        self.rejected = channel['rejected']
        self.points['original'] = channel[query_type]
        self.points['user'] = deepcopy(channel[query_type])
        self.axes.set_ylabel(r'$\mu$V')
        self.axes.set_xlabel('time (ms)')
        self.axes.set_xticks(range(0,len(self.data[0]) + 1,100))
        self.axes.set_xticklabels([str((x * header['sample_rate']) - header['pretrigger_time']) for x in range(0,len(self.data[0]) + 1,100)])
        self.axes.set_xlim([0, len(self.data[0])])
        self.axes.format_coord = lambda x, y: ''
        self.fig.tight_layout()
        self.fig.canvas.manager.window.showMaximized()
        self.background = self.fig.canvas.copy_from_bbox(self.axes.bbox)
        self.fig.canvas.mpl_connect('key_release_event', self.key_release)
        self.fig.canvas.mpl_connect('button_press_event', self.mouse_press)
        self.fig.canvas.mpl_connect('scroll_event', self.mouse_scroll)
        self.current_plot, = self.axes.plot(self.data[0],color=queryManager.colours['good_plot'])
        self.original_markers = {}
        self.original_markers['peak1'] = self.axes.annotate(self.peak1_label,color='white',xy=(0,0),xytext=self.peak1_offset,textcoords="offset points",arrowprops={'facecolor':queryManager.colours['original'],'shrink':0.05,'headlength':self.peak1_arrowlength}, horizontalalignment='center',visible=False)
        self.original_markers['peak2'] = self.axes.annotate(self.peak2_label,color='white',xy=(0,0),xytext=self.peak2_offset,textcoords="offset points",arrowprops={'facecolor':queryManager.colours['original'],'shrink':0.05,'headlength':self.peak2_arrowlength}, horizontalalignment='center',visible=False)
        self.user_markers = {}
        self.user_markers['peak1'] = self.axes.annotate(self.peak1_label,xy=(0,0),xytext=self.peak1_offset,textcoords="offset points",arrowprops={'facecolor':queryManager.colours['user'],'shrink':0.05,'headlength':self.peak1_arrowlength}, horizontalalignment='center',visible=False)
        self.user_markers['peak2'] = self.axes.annotate(self.peak2_label,xy=(0,0),xytext=self.peak2_offset,textcoords="offset points",arrowprops={'facecolor':queryManager.colours['user'],'shrink':0.05,'headlength':self.peak2_arrowlength}, horizontalalignment='center',visible=False)
        self.cursor = dataCursor(self.axes, range(len(self.data[0])), self.data[0],header['sample_rate'], -header['pretrigger_time'])
        self.update_display()
        plt.show()

    def update_points(self):
        for idx in xrange(len(self.points['original'])):
            self.points['original'][idx] = self.points['user'][idx]

    def update_display(self):
        self.fig.suptitle(self.channel_name + ': Trial ' + str(self.current_trial + 1))
        self.current_plot.set_ydata(self.data[self.current_trial])
        if self.rejected[self.current_trial]:
            self.current_plot.set_c(queryManager.colours['bad_plot'])
        else:
            self.current_plot.set_c(queryManager.colours['good_plot'])
        self.cursor.dot.visible=False
        self.cursor.annotation.visible=False
        ylim = max(250,((max(map(abs,self.data[self.current_trial])) // 100) + 2) * 100)
        if self.query_type == 'time_window':
            self.axes.set_ylim([0,ylim])
        elif self.query_type == 'ptp':
            self.axes.set_ylim([-ylim,ylim])
        if self.fill:
            self.fill.remove()
            del self.fill
            self.fill = None
        if self.points['original'][self.current_trial][1]:
            self.original_markers['peak1'].xy = (self.points['original'][self.current_trial][1],self.data[self.current_trial][self.points['original'][self.current_trial][1]])
            self.original_markers['peak2'].xy = (self.points['original'][self.current_trial][2],self.data[self.current_trial][self.points['original'][self.current_trial][2]])
            self.original_markers['peak1'].set_visible(True)
            self.original_markers['peak2'].set_visible(True)
        else:
            self.original_markers['peak1'].set_visible(False)
            self.original_markers['peak2'].set_visible(False)
        if self.points['user'][self.current_trial][1]:
            self.user_markers['peak1'].xy = (self.points['user'][self.current_trial][1],self.data[self.current_trial][self.points['user'][self.current_trial][1]])
            self.user_markers['peak1'].set_visible(True)
        else:
            self.user_markers['peak1'].set_visible(False)
        if self.points['user'][self.current_trial][2]:
            self.user_markers['peak2'].xy = (self.points['user'][self.current_trial][2],self.data[self.current_trial][self.points['user'][self.current_trial][2]])
            self.user_markers['peak2'].set_visible(True)
            if (self.query_type == 'time_window') and self.points['user'][self.current_trial][1]:
                self.fill = self.axes.fill_between(range(*self.points['user'][self.current_trial][1:]),self.data[self.current_trial][self.points['user'][self.current_trial][1]:self.points['user'][self.current_trial][2]],facecolor=queryManager.colours['fill'])
                self.axes.draw_artist(self.fill)
        else:
            self.user_markers['peak2'].set_visible(False)
        self.axes.draw_artist(self.original_markers['peak1'])
        self.axes.draw_artist(self.original_markers['peak2'])
        self.axes.draw_artist(self.user_markers['peak1'])
        self.axes.draw_artist(self.user_markers['peak2'])
        self.axes.draw_artist(self.current_plot)

    def mouse_scroll(self,event):
            new_trial = max(0,self.current_trial - 1) if event.step > 0 else min(len(self.data) - 1,self.current_trial + 1)
            if new_trial != self.current_trial:
                self.fig.canvas.restore_region(self.background)
                self.current_trial = new_trial
                self.update_display()
                self.fig.canvas.draw()
                self.cursor.background = self.fig.canvas.copy_from_bbox(self.axes.bbox)
                self.cursor.update(self.data[new_trial])

    def key_release(self,event):
        if event.key == 'enter':
            if askokcancel('Accept?','Accept all trials?'):
                self.update_points()
                plt.close()
        elif event.key == 'backspace':
            if self.rejected[self.current_trial]:
                self.current_plot.set_c(queryManager.colours['good_plot'])
                self.rejected[self.current_trial] = False
            else:
                self.current_plot.set_c(queryManager.colours['bad_plot'])
                self.rejected[self.current_trial] = True
            self.fig.canvas.draw()
            self.cursor.background = self.fig.canvas.copy_from_bbox(self.axes.bbox)

        elif event.key == 'delete':
            self.user_markers['peak1'].set_visible(False)
            self.user_markers['peak2'].set_visible(False)
            self.axes.draw_artist(self.user_markers['peak1'])
            self.axes.draw_artist(self.user_markers['peak2'])
            self.points['user'][self.current_trial] = [None,None]
            if (self.query_type == 'time_window') and self.fill:
                self.fill.remove()
                del self.fill
                self.fill = None
            self.fig.canvas.draw()
            self.cursor.background = self.fig.canvas.copy_from_bbox(self.axes.bbox)
        elif event.key == 'escape':
            if askokcancel('Quit?','Quit MEP analysis?'):
                exit()
        elif event.key in {'left','right'}:
            new_trial = max(0,self.current_trial - 1) if event.key == 'left' else min(len(self.data) - 1,self.current_trial + 1)
            if new_trial != self.current_trial:
                self.current_trial = new_trial
                self.update_display()
                self.fig.canvas.draw()
                self.cursor.background = self.fig.canvas.copy_from_bbox(self.axes.bbox)
                self.cursor.update(self.data[new_trial])

    def mouse_press(self,event):
        if (event.button == 1) and event.xdata:
            self.points['user'][self.current_trial][1] = int(round(self.cursor.annotation.xy[0]))
            self.user_markers['peak1'].xy = (self.points['user'][self.current_trial][1],self.data[self.current_trial][self.points['user'][self.current_trial][1]])
            self.user_markers['peak1'].set_visible(True)
            self.axes.draw_artist(self.user_markers['peak1'])
            if self.points['user'][self.current_trial][1] and (self.query_type == 'time_window'):
                if self.fill:
                    self.fill.remove()
                    del self.fill
                    self.fill = None
                self.fill = self.axes.fill_between(range(*self.points['user'][self.current_trial][1:]),self.data[self.current_trial][self.points['user'][self.current_trial][1]:self.points['user'][self.current_trial][2]],facecolor=queryManager.colours['fill'])
            self.fig.canvas.draw()
            self.cursor.background = self.fig.canvas.copy_from_bbox(self.axes.bbox)
        elif (event.button == 3) and event.xdata:
            self.points['user'][self.current_trial][2] = int(round(self.cursor.annotation.xy[0]))
            self.user_markers['peak2'].xy = (self.points['user'][self.current_trial][2],self.data[self.current_trial][self.points['user'][self.current_trial][2]])
            self.user_markers['peak2'].set_visible(True)
            self.axes.draw_artist(self.user_markers['peak2'])
            if self.points['user'][self.current_trial][1] and (self.query_type == 'time_window'):
                if self.fill:
                    self.fill.remove()
                    del self.fill
                    self.fill = None
                self.fill = self.axes.fill_between(range(*self.points['user'][self.current_trial][1:]),self.data[self.current_trial][self.points['user'][self.current_trial][1]:self.points['user'][self.current_trial][2]],facecolor=queryManager.colours['fill'])
            self.fig.canvas.draw()
            self.cursor.background = self.fig.canvas.copy_from_bbox(self.axes.bbox)

def plot_data(mep_dataset,query_type):
    # Cycle through channels
    for channel in mep_dataset.channels:
        if query_type is None:
            query_type = ['ptp'] if channel['ptp'] else []
            query_type += ['time_window'] if channel['time_window'] else []
        # Pass all arguments to queryManager
        if 'ptp' in query_type:
            if channel['ptp']:
                queryManager(channel,mep_dataset.header,'ptp')
            else:
                raise ValueError('No peak-to-peak data found in MEPDataset.')
        if 'time_window' in query_type:
            if channel['time_window']:
                queryManager(channel,mep_dataset.header,'time_window')
            else:
                raise ValueError('No peak-to-peak data found in MEPDataset.')