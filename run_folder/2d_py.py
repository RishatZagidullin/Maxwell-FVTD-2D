from PyQt5 import QtGui  # (the example applies equally well to PySide)
import numpy as np
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui
from PyQt5.QtCore import Qt
import pyqtgraph.ptime as ptime
import pyqtgraph.opengl as gl
from pyqtgraph.pgcollections import OrderedDict

Gradients = OrderedDict([
	('bw', {'ticks': [(0.0, (0, 0, 0, 255)), (1, (255, 255, 255, 255))], 'mode': 'rgb'}),
    ('rgb', {'ticks': [(0.0, (0, 0, 255, 255)), (0.5, (0, 255, 0, 255)), (1.0, (255, 0, 0, 255))], 'mode': 'rgb'}),
    ('hot', {'ticks': [(0.3333, (185, 0, 0, 255)), (0.6666, (255, 220, 0, 255)), (1, (255, 255, 255, 255)), (0, (0, 0, 0, 255))], 'mode': 'rgb'}),
    ('jet', {'ticks': [(1, (166, 0, 0, 255)), (0.32247191011235954, (0, 255, 255, 255)), (0.11348314606741573, (0, 68, 255, 255)), (0.6797752808988764, (255, 255, 0, 255)), (0.902247191011236, (255, 0, 0, 255)), (0.0, (0, 0, 166, 255)), (0.5022471910112359, (0, 255, 0, 255))], 'mode': 'rgb'}),
    ('summer', {'ticks': [(1, (255, 255, 0, 255)), (0.0, (0, 170, 127, 255))], 'mode': 'rgb'} ),
    ('space', {'ticks': [(0.562, (75, 215, 227, 255)), (0.087, (255, 170, 0, 254)), (0.332, (0, 255, 0, 255)), (0.77, (85, 0, 255, 255)), (0.0, (255, 0, 0, 255)), (1.0, (255, 0, 127, 255))], 'mode': 'rgb'}),
    ('winter', {'ticks': [(1, (0, 255, 127, 255)), (0.0, (0, 0, 255, 255))], 'mode': 'rgb'})
])
animation_on = False
def stop_visual():
    global animation_on
    animation_on = False

def start_visual():
    global animation_on
    animation_on = True

def reset_visual():
   global i, animation_on
   i = 0
   animation_on = False
   iterationvalue.setText('T= {}'.format(i))
   val1 = slider.value()*4
   colors  = np.genfromtxt("colours.txt", skip_header=i*len(faces), max_rows=len(faces), usecols = (val1, val1+1, val1+2, val1+3))
   m1.setMeshData(vertexes=verts, faces=faces, faceColors=colors, smooth=False, drawEdges=True, edgeColor=(1,1,1,1)) 

def change_particle_size():
    global slidervalue, i, slider, colors, verts, faces
    slidervalue.setText('Particle size = {:.1f}'.format(slider.value()/10.0+0.1))
    val1 = slider.value()*4
    colors  = np.genfromtxt("colours.txt", skip_header=i*len(faces), max_rows=len(faces), usecols = (val1, val1+1, val1+2, val1+3))
    m1.setMeshData(vertexes=verts, faces=faces, faceColors=colors, smooth=False, drawEdges=True, edgeColor=(1,1,1,1))   

## Always start by initializing Qt (only once per application)
app = QtGui.QApplication([])

## Define a top-level widget to hold everything
w = QtGui.QWidget()
w.setWindowTitle("advection-coagulation")
## Create some widgets to be placed inside
start_btn = QtGui.QPushButton('start')
start_btn.clicked.connect(start_visual)
stop_btn = QtGui.QPushButton('stop')
stop_btn.clicked.connect(stop_visual)
reset_btn = QtGui.QPushButton('reset')
reset_btn.clicked.connect(reset_visual)
## Create a grid layout to manage the widgets size and position
layout = QtGui.QGridLayout()
w.setLayout(layout)
w.resize(875, 800)
w.show()
## Add widgets to the layout in their proper positions
i = 0
iterationvalue = QtGui.QLabel('T= {}'.format(i))
layout.addWidget(iterationvalue, 1, 1)
layout.addWidget(start_btn, 2, 1)
layout.addWidget(stop_btn, 3, 1)
layout.addWidget(reset_btn, 4, 1)

slider = QtGui.QSlider(Qt.Horizontal)
slider.setTickPosition(QtGui.QSlider.TicksBelow)
slider.setSingleStep(1)
slider.setRange(0, 9)
slider.setValue(0)
slider.valueChanged.connect(change_particle_size)

slidervalue = QtGui.QLabel('Particle size = {:.1f}'.format(slider.value()/10.0+0.1))

layout.addWidget(slider, 2, 0)
layout.addWidget(slidervalue, 1, 0)
grid = gl.GLViewWidget()
grid.setCameraPosition(elevation=90, azimuth=0) ##, distance=0)

verts = np.loadtxt("vertices.txt")
faces = np.loadtxt("faces.txt", dtype=int)
colors = np.genfromtxt("colours.txt", max_rows=len(faces), usecols = (0, 1, 2, 3))
## Mesh item will automatically compute face normals.
m1 = gl.GLMeshItem(vertexes=verts, faces=faces, faceColors=colors, smooth=False, drawEdges=True, edgeColor=(1,1,1,1))
m1.translate(-0.75, -0.25, 0.0)
m1.setGLOptions('additive')
m1.rotate(90, 0,0,1)
grid.addItem(m1)

gw = pg.GradientEditorItem(orientation='right')
# load predefined color gradient
GradiendMode = Gradients["rgb"]
gw.restoreState(GradiendMode)

ax = pg.AxisItem('left')
ax.setRange(-1.0, 1.0)
cb = pg.GraphicsLayoutWidget()

cb.addItem(ax)
cb.addItem(gw)
cb.resize(100, 700)
cb.show()
layout.addWidget(cb, 0, 1)
layout.addWidget(grid, 0, 0)
layout.setColumnStretch(0, 1)

def updateData():
    global colors, i, animation_on, slider, verts, faces, total_frames

    QtCore.QTimer.singleShot(100, updateData)

    if (animation_on):
        iterationvalue.setText('T= {}'.format(i))
        val1 = slider.value()*4
        colors  = np.genfromtxt("colours.txt", skip_header=i*len(faces), max_rows=len(faces), usecols = (val1, val1+1, val1+2, val1+3))
        i = (i+1)%total_frames
        m1.setMeshData(vertexes=verts, faces=faces, faceColors=colors, smooth=False, drawEdges=True, edgeColor=(1,1,1,1))   
        grid.grabFrameBuffer().save('fileName_%d.png' % i)
    

updateData()

## Start Qt event loop unless running in interactive mode.

if __name__ == '__main__':
    import sys
    total_frames = int(sys.argv[1])
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()
