import h5py
import matplotlib
matplotlib.use('GTKAgg') # Change this as desired.
import gobject
from pylab import *

### from enthought.mayavi import mlab # Uncomment for Linux
### from enthought.mayavi.sources.array_source import ArraySource # Uncomment for Linux
from mayavi import mlab  # For OSX

# if len(sys.argv) != 4:
#   print( "Error: invalid arguments" )
#   print( "usage: heat <h5filename> <Value of z-plane slab> <time>" )
#   exit()
  
#### Obtain filename from command-line parameters
# filename = sys.argv[1]
# section = sys.argv[2]
# time = sys.argv[3]

#### Value of z-plane slab
section = 50    

#### Time step
time = 50

# Open file
#### file = h5py.File( filename, "r" )
file = h5py.File( "data.h5", "r" )

# Extract temperature data
temperature = file["temperature"]

# Function called for updating the figure
def updatefig(*args):
  global temperature, frame
  frame+=1
  print "Rendering timestep t=",frame
  if(frame>=len(temperature)):
    return False
  return True

frame = 0
while( updatefig(True) ):
  if( frame % 10 == 0 ):
    fig = mlab.contour3d( temperature[time,:,:,:], colormap='hot', opacity=0.3 )
fig = mlab.contour3d( temperature[time,:,:,:], colormap='hot', opacity=0.3 )
mlab.show()

# Create figure and plot initial data
fig = plt.figure(1)
img = subplot(111)
im  = img.imshow( temperature[time,:,:,section], cmap=cm.hot, interpolation="nearest", origin="lower", vmax=1.01 )
manager = get_current_fig_manager()

frame = 1
fig.colorbar(im)

# Whenever idle, update the figure while updatefig returns True
gobject.idle_add( updatefig )
show()







