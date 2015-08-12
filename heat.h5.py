import h5py
import matplotlib
matplotlib.use("GTKAgg")	# Change this as desired.
import gobject
from pylab import *

# if len(sys.argv) != 3:
#   print( "Error: invalid arguments" )
#   print( "usage: heat <h5filename> <Value of z-plane slab>" )
#   exit()
  
# # Obtain filename from command-line parameters
# filename = sys.argv[1]
# section = sys.argv[2]

#### Value of z-plane slab
section = 50  	

# Open file
#### file = h5py.File( filename, "r" )
file = h5py.File( "data.h5", "r" )

# Extract temperature data
temperature = file["temperature"]

# Function called for updating the figure
def updatefig(*args):
  global temperature, frame
  im.set_array(temperature[frame,:,:,section])
  manager.canvas.draw()
  frame+=1
  print "Rendering timestep t=",frame
  if(frame>=len(temperature)):
    return False
  return True

# Create figure and plot initial data
fig = plt.figure(1)
img = subplot(111)
im  = img.imshow( temperature[0,:,:,section], cmap=cm.hot, interpolation="nearest", origin="lower", vmax=1.01 )
manager = get_current_fig_manager()

frame = 1
fig.colorbar(im)

# Whenever idle, update the figure while updatefig returns True
gobject.idle_add( updatefig )
show()
