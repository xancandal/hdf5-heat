import Scientific.IO.NetCDF as nc
import matplotlib
matplotlib.use('GTKAgg') # Change this as desired.
import gobject
from pylab import *

# if len(sys.argv) != 2:
#   print "Error: invalid arguments"
#   print "Usage: heat <ncfilename>"
#   exit()

# # Obtain filename from command-line parameters
# filename = sys.argv[1]

# Open file
#### ncfile = nc.NetCDFFile( filename, "r" )
ncfile = nc.NetCDFFile( "data.nc", "r" )

# Extract temperature data
temperature = ncfile.variables["temperature"]

frame = 0

# Function called for updating the figure
def updatefig(*args):
  global temperature, frame
  im.set_array(temperature[frame,:,:])
  manager.canvas.draw()
  frame+=1
  print "Rendering timestep t=",frame
  if(frame>=len(temperature)):
    return False
  return True
  
# Create figure and plot initial data
fig = plt.figure(1)
img = subplot(111)
im  = img.imshow( temperature[1,:,:], cmap=cm.hot, interpolation="nearest", origin="lower", vmax=1.01 )
manager = get_current_fig_manager()

frame=1
fig.colorbar(im)

# Whenever idle, update the figure while updatefig returns True
gobject.idle_add(updatefig)
show()