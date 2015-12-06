import math, numpy, copy, sys
from optparse import OptionParser
import operator 
import json
import matplotlib.pyplot as plt

parser = OptionParser()
parser.add_option("-e", "--event",type=int, help='Which event to examine', default=100)
options, args = parser.parse_args()

particle_vars = numpy.load("../Data/particle_vars.npy")

etas = []
phis = []
pts = []

for i, particles in enumerate(particle_vars):
  if i!=options.event: continue 
  for particle in particles:
    eta_jet = -2.48893314399
    phi_jet = 4.56855795497
    eta = particle['eta']
    phi = particle['phi']
    etas.append(particle['eta'])
    phis.append(particle['phi'])
    pts.append(particle['pt'])

nbins = 50 
H,xedges,yedges = numpy.histogram2d(etas,phis,bins=nbins,weights=pts)

# H needs to be rotated and flipped
H = numpy.rot90(H)
H = numpy.flipud(H)
 
# Mask zeros
Hmasked = numpy.ma.masked_where(H==0,H) # Mask pixels with a value of zero
 
# Plot 2D histogram using pcolor
fig = plt.figure()
plt.pcolormesh(xedges,yedges,Hmasked)
plt.xlabel('Eta')
plt.ylabel('Phi')
cbar = plt.colorbar()
cbar.ax.set_ylabel('pT')

jets = numpy.load("../Data/jet_vars/our_jet_vars_"+str(options.event)+".npy")[0]
circles = []
for i, jet in enumerate(jets):
  if jet['pt'] < 20: continue
  circles.append(plt.Circle((jet['eta'],jet['phi']),0.4,color='r',fill=False))
  print jet['eta'],jet['phi'],jet['pt']

for circle in circles:
  fig.gca().add_artist(circle)

fig.savefig('../Plots/event_display_example.png')
