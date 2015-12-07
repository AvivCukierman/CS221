import math, numpy, copy, sys
from optparse import OptionParser
import operator 
import json
import matplotlib.pyplot as plt

parser = OptionParser()
parser.add_option("-e", "--event",type=int, help='Which event to examine', default=100)
parser.add_option("-p", "--truth_particles",action="store_true", help='Only plot truth particles', default=False)
parser.add_option("-j", "--truth_jets",action="store_true", help='Only plot truth jets', default=False)
options, args = parser.parse_args()

particle_vars = numpy.load("../Data/particle_vars.npy")

etas = []
phis = []
pts = []

for i, particles in enumerate(particle_vars):
  if i!=options.event: continue 
  for particle in particles:
    eta = particle['eta']
    phi = particle['phi']
    if options.truth_particles:
      if not particle['truth']: continue
    etas.append(particle['eta'])
    phis.append(particle['phi'])
    pts.append(particle['pt'])

if options.truth_jets: jets = numpy.load("../Data/tjet_vars/our_tjet_vars_"+str(options.event)+".npy")[0]
else: jets = numpy.load("../Data/jet_vars/our_jet_vars_"+str(options.event)+".npy")[0]
circles = []
for i, jet in enumerate(jets):
  if jet['pt'] < 20: continue
  circles.append(plt.Circle((jet['eta'],jet['phi']),0.4,color='r',fill=False))
  print jet['eta'],jet['phi'],jet['pt']

fig2 = plt.figure(1)
plt.scatter(etas,phis,c=pts)
plt.xlabel('Eta')
plt.ylabel('Phi')
plt.xlim(-3,3)
plt.ylim(0,6.29)
for circle in circles:
  fig2.gca().add_artist(circle)
figname = 'event_display_scatter_e'+str(options.event)
if options.truth_particles: figname+='_tp'
fig2.savefig('../Plots/'+figname+'.png')

nbins = 50 
H,xedges,yedges = numpy.histogram2d(etas,phis,bins=nbins,weights=pts)

# H needs to be rotated and flipped
H = numpy.rot90(H)
H = numpy.flipud(H)
 
# Mask zeros
Hmasked = numpy.ma.masked_where(H==0,H) # Mask pixels with a value of zero
 
# Plot 2D histogram using pcolor
fig = plt.figure(2)
plt.pcolormesh(xedges,yedges,Hmasked)
plt.xlabel('Eta')
plt.ylabel('Phi')
plt.xlim(-3,3)
plt.ylim(0,6.29)
cbar = plt.colorbar()
cbar.ax.set_ylabel('pT')

circles = []
for i, jet in enumerate(jets):
  if jet['pt'] < 20: continue
  circles.append(plt.Circle((jet['eta'],jet['phi']),0.4,color='r',fill=False))
  print jet['eta'],jet['phi'],jet['pt']
for circle in circles:
  fig.gca().add_artist(circle)

figname = 'event_display_hist_e'+str(options.event)
if options.truth_particles: figname+='_tp'
fig.savefig('../Plots/'+figname+'.png')
