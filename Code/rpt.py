from optparse import OptionParser
import numpy as np
from operator import itemgetter

parser = OptionParser()

# job configuration
parser.add_option("--submitDir", help="dir to store the output", default="output")
parser.add_option("--inputParticles", help="input file containing particle vars", default="particle_vars.npy")
parser.add_option("--inputJets", help="input file containing jet vars", default="jet_vars.npy")
parser.add_option("--events", help="input file containing event vars", default="event_vars.npy")
parser.add_option("--nevents", help="number of events to do", default=1)

(options, args) = parser.parse_args()

events = np.load(options.events)
all_particles = np.load(options.inputParticles)
all_jets = np.load(options.inputJets)

i=0
charged_rpts = []
all_rpts = []


for event,particles,jets in zip(events,all_particles,all_jets):
  for jet_index,jet in enumerate(jets):
    charged_pt_from_HS = 0
    charged_pt = 0
    all_pt_from_HS = 0
    particle_indices = jet['particle_indices']
    for particle_index in particle_indices:
      particle = particles[particle_index]
      if particle['charge']==1: charged_pt += particle['pt']
      if particle['truth']==1:
        if particle['charge']==1: charged_pt_from_HS += particle['pt']
        all_pt_from_HS += particle['pt']
    charged_rpt = charged_pt_from_HS/charged_pt if charged_pt>0 else -1
    charged_rpts.append(charged_rpt)
    all_rpts.append(all_pt_from_HS/jet['pt'])

  i+=1
  if i>=options.nevents: break 

np.save(options.submitDir+'/charged_rpts.npy',np.array(charged_rpts))
np.save(options.submitDir+'/all_rpts.npy',np.array(all_rpts))
