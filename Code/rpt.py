from optparse import OptionParser
import numpy as np
from operator import itemgetter

parser = OptionParser()

# job configuration
parser.add_option("--submitDir", help="dir to store the output", default="../Output")
parser.add_option("--inputParticles", help="input file containing particle vars", default="../Data/particle_vars.npy")
parser.add_option("--inputJets", help="input file containing jet vars", default="../Data/jet_vars")
parser.add_option("--events", help="input file containing event vars", default="../Data/event_vars.npy")
parser.add_option("--nevents", help="number of events to do", type=int, default=1)
(options, args) = parser.parse_args()

events = np.load(options.events)
all_particles = np.load(options.inputParticles)

charged_rpts = []
all_rpts = []

jet_filenames = [options.inputJets+'/our_jet_vars_'+str(i)+'.npy' for i in range(options.nevents)]

for jet_filename,event,particles in zip(jet_filenames,events,all_particles):

  print jet_filename

  jets = np.load(jet_filename)[0]

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

np.save(options.submitDir+'/charged_rpts.npy',np.array(charged_rpts))
np.save(options.submitDir+'/all_rpts.npy',np.array(all_rpts))
