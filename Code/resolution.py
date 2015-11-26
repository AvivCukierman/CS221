from optparse import OptionParser
import numpy as np
from operator import itemgetter

parser = OptionParser()

# job configuration
parser.add_option("--submitDir", help="dir to store the output", default="../Output")
parser.add_option("--inputParticles", help="input file containing particle vars", default="../Data/particle_vars.npy")
parser.add_option("--inputJets", help="input file containing jet vars", default="../Data/jet_vars")
parser.add_option("--inputTJets", help="input file contaning truth jet vars", default="../Data/tjet_vars")
parser.add_option("--events", help="input file containing event vars", default="../Data/event_vars.npy")
parser.add_option("--nevents", help="number of events to do", type=int, default=1)

(options, args) = parser.parse_args()

events = np.load(options.events)
all_particles = np.load(options.inputParticles)

efficiencies = []
multiple_matches_rates = []
jetpts = []
tjetpts = []
offsets = []

jet_filenames = [options.inputJets+'/our_jet_vars_'+str(i)+'.npy' for i in range(options.nevents)]
tjet_filenames = [options.inputTJets+'/our_tjet_vars_'+str(i)+'.npy' for i in range(options.nevents)]

for jet_filename,tjet_filename,event,particles in zip(jet_filenames,tjet_filenames,events,all_particles):
  print jet_filename,tjet_filename
  print len(particles)

  jets = np.load(jet_filename)[0]
  tjets = np.load(tjet_filename)[0]

  primary_particles = []
  for tjet_index,tjet in enumerate(tjets):
    particle_indices = tjet['particle_indices']
    particlepts = [(i,particles[i]['pt']) for i in particle_indices]
    primary_particle = max(particlepts,key=itemgetter(1))[0]
    primary_particles.append((tjet_index,primary_particle))

  tjet_matches = {j:[] for j in range(len(jets))}
  for jet_index,jet in enumerate(jets):
    particle_indices = jet['particle_indices']
    for tjet_index,primary_particle in primary_particles:
      if primary_particle in particle_indices: tjet_matches[jet_index].append(tjet_index) 

  efficiency = float(sum([len(tjet_matches[j]) for j in range(len(jets))]))/len(tjets)
  multiple_matches_rate = float(len([j for j in range(len(jets)) if len(tjet_matches[j])>1]))/len(jets)
  efficiencies.append(efficiency)
  multiple_matches_rates.append(multiple_matches_rate)
  print efficiencies
  print multiple_matches_rate

  for j in range(len(jets)):
    if not len(tjet_matches[j]) == 1: continue
    tjet_match = tjet_matches[j][0]
    tjetpt = tjets[tjet_match]['pt']
    jetpt = jets[j]['pt']
    jetpts.append(jetpt)
    tjetpts.append(tjetpt)
    offsets.append(jetpt-tjetpt)

np.save(options.submitDir+'/efficiencies.npy',np.array(efficiencies))
np.save(options.submitDir+'/multiple_matches.npy',np.array(multiple_matches_rates))
np.save(options.submitDir+'/jetpts.npy',np.array(jetpts))
np.save(options.submitDir+'/tjetpts.npy',np.array(tjetpts))
np.save(options.submitDir+'/offsets.npy',np.array(offsets))
