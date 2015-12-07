from optparse import OptionParser
import numpy as np
from operator import itemgetter

parser = OptionParser()

# job configuration
parser.add_option("--submitDir", help="dir to store the output", default="../Output")
parser.add_option("--inputParticles", help="input file containing particle vars", default="particle_vars.npy")
parser.add_option("--inputData", help="input folder containing data", default="../Data")
parser.add_option("--inputJets", help="input folder containing jet vars", default="jet_vars")
parser.add_option("--inputTJets", help="input folder contaning truth jet vars", default="tjet_vars")
parser.add_option("--events", help="input file containing event vars", default="event_vars.npy")
parser.add_option("--nevents", help="number of events to do", type=int, default=1)
parser.add_option("-x","--x", help="hard scatter factor", type=int, default=1)

(options, args) = parser.parse_args()

events = np.load(options.inputData+'/'+options.events)
all_particles = np.load(options.inputData+'/'+options.inputParticles)

efficiencies = []
multiple_matches_rates = []
jetpts = []
tjetpts = []
offsets = []
responses = []

if options.x==1: jet_filenames = [options.inputData+'/'+options.inputJets+'/our_jet_vars_'+str(i)+'.npy' for i in range(options.nevents)]
else: jet_filenames = [options.inputData+'/'+options.inputJets+'/jet_vars_'+str(i)+'_x'+str(options.x)+'.npy' for i in range(options.nevents)]
tjet_filenames = [options.inputData+'/'+options.inputTJets+'/our_tjet_vars_'+str(i)+'.npy' for i in range(options.nevents)]

i=0
for jet_filename,tjet_filename,event,particles in zip(jet_filenames,tjet_filenames,events,all_particles):
  #print i
  i+=1
  if i<=81: continue
  #print jet_filename,tjet_filename
  #print len(particles)

  if options.x==1: jets = np.load(jet_filename)[0]
  else: jets = np.load(jet_filename)
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
  #print efficiencies
  #print multiple_matches_rate

  for j in range(len(jets)):
    if not len(tjet_matches[j]) == 1: continue
    tjet_match = tjet_matches[j][0]
    tjetpt = tjets[tjet_match]['pt']
    jetpt = jets[j]['pt']
    if options.inputJets=='jet_vars':
      if jetpt<20: continue
    jetpts.append(jetpt)
    tjetpts.append(tjetpt)
    offsets.append(jetpt-tjetpt)
    responses.append(jetpt/tjetpt)

np.save(options.submitDir+'/'+options.inputJets+'_x'+str(options.x)+'_efficiencies.npy',np.array(efficiencies))
np.save(options.submitDir+'/'+options.inputJets+'_x'+str(options.x)+'_multiple_matches.npy',np.array(multiple_matches_rates))
np.save(options.submitDir+'/'+options.inputJets+'_x'+str(options.x)+'_jetpts.npy',np.array(jetpts))
np.save(options.submitDir+'/'+options.inputJets+'_x'+str(options.x)+'_tjetpts.npy',np.array(tjetpts))
np.save(options.submitDir+'/'+options.inputJets+'_x'+str(options.x)+'_offsets.npy',np.array(offsets))
np.save(options.submitDir+'/'+options.inputJets+'_x'+str(options.x)+'_responses.npy',np.array(responses))
