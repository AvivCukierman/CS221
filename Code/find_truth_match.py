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
parser.add_option("--nevents", help="number of events to do", type=int, default=1000)

(options, args) = parser.parse_args()

all_particles = np.load(options.inputData+'/'+options.inputParticles)


jet_filenames = [options.inputData+'/'+options.inputJets+'/our_jet_vars_'+str(i)+'.npy' for i in range(options.nevents)]
tjet_filenames = [options.inputData+'/'+options.inputTJets+'/our_tjet_vars_'+str(i)+'.npy' for i in range(options.nevents)]

all_HS_tjet_matches = []
for jet_filename,tjet_filename,particles in zip(jet_filenames,tjet_filenames,all_particles):
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

    HS_tjet_matches = [] 
    for jet_index,jet in enumerate(jets):
        tjets_matched = [tjets[m] for m in tjet_matches[jet_index]]
        if not len(tjets_matched)==1: continue
        if tjets_matched[0]['pt']<20: continue
        HS_tjet_matches.append((jet_index,tjet_matches[jet_index][0]))
    print HS_tjet_matches

    all_HS_tjet_matches.append(HS_tjet_matches)

np.save(options.submitDir+'/tjet_matches.npy',np.array(all_HS_tjet_matches))
