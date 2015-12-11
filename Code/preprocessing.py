import math, numpy, copy, sys
from optparse import OptionParser
import operator 
import json
from numpy.linalg import inv
import numpy as np
from sklearn.lda import LDA
import pickle


numParticles = 1000

# Each file is a numpy array that should have N = 1000 entries, which is the number of events we can use for now
# 
# event_vars.npy:
# Each element of the array is [NPV,rho,sigma_rho] for each event.
# 
# particle_vars.npy:
# Each element of the array is an array which is N_particles long (which varies from event to event). Each element of that array is a dictionary of the following form:
# {'pt':_,'eta':_,'phi':_,'m':_,'truth':_,'pu':_,'charge':_}
# 'pt','eta','phi','m' are self-explanatory.
# 'truth' is 1 if particle is from the HS vertex and 0 if from a pileup vertex
# 'pu' indicates which vertex the particle came from - 0 is HS vertex, > 0 is a pileup vertex (we probably won't use this variable)
# 'charge' indicates the charge of the particle, a member of {-1,0,+1}. Keep in mind that at ATLAS we have access to the variable 'pu' (and therefore 'truth') for charged particles. I.e., the tracks are the charged particles.
# 
# I ran standard anti-kt4 jet finding on the particles in each event. 
# jet_vars.npy:
# The jets from anti-kt4. Each element of the array is an array which is N_jets long. Each element of that array is a dictionary of the following form:
# {'pt':_,'eta':_,'phi':_,'m':_,'area':_,'width':_}
# 'area' is the area associated with the jet
# 'width' is the sum of (pt_particle)*(dr_particle) for each particle in the jet, where dr is measured from the jet axis (bigger is more spread out)
# 
# sjet_vars.npy:
# The same as jet_vars.npy, but this time with pt -> pt-rho*area.
# 
# tjet_vars.npy:
# The same as jet_vars.npy, but this time jet finding was only run on particles with 'truth'==1.

#event_vars = numpy.load("../Data/event_vars.npy")
#particle_vars = numpy.load("../Data/particle_vars.npy")
#sjet_vars = numpy.load("../Data/sjet_vars.npy")
#tjet_vars = numpy.load("../Data/tjet_vars.npy")

#if debug: print len(particle_vars)


with open('feature_normalization_tjets.json') as f:
  norm = json.load(f) 
features = norm.keys()

data = []
truth = []

for i in range(numParticles):
    particles_features = numpy.load("../Data/particle_features_tjets/particle_features_"+str(i)+".npy") #only look at particles within high pT jets
    if len(particles_features)==0: continue

    for particle_features in particles_features:
        y = 1 if particle_features['truth']==1 else -1
        norm_particle_features = [float(particle_features[features[k]])/norm[features[k]] for k in range(len(features))]
#        pt = particle_features['pt']
        data.append(norm_particle_features)
        truth.append(y)

data = np.array(data)
truth = np.array(truth)
np.save("../Data/particle_features_tjets/particle_features.npy", data)
np.save("../Data/particle_features_tjets/particle_features_truth.npy", truth)

print len(truth)