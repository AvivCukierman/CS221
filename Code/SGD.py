import math, numpy, copy, sys
from optparse import OptionParser
import operator 

parser = OptionParser()
parser.add_option("-e", "--events",type=int, help='How many events to do for training and testing', default=10)
parser.add_option("-i", "--iterations",type=int, help='How many iterations to do SGD', default=10)
parser.add_option("-d", "--debug", action="store_true", default=False)
options, args = parser.parse_args()

numEvents = int(options.events)
debug = options.debug

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
particle_vars = numpy.load("../Data/particle_vars.npy")
#sjet_vars = numpy.load("../Data/sjet_vars.npy")
#tjet_vars = numpy.load("../Data/tjet_vars.npy")

if debug: print len(particle_vars)

def dotProduct(v1,v2):
  result = 0
  for k in v1.keys():
    result+=v1[k]*v2[k]
  return result

def SGDupdate(w,eta,y,phi):
  for k in phi.keys():
    w[k] += eta*y*phi[k]


particles_features = numpy.load("../Data/particle_features/particle_features_0.npy")
features = particles_features[0].keys()
maxvalue = {k:max([abs(particle_features[k]) for particle_features in particles_features]) for k in features}
for k in features:
  if maxvalue[k] == 0: maxvalue[k]=1 #if they're all zero, doesn't matter
w = {k:0 for k in features}

import pdb
eta = 0.1
for it in range(options.iterations):
  for i, particles in enumerate(particle_vars):
    if i==numEvents: break
    particles_features = numpy.load("../Data/particle_features/particle_features_"+str(i)+".npy") #only look at particles within high pT jets
    trainError = 0
    #print float(len([particle_features for particle_features in particles_features if particles[particle_features['index']]['truth']==1]))/len(particles_features)
    for particle_features in particles_features:
      norm_particle_features = {k: float(particle_features[k])/maxvalue[k] for k in features}
      particle = particles[particle_features['index']]
      prediction = dotProduct(w,particle_features)
      y = 5 if particle['truth']==1 else -1
      margin = prediction*y
      if margin<0: trainError+=1
      #hinge loss
      if margin>1: continue
      else: SGDupdate(w,eta,y,particle_features)
    trainError = float(trainError)/len(particles_features)
    print trainError

print '\n'

for i,_ in enumerate(particle_vars):
  i+=numEvents
  if i==2*numEvents: break
  particles = particle_vars[i]

  particles_features = numpy.load("../Data/particle_features/particle_features_"+str(i)+".npy") #only look at particles within high pT jets
  testError = 0
  truthFraction = 0
  testHSError = 0
  testPUError = 0
  for particle_features in particles_features:
    particle = particles[particle_features['index']]
    prediction = dotProduct(w,particle_features)
    y = 1 if particle['truth']==1 else -1
    if y>0: truthFraction +=1
    margin = prediction*y
    if margin<0: 
      testError+=1
      if y>0: testHSError+=1
      else: testPUError+=1
  numParticles = len(particles_features)
  testError = float(testError)/numParticles
  testHSError = float(testHSError)/truthFraction
  testPUError = float(testPUError)/(numParticles - truthFraction)
  truthFraction = float(truthFraction)/numParticles
  print "testError:" + str(testError)
  print "testHSError:" + str(testHSError)
  print "testPUError:" + str(testPUError)
  print "truthFraction:" + str(truthFraction)
  print '\n'


#numpy.save("../Data/particle_features/particle_features_"+str(eventNum)+".npy", newparticles) 
