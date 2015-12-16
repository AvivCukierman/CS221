import math, numpy, copy, sys
from optparse import OptionParser
import operator 
import json

parser = OptionParser()
parser.add_option("-e", "--events",type=int, help='How many events to do for training and testing', default=10)
parser.add_option("-i", "--iterations",type=int, help='How many iterations to do SGD', default=10)
parser.add_option("-x", "--hs_factor",type=float, help='How much more we care about getting HS wrong than PU', default=1.0)
parser.add_option("-c", "--category",type=str, help='Which category of features to use', default='eta')
parser.add_option("-d", "--debug", action="store_true", default=False)
options, args = parser.parse_args()

numEvents = int(options.events)
debug = options.debug
category = options.category

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

def dotProduct(v1,v2):
  if len(v1.keys()) != len(v2.keys()): print 'Dot Product isn\'t right'
  result = 0
  for k in v1.keys():
    result+=v1[k]*v2[k]
  return result

def SGDupdate(w,eta,y,phi):
  if y==1: eta*=options.hs_factor
  for k in phi.keys():
    w[k] += eta*y*phi[k]

with open('feature_normalization_tjets_categories.json') as f:
  norm = json.load(f) 

import sys
if category not in norm:
    print 'Category doesn\'t exist!'
    sys.exit(0)

norm = norm[category]
norm['1'] = 1
features = norm.keys()
w = {k:0 for k in features}

import pdb
for it in range(options.iterations):
  eta = 2./numEvents/(it+1)
  #for i, particles in enumerate(particle_vars):
  for i in range(numEvents):
    #if i==numEvents: break
    particles_features = numpy.load("../Data/particle_features_tjets/particle_features_"+str(i)+".npy") #only look at particles within high pT jets
    if len(particles_features)==0: continue
    trainError = 0
    #print float(len([particle_features for particle_features in particles_features if particles[particle_features['index']]['truth']==1]))/len(particles_features)
    for particle_features in particles_features:
      y = 1 if particle_features['truth']==1 else -1
      norm_particle_features = {k: float(particle_features[k])/norm[k] for k in features}
      pt = particle_features['pt']
      prediction = dotProduct(w,norm_particle_features)
      margin = prediction*y
      if margin<0: trainError+=1
      #hinge loss
      if margin>1: continue
      else: SGDupdate(w,pt*eta,y,norm_particle_features)
    trainError = float(trainError)/len(particles_features)
    #print trainError

with open('../Output/weights_e'+str(options.events)+'_i'+str(options.iterations)+'_x'+str(options.hs_factor)+'_tjets_'+category+'.json','w') as g:
  json.dump(w,g)

print '../Output/SGD_results_e'+str(options.events)+'_i'+str(options.iterations)+'_x'+str(options.hs_factor)+'_tjets_'+category+'.txt'
f = open('../Output/SGD_results_e'+str(options.events)+'_i'+str(options.iterations)+'_x'+str(options.hs_factor)+'_tjets_'+category+'.txt','w')
for i in range(numEvents,2*numEvents):
  #particles = particle_vars[i]

  particles_features = numpy.load("../Data/particle_features_tjets/particle_features_"+str(i)+".npy") #only look at particles within high pT jets
  if len(particles_features)==0: continue
  testError = 0
  truthFraction = 0
  testHSError = 0
  testPUError = 0
  truthpt = 0
  totalpt = 0
  HSpt = 0
  totalPUpt = 0
  totalPUptretained = 0
  for particle_features in particles_features:
    y = 1 if particle_features['truth']==1 else -1
    pt = particle_features['pt'] 
    norm_particle_features = {k: float(particle_features[k])/norm[k] for k in features}
    prediction = dotProduct(w,norm_particle_features)
    margin = prediction*y
    if y>0:
        truthFraction +=1
        truthpt+=pt
    if y<0:
        totalPUpt+=pt
        if prediction>0: totalPUptretained+=pt
    if prediction>0:
        totalpt+=pt
        if y>0: HSpt+=pt
    if margin<0: 
      testError+=1
      if y>0: testHSError+=1
      else: testPUError+=1
  numParticles = len(particles_features)
  testError = float(testError)/numParticles
  testHSError = float(testHSError)/truthFraction
  testPUError = float(testPUError)/(numParticles - truthFraction)
  truthFraction = float(truthFraction)/numParticles
  totalpt /= truthpt
  HSpt /= truthpt
  if options.debug:
      print "testError:" + str(testError)
      print "testHSError:" + str(testHSError)
      print "testPUError:" + str(testPUError)
      print "truthFraction:" + str(truthFraction)
      print "response:" + str(totalpt)
      print "HS retained:" + str(HSpt)
      print "PU retained:" + str(totalPUptretained/totalPUpt)
      print '\n'

  f.write('testError:' + str(testError)+'\n')
  f.write('testHSError:' + str(testHSError)+'\n')
  f.write('testPUError:' + str(testPUError)+'\n')
  f.write('truthFraction:' + str(truthFraction)+'\n')
  f.write('response:' + str(totalpt)+'\n')
  f.write('HS retained:' + str(HSpt)+'\n')
  f.write('PU retained:' + str(totalPUptretained/totalPUpt)+'\n\n')

f.close()
