import math, numpy, copy, sys
from optparse import OptionParser
import operator 
import json
from anti_kt import antikt_algorithm

parser = OptionParser()
parser.add_option("-e", "--event",type=int, help='Which event to examine', default=100)
parser.add_option("-w", "--weighted",action="store_true", help='Use classification as weight on pT rather than binary', default=False)
options, args = parser.parse_args()

def dotProduct(v1,v2):
  result = 0
  for k in v1.keys():
    result+=v1[k]*v2[k]
    #print k,v1[k],v2[k],result
  return result

def weight(z):
  if z>0: return 1
  return 2/(1+math.exp(-2*z))

with open('feature_normalization.json') as f:
  norm = json.load(f) 
features = norm.keys()
with open('../Output/weights_e80_i10_x2.json') as file:
  w = json.load(file)

ghost_pt = 1e-100

particle_vars = numpy.load("../Data/particle_vars.npy")

for i, particles in enumerate(particle_vars):
  if i!=options.event: continue 
  particles_features = numpy.load("../Data_notDropBox/particle_features/particle_features_"+str(i)+".npy") #only look at particles within high pT jets
  trainError = 0
  #print float(len([particle_features for particle_features in particles_features if particles[particle_features['index']]['truth']==1]))/len(particles_features)
  testError = 0
  truthFraction = 0
  testHSError = 0
  testPUError = 0

  new_particles = []
  for particle_features in particles_features:
    norm_particle_features = {k: float(particle_features[k])/norm[k] for k in features}
    index = particle_features['index']
    particle = particles[index]
    prediction = dotProduct(w,norm_particle_features)

    newparticle = particle
    newparticle['index'] = index 
    if not options.weighted:
      if prediction<0: newparticle['pt'] = ghost_pt #binary classification
    else:
      if prediction<0: newparticle['pt'] *= weight(prediction) #weighted classification
    new_particles.append(newparticle)

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

  jets = antikt_algorithm(new_particles)

  if options.weighted: outfilename = '../Data/weighted_classifier_jet_vars/jet_vars_'+str(options.event)+'.npy'
  else: outfilename = '../Data/binary_classifier_jet_vars/jet_vars_'+str(options.event)+'.npy'
  print outfilename
  numpy.save(outfilename,jets)
