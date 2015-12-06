import math, numpy, copy, sys
from optparse import OptionParser
import operator 
import json

parser = OptionParser()
parser.add_option("-e", "--event",type=int, help='Which event to examine', default=50)
options, args = parser.parse_args()

def dotProduct(v1,v2):
  result = 0
  for k in v1.keys():
    result+=v1[k]*v2[k]
    print k,v1[k],v2[k],result
  return result

particles_features = numpy.load("../Data/particle_features/particle_features_0.npy")
features = particles_features[0].keys()
maxvalue = {k:max([abs(particle_features[k]) for particle_features in particles_features]) for k in features}
for k in features:
  if maxvalue[k] == 0: maxvalue[k]=1 #if they're all zero, doesn't matter
with open('../Output/weights_e40_i10_x200.json') as file:
  w = json.load(file)

particle_vars = numpy.load("../Data/particle_vars.npy")

for i, particles in enumerate(particle_vars):
  if i!=options.event: continue 
  particles_features = numpy.load("../Data/particle_features/particle_features_"+str(i)+".npy") #only look at particles within high pT jets
  trainError = 0
  #print float(len([particle_features for particle_features in particles_features if particles[particle_features['index']]['truth']==1]))/len(particles_features)
  for particle_features in particles_features:
    norm_particle_features = {k: float(particle_features[k])/maxvalue[k] for k in features}
    particle = particles[particle_features['index']]
    prediction = dotProduct(w,particle_features)
    print particle['truth']
    break

