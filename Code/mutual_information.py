import math, numpy, copy, sys
from optparse import OptionParser
import operator 
import json
from anti_kt_algorithm import antikt_algorithm

parser = OptionParser()
parser.add_option("-s", "--start_event",type=int, help='Which event to start on', default=200)
parser.add_option("-e", "--end_event",type=int, help='Which event to end on', default=201)
parser.add_option( "--submitDir",type=str, help='Where to store the output', default='../Output')
options, args = parser.parse_args()

with open('feature_normalization_tjets.json') as f:
  norm = json.load(f) 
features = norm.keys()

#particle_vars = numpy.load("../Data/particle_vars.npy")

def mutual_info(event):
    particles_features = numpy.load("../Data/particle_features_tjets/particle_features_"+str(event)+".npy") #only look at particles within high pT jets

    jetpts = []
    tjetpts = []
    offsets = []
    responses = []

    N = len(particles_features)
    N00 = {k:0 for k in features if 'gt' in k} 
    N01 = {k:0 for k in features if 'gt' in k} #truth = 0, indicator = 1
    N10 = {k:0 for k in features if 'gt' in k} #truth = 1, indicator = 0
    N11 = {k:0 for k in features if 'gt' in k}
    for particle_features in particles_features:
        y = 1 if particle_features['truth'] else 0
        norm_particle_features = {k: float(particle_features[k])/norm[k] for k in features}
        #print norm_particle_features
        for k in features:
          if 'gt' not in k: continue
          if norm_particle_features[k]:
              if y: N11[k]+=1
              else: N01[k]+=1
          else:
              if y: N10[k]+=1
              else: N00[k]+=1

    return N,N00,N01,N10,N11 

N = 4 
N00 = {k:1 for k in features if 'gt' in k} #lambda smoothing
N01 = {k:1 for k in features if 'gt' in k} #truth = 0, indicator = 1
N10 = {k:1 for k in features if 'gt' in k} #truth = 1, indicator = 0
N11 = {k:1 for k in features if 'gt' in k}
for event in range(options.start_event,options.end_event):
    n,n00,n01,n10,n11 = mutual_info(event)
    N+=n
    for k in features:
        if 'gt' not in k: continue
        N00[k]+=n00[k]
        N01[k]+=n01[k]
        N10[k]+=n10[k]
        N11[k]+=n11[k]

mi = {k:0 for k in features if 'gt' in k}
for k in features:
    if 'gt' in k:
        n00 = float(N00[k])
        n01 = float(N01[k])
        n10 = float(N10[k])
        n11 = float(N11[k])
        if n00+n10==0:
            print 'all 0:'+k
            continue
        if n01+n11==0:
            print 'all 1:'+k
            continue
        term11 = n11/N*math.log((N*n11)/((n11+n10)*(n11+n01)))
        term01 = n01/N*math.log((N*n01)/((n01+n00)*(n11+n01)))
        term10 = n10/N*math.log((N*n10)/((n10+n11)*(n00+n10)))
        term00 = n00/N*math.log((N*n00)/((n00+n01)*(n00+n10)))
        mi[k] = term11+term01+term10+term00

#print mi
f=open('mutual_information.txt','w')
for k in sorted(mi,key=mi.get,reverse=True):
    f.write(str(k)+','+str(mi[k])+'\n')
f.close()

#print options.submitDir+'/'+filename+'_x'+str(options.hs_factor)+'_jetpts.npy'
#numpy.save(options.submitDir+'/'+filename+'_x'+str(options.hs_factor)+'_jetpts.npy',numpy.array(jetpts))
