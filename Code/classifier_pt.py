import math, numpy, copy, sys
from optparse import OptionParser
import operator 
import json
from anti_kt_algorithm import antikt_algorithm

parser = OptionParser()
parser.add_option("-s", "--start_event",type=int, help='Which event to start on', default=0)
parser.add_option("-e", "--end_event",type=int, help='Which event to end on', default=1)
parser.add_option('-w', action='store_true', default=False, dest='weighted',help = 'Use weighted classification as opposed to binary classification')
parser.add_option("-x", "--hs_factor",type=float, help='How much more we value identifying HS correctly than PU (hyperparameter)', default=1.0)
parser.add_option( "--submitDir",type=str, help='Where to store the output', default='../Output')
options, args = parser.parse_args()

def dotProduct(v1,v2):
  result = 0
  for k in v1.keys():
    result+=v1[k]*v2[k]
    #print k,v1[k],v2[k],result
  return result

def weight(z):
  return 1/(1+math.exp(-z))


with open('feature_normalization_tjets_categories.json') as f:
  norm = json.load(f) 
norm = norm['pt']
norm['1']=1
features = norm.keys()

#particle_vars = numpy.load("../Data/particle_vars.npy")

def makejets(hs_factor,event,weighted):
    if hs_factor>0:
      with open('../Output/weights_e200_i10_x'+str(hs_factor)+'_tjets_pt.json') as file:
          w = json.load(file)
    particles_features = numpy.load("../Data/particle_features_tjets/particle_features_"+str(event)+".npy") #only look at particles within high pT jets

    #jetpts = []

    for particle_features in particles_features:
        norm_particle_features = {k: float(particle_features[k])/norm[k] for k in features}
        if hs_factor>0:
          prediction = dotProduct(w,norm_particle_features)
        print particle_features['pt'],weight(prediction)


    return 

          #y = 1 if particle['truth']==1 else -1
          #if y>0: truthFraction +=1
          #margin = prediction*y

          #if margin<0: 
            #testError+=1
            #if y>0: testHSError+=1
            #else: testPUError+=1
        #numParticles = len(particles_features)
        #testError = float(testError)/numParticles
        #testHSError = float(testHSError)/truthFraction
        #testPUError = float(testPUError)/(numParticles - truthFraction)
        #truthFraction = float(truthFraction)/numParticles
        #print "testError:" + str(testError)
        #print "testHSError:" + str(testHSError)
        #print "testPUError:" + str(testPUError)
        #print "truthFraction:" + str(truthFraction)

    #if weighted: outfilename = '../Data/weighted_classifier_jet_vars/jet_vars_'+str(event)+'_x'+str(hs_factor)+'.npy'
    #else: outfilename = '../Data/binary_classifier_jet_vars/jet_vars_'+str(event)+'_x'+str(hs_factor)+'.npy'
#print outfilename
    #numpy.save(outfilename,jets)

for event in range(options.start_event,options.end_event):
    makejets(options.hs_factor,event,options.weighted)


#print options.submitDir+'/'+filename+'_x'+str(options.hs_factor)+'_jetpts.npy'
#numpy.save(options.submitDir+'/'+filename+'_x'+str(options.hs_factor)+'_jetpts.npy',numpy.array(jetpts))
