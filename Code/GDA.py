import math, numpy, copy, sys
from optparse import OptionParser
import operator 
import json
from numpy.linalg import inv
import numpy as np
from sklearn.lda import LDA
from sklearn.qda import QDA
import pickle
import random

parser = OptionParser()
parser.add_option("-p", "--particles",type=int, help='How many particles to use for training and testing', default=100000)
parser.add_option("-d", "--debug", action="store_true", default=False)
parser.add_option("-L", "--Lambda", type=float, help="Laplace smoothing", default = 0.5)
parser.add_option("-m", "--multiplier", type=float, help="multiply probabilities", default = 10)
options, args = parser.parse_args()

numParticles = int(options.particles)
debug = options.debug
Lambda = options.Lambda
multiplier = options.multiplier

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
for k in range(len(features)):
    if features[k] == 'pt':
        ptIndex = k

# train_data = []
# train_truth = []
# 
# for i in range(numParticles):
#     particles_features = numpy.load("../Data/particle_features_tjets/particle_features_"+str(i)+".npy") #only look at particles within high pT jets
#     if len(particles_features)==0: continue
# 
#     for particle_features in particles_features:
#         y = 1 if particle_features['truth']==1 else -1
#         norm_particle_features = [float(particle_features[features[k]])/norm[features[k]] for k in range(len(features))]
# #        pt = particle_features['pt']
#         train_data.append(norm_particle_features)
#         train_truth.append(y)
# 
# train_data = np.array(train_data)
# train_truth = np.array(train_truth)
# np.save("../Data/particle_features_tjets/particle_features.npy", train_data)
# np.save("../Data/particle_features_tjets/particle_features_truth.npy", train_truth)

data = np.load("../Data/particle_features_tjets/particle_features.npy")
truth = np.load("../Data/particle_features_tjets/particle_features_truth.npy")

print len(data),len(truth)

train_data = data[0:numParticles]
train_truth = truth[0:numParticles]
test_data = data[numParticles:2*numParticles]
test_truth = truth[numParticles:2*numParticles]

maxPt = max(particle[ptIndex] for particle in train_data)
keptParticleIndices = []
for i in range(numParticles):
#    if train_truth[i]==1 or random.uniform(0,1) < multiplier*(train_data[i][ptIndex]+Lambda)/(maxPt + Lambda):
    if random.uniform(0,1) < multiplier*(train_data[i][ptIndex]+Lambda)/(maxPt + Lambda):
        keptParticleIndices.append(i)
        
train_data = train_data[keptParticleIndices]
train_truth = train_truth[keptParticleIndices]
print len(train_truth)

clf = LDA()
clf.fit(train_data,train_truth)

# output=open("lda.pkl",'wb')
# pickle.dump(clf,output)
# output.close()

# test_data = []
# test_truth = []
# 
# for i in range(numParticles,2*numParticles):
#     particles_features = numpy.load("../Data/particle_features_tjets/particle_features_"+str(i)+".npy") #only look at particles within high pT jets
#     if len(particles_features)==0: continue
# 
#     for particle_features in particles_features:
#         y = 1 if particle_features['truth']==1 else -1
#         norm_particle_features = [float(particle_features[features[k]])/norm[features[k]] for k in range(len(features))]
# #        pt = particle_features['pt']
#         test_data.append(norm_particle_features)
#         test_truth.append(y)
# 
# test_data = np.array(test_data)
# test_truth = np.array(test_truth)
# 
# np.save("../Data/particle_features_tjets/particle_features_"+str(numParticles)+"_upto_"+str(2*numParticles)+".npy", test_data)
# np.save("../Data/particle_features_tjets/particle_features_"+str(numParticles)+"_upto_"+str(2*numParticles)+"_truth.npy", test_truth)

#print train_data.shape, train_truth.shape
#print test_data.shape, test_truth.shape
print(clf.score(train_data,train_truth))
print(clf.score(test_data,test_truth))

p = np.where(train_truth == -1)[0]
p_truth = np.where(test_truth == -1)[0]
hs = np.where(train_truth == 1)[0]
hs_truth = np.where(test_truth == 1)[0]
print("Pileup")
print(clf.score(train_data[p],train_truth[p]))
print(clf.score(test_data[p_truth],test_truth[p_truth]))
print("Hard Scatter")
print(clf.score(train_data[hs],train_truth[hs]))
print(clf.score(test_data[hs_truth],test_truth[hs_truth]))

print '\n'

"""
numTrue = 0
numFalse = 0
numTotal = 0
#these three are unnormalized; we normalize at the end
muFalse = {k:0 for k in features}
muTrue = {k:0 for k in features}
sigma = {(k1,k2):0 for k1 in features for k2 in features}

for i in range(numParticles):
    if debug: print i
    
    particles_features = numpy.load("../Data/particle_features_tjets/particle_features_"+str(i)+".npy") #only look at particles within high pT jets
    if len(particles_features)==0: continue

    for particle_features in particles_features:
        y = 1 if particle_features['truth']==1 else -1
        norm_particle_features = {k: float(particle_features[k])/norm[k] for k in features}
        pt = particle_features['pt']

        numTotal += 1
        if y == 1:
            numTrue += 1
            for k in norm_particle_features.keys():
                muTrue[k] += norm_particle_features[k]
        else:
            numFalse += 1
            for k in norm_particle_features.keys():
                muFalse[k] += norm_particle_features[k]

for k in muTrue.keys():
    muTrue[k] /= float(numTrue)
for k in muFalse.keys():
    muFalse[k] /= float(numFalse)
phi = float(numTrue) / numTotal

good_features=[]
for i in range(numParticles):
    particles_features = numpy.load("../Data/particle_features_tjets/particle_features_"+str(i)+".npy") #only look at particles within high pT jets
    if len(particles_features)==0: continue

    for particle_features in particles_features:
        y = 1 if particle_features['truth']==1 else -1
        norm_particle_features = {k: float(particle_features[k])/norm[k] for k in features}
        pt = particle_features['pt']

        residual = {k: norm_particle_features[k] - muTrue[k] for k in features}
        for k1 in features:
            alwaysZero = True
            for k2 in features:
                sigma[(k1,k2)] += residual[k1]*residual[k2]
                if sigma[(k1,k2)] != 0:
                    alwaysZero = False
            if not alwaysZero:
                good_features.append(k1)

for k1 in good_features:
    for k2 in good_features:
        sigma[(k1,k2)] /= float(numTotal)
        
sigmaList = [[0 for good_feature in good_features] for good_feature in good_features]
for i in range(len(good_features)):
    for j in range(len(good_features)):
        sigmaList[i][j] = sigma[(good_features[i],good_features[j])]
sigmaMat = np.matrix(sigmaList)
print "blah"
sigmaInv = inv(sigmaMat)
print "blahblah"

with open('../Output/gda_e'+str(options.events)+'_x'+str(options.hs_factor)+'_muTrue_tjets.json','w') as g:
  json.dump(muTrue,g)
with open('../Output/gda_e'+str(options.events)+'_x'+str(options.hs_factor)+'_muFalse_tjets.json','w') as g:
  json.dump(muFalse,g)
with open('../Output/gda_e'+str(options.events)+'_x'+str(options.hs_factor)+'_sigma_tjets.json','w') as g:
  json.dump(sigma,g)
with open('../Output/gda_e'+str(options.events)+'_x'+str(options.hs_factor)+'_phi_tjets.json','w') as g:
  json.dump(phi,g)
"""

print '../Output/gda_results_p'+str(numParticles)+'_L'+str(Lambda)+'_m'+str(multiplier)+'_tjets.txt'
f = open('../Output/gda_results_p'+str(numParticles)+'_L'+str(Lambda)+'_m'+str(multiplier)+'_tjets.txt','w')

testError = 0
truthFraction = 0
testHSError = 0
testPUError = 0
truthpt = 0
totalpt = 0
HSpt = 0
totalPUpt = 0
totalPUptretained = 0
for i in range(len(test_data)):
    predictFalse, predictTrue = clf.predict_proba(test_data[i])[0]

    pt = test_data[i][ptIndex]
    y = test_truth[i]
    sign = 1 if predictTrue > 0.5 else -1
    prediction = predictTrue * sign
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

numParticles = len(test_data)
testError = float(testError)/numParticles
testHSError = float(testHSError)/truthFraction
testPUError = float(testPUError)/(numParticles - truthFraction)
truthFraction = float(truthFraction)/numParticles
totalpt /= truthpt
HSpt /= truthpt
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
