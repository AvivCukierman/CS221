import math, numpy, copy, sys
from optparse import OptionParser
import operator 
import json
from anti_kt_algorithm import antikt_algorithm

parser = OptionParser()
parser.add_option("-s", "--start_event",type=int, help='Which event to start on', default=200)
parser.add_option("-e", "--end_event",type=int, help='Which event to end on', default=201)
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

def combine(ptI,ptJ,etaI,etaJ,phiI,phiJ,mI,mJ):
  pI = ptI * math.cosh(etaI)
  pJ = ptJ * math.cosh(etaJ)
  #assume massless or don't
  EI = math.sqrt(pI**2 + mI**2)
  EJ = math.sqrt(pJ**2 + mJ**2)

  #E-scheme: 4-momenta add. Other schemes (like weighted averages of phi and eta)
  #exist and perform similarly -- see kT paper by Ellis and Soper.
  #For list and description of other schemes, see FastJet manual, section 3.4
  pxI = ptI * math.cos(phiI)
  pyI = ptI * math.sin(phiI)
  pzI = ptI * math.sinh(etaI)
  pxJ = ptJ * math.cos(phiJ)
  pyJ = ptJ * math.sin(phiJ)
  pzJ = ptJ * math.sinh(etaJ)
  newPx = pxI + pxJ
  newPy = pyI + pyJ
  newPz = pzI + pzJ
  newPT = math.sqrt(newPx**2 + newPy**2)
  newP = math.sqrt(newPT**2 + newPz**2)
  newE = EI + EJ
  #note: if we use acosh(newP/newPT) we lose sign info
  newEta = math.asinh(newPz/newPT)
  newM = math.sqrt(max(0.0, newE**2 - newP**2))
  
  c = newPx / newPT
  s = newPy / newPT
  if c >= 0:
      if s >= 0:
          newPhi = math.asin(s)
      else:
          newPhi = 2*math.pi + math.asin(s)
  else:
      if s >= 0:
          newPhi = math.acos(c)
      else:
          newPhi = 2*math.pi - math.acos(c)
  return newPT,newEta,newPhi,newM

def combine_all(particles):
  newPT = particles[0]['pt']
  newEta = particles[0]['eta']
  newPhi = particles[0]['phi']
  newM = particles[0]['m']

  for i,particle in enumerate(particles):
    if i==0: continue
    newPT,newEta,newPhi,newM = combine(newPT,particle['pt'],newEta,particle['eta'],newPhi,particle['phi'],newM,particle['m'])

  return {'pt':newPT,'eta':newEta,'phi':newPhi,'m':newM}


with open('feature_normalization_tjets.json') as f:
  norm = json.load(f) 
features = norm.keys()

particle_vars = numpy.load("../Data/particle_vars.npy")

def makejets(hs_factor,event,weighted):
    particles = particle_vars[event]
    tjet_matches = numpy.load('../Output/tjet_matches.npy')[event]
    jet_vars = numpy.load('../Data/jet_vars/our_jet_vars_'+str(event)+'.npy')[0]
    tjet_vars = numpy.load('../Data/tjet_vars/our_tjet_vars_'+str(event)+'.npy')[0]
    if hs_factor>0:
      with open('../Output/weights_e200_i10_x'+str(hs_factor)+'_tjets.json') as file:
          w = json.load(file)
    particles_features = numpy.load("../Data/particle_features_tjets/particle_features_"+str(event)+".npy") #only look at particles within high pT jets

    ghost_pt = 1e-100

    #trainError = 0
#print float(len([particle_features for particle_features in particles_features if particles[particle_features['index']]['truth']==1]))/len(particles_features)
    #testError = 0
    #truthFraction = 0
    #testHSError = 0
    #testPUError = 0

    jetpts = []
    tjetpts = []
    offsets = []
    responses = []

    for jet_index,tjet_index in tjet_matches:
        jet = jet_vars[jet_index]
        indices = jet['particle_indices']

        new_particles = []
        for particle_features in particles_features:
          index = particle_features['index']
          if index not in indices: continue
          norm_particle_features = {k: float(particle_features[k])/norm[k] for k in features}
          particle = particles[index]
          if hs_factor>0:
            prediction = dotProduct(w,norm_particle_features)
          else:
            prediction = 1

          newparticle = particle
          newparticle['index'] = index 
          if not weighted:
            if prediction<0: newparticle['pt'] = ghost_pt #binary classification
          else:
            if prediction<0: newparticle['pt'] *= weight(prediction) #weighted classification
          new_particles.append(newparticle)

        newjet = combine_all(new_particles)
        print newjet,jet
        tjet = tjet_vars[tjet_index]

        jetpts.append(newjet['pt'])
        tjetpts.append(tjet['pt'])
        offsets.append(newjet['pt']-tjet['pt'])
        responses.append(newjet['pt']/tjet['pt'])

    return jetpts,tjetpts,offsets,responses

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

jetpts = []
tjetpts = []
offsets = []
responses = []
for event in range(options.start_event,options.end_event):
    e_jetpts,e_tjetpts,e_offsets,e_responses = makejets(options.hs_factor,event,options.weighted)
    jetpts+=e_jetpts
    tjetpts+=e_tjetpts
    offsets+=e_offsets
    responses+=e_responses

if options.weighted: filename = 'weighted_classifier_tjet_vars' 
else: filename = 'binary_classifier_tjet_vars'

print options.submitDir+'/'+filename+'_x'+str(options.hs_factor)+'_jetpts.npy'
numpy.save(options.submitDir+'/'+filename+'_x'+str(options.hs_factor)+'_jetpts.npy',numpy.array(jetpts))
numpy.save(options.submitDir+'/'+filename+'_x'+str(options.hs_factor)+'_tjetpts.npy',numpy.array(tjetpts))
numpy.save(options.submitDir+'/'+filename+'_x'+str(options.hs_factor)+'_offsets.npy',numpy.array(offsets))
numpy.save(options.submitDir+'/'+filename+'_x'+str(options.hs_factor)+'_responses.npy',numpy.array(responses))
