import math, numpy, copy, sys
from optparse import OptionParser
import operator 

parser = OptionParser()
parser.add_option("-e", "--event", help='Which event to do', default=0)
parser.add_option("-d", "--debug", action="store_true", default=False)
options, args = parser.parse_args()

eventNum = int(options.event)
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

def y(eta, m, pt):
  #pz is signed, pL unsigned
  pz = pt * math.sinh(eta)
  pL = abs(pz)
  p = pt * math.cosh(eta)
  energy = math.sqrt(p**2 + m**2)
#    return 0.5 * math.log((energy + pz) / (energy - pz))
  #computational approximation to above line:
  return (0.5 * math.log((pt**2 + m**2) / (energy + pL)**2)) * (-1 if pz > 0 else 1)

def deltaIJSquared(p1,p2):
  return deltaY(p1,p2)**2 + deltaPhi(p1,p2)**2 

def deltaPhi(p1,p2):
  phiI = p1['phi']
  phiJ = p2['phi']
  phiDiff = abs(phiI - phiJ)
  return min(phiDiff, 2*math.pi - phiDiff)**2

def deltaY(p1,p2):
  YI = p1['y']
  YJ = p2['y']
  return abs(YI-YJ)

for i, particles in enumerate(particle_vars):
  if not i==eventNum: continue
  for particle in particles:
    particle['y']=y(particle['eta'],particle['m'],particle['pt']) #saves calculation later

import pdb
for i, particles in enumerate(particle_vars):
  if not i==eventNum: continue
  jets = numpy.load("../Data/jet_vars/our_jet_vars_"+str(i)+".npy")[0] #only look at particles within high pT jets
  newparticles = []
  for jet in jets:
    #print jet
    if jet['pt']<20: continue
    if len(jet['particle_indices'])==0: continue
    for particle_index in jet['particle_indices']:
      particle = particles[particle_index]
      newparticle = {'index':particle_index}

      newparticle['eta'] = particle['eta']
      newparticle['phi'] = particle['phi']

      jet['y']=y(jet['eta'],jet['m'],jet['pt'])
      newparticle['jetpt']=jet['pt']
      newparticle['drjet']=math.sqrt(deltaIJSquared(jet,particle))

      for rmin in numpy.linspace(0.1,0.4,4): 

        sumpt=0
        sumchargept = 0
        sumtruthchargept = 0
        for p in particles:
          ydiff = deltaY(p,particle)
          if ydiff>rmin: continue
          phidiff = deltaPhi(p,particle)
          if phidiff>rmin: continue
          if phidiff**2+ydiff**2<rmin**2:
            pt = p['pt']
            sumpt+=pt
            if p['charge']==1:
              sumchargept+=pt
              if p['truth']==1:
                sumtruthchargept+=pt
        sumptfeature = 'sumpt'+str(int(rmin*10))
        newparticle[sumptfeature] = sumpt
        sumchargeptfeature = 'sumchargept'+str(int(rmin*10))
        newparticle[sumchargeptfeature] = sumchargept
        sumtruthchargeptfeature = 'sumtruthchargept'+str(int(rmin*10))
        newparticle[sumtruthchargeptfeature] = sumtruthchargept
        for ptmin in numpy.linspace(0,5,6):
          newparticle[sumptfeature+'gt'+str(ptmin)] = int(sumpt>ptmin)
          newparticle[sumchargeptfeature+'gt'+str(ptmin)] = int(sumchargept>ptmin)
          newparticle[sumtruthchargeptfeature+'gt'+str(ptmin)] = int(sumtruthchargept>ptmin)
        rpt = sumchargept/sumtruthchargept if sumtruthchargept>0 else -1 
        rptfeature = 'rpt'+str(int(rmin*10))
        for rptmin in numpy.linspace(0,0.8,5):
          newparticle[rptfeature+'gt'+str(rptmin)] = int(rpt>rptmin)

        for rmax in numpy.arange(rmin+0.1,0.41,0.1):
          puppi=0
          for i,p in enumerate(particles):
            ydiff = deltaY(p,particle)
            if ydiff>rmax: continue
            phidiff = deltaPhi(p,particle)
            if phidiff>rmax: continue
            dr = math.sqrt(phidiff**2+ydiff**2)
            if dr>rmin and dr<rmax:
              puppi+=p['pt']/dr
              print dr,p['pt'],rmin,rmax,i
          featurename = 'puppi'+str(int(rmin*10))+str(int(rmax*10))
          newparticle[featurename] = puppi
          alphapuppi = math.log(puppi) if puppi>0 else -float('Inf')
          for alphamin in numpy.linspace(0,10,6):
            newparticle['alpha'+featurename+'gt'+str(int(alphamin))] = int(alphapuppi>alphamin)

      newparticles.append(newparticle)
      features = newparticle.keys()
      features.sort()
      for k in features: print k,newparticle[k]
      break

#print newparticles 

#numpy.save("../Data/particle_features/particle_features_"+str(eventNum)+".npy", newparticles) 
