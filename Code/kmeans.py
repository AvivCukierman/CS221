import math, numpy, copy, sys
from optparse import OptionParser
import operator 

parser = OptionParser()
parser.add_option("-t", "--truth", action="store_true", default=False, help="restrict to only truth = TRUE")
parser.add_option("-e", "--event", default=0)
parser.add_option("-d", "--debug", action="store_true", default=False)
options, args = parser.parse_args()

truth = options.truth
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

def modifiedPT(pt):
    return 1.0/pt**2

R = 0.4

def y(eta, m, pt):
    #pz is signed, pL unsigned
    pz = pt * math.sinh(eta)
    pL = abs(pz)
    p = pt * math.cosh(eta)
    energy = math.sqrt(p**2 + m**2)
#    return 0.5 * math.log((energy + pz) / (energy - pz))
    #computational approximation to above line:
    return (0.5 * math.log((pt**2 + m**2) / (energy + pL)**2)) * (-1 if pz > 0 else 1)

def deltaIJSquared(yI, phiI, yJ, phiJ):
    phiDiff = abs(phiI - phiJ)
    return (yI-yJ)**2 + min(phiDiff, 2*math.pi - phiDiff)**2

def distIJ(ptI, etaI, phiI, mI, ptJ, etaJ, phiJ, mJ):
    #for now, using eta's instead of y's. maybe not good approximation if m's are non-negligible
#    return min(modifiedPT(ptI), modifiedPT(ptJ)) * deltaIJSquared(etaI, phiI, etaJ, phiJ)/R**2
    return min(modifiedPT(ptI), modifiedPT(ptJ)) * deltaIJSquared(y(etaI, mI, ptI), phiI, y(etaJ, mJ, ptJ), phiJ)/R**2

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

  return newPT,newEta,newPhi,newM

import pdb
num_events = len(particle_vars)
jets = []
for i, event in enumerate(particle_vars):
    if i != eventNum:
        continue
    if not truth: jets = numpy.load("../Data/jet_vars/our_jet_vars_"+str(i)+".npy")[0] #use antikt jets as wet start
    else: jets = numpy.load("../Data/tjet_vars/our_tjet_vars_"+str(i)+".npy")[0] #use antikt jets as wet start

    for it in range(5):
      #E? step
      newjets = []
      for jet in jets:
        if len(jet['particle_indices'])>0:
          newjet_pt,newjet_eta,newjet_phi,newjet_m = combine_all([particle for i,particle in enumerate(event) if i in jet['particle_indices']])
          newjet = {'pt':newjet_pt,'eta':newjet_eta,'phi':newjet_phi,'m':newjet_m,'particle_indices':[]}
          newjets.append(newjet)
      jets = newjets
      print jets[0:2]

      #M? step
      for i, particle in enumerate(event):
          if truth and not particle['truth']:
              continue

          dists = [deltaIJSquared(y(particle['eta'],particle['m'],particle['pt']),particle['phi'],y(jet['eta'],jet['m'],jet['pt']),jet['phi']) for jet in jets]
          min_index = min(enumerate(dists), key=operator.itemgetter(1))[0]
          jets[min_index]['particle_indices'].append(i)



    
'''if truth:
    numpy.save("../Data/tjet_vars/our_tjet_vars_"+str(eventNum)+".npy", jets)                
else:
    numpy.save("../Data/jet_vars/our_jet_vars_"+str(eventNum)+".npy", jets) 
'''
