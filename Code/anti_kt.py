import math, numpy, copy, sys, random
from Queue import PriorityQueue
from optparse import OptionParser

#Can get help messages with -h or --help
parser = OptionParser()
parser.add_option("-t", "--truth", action="store_true", default=False, help="restrict to only truth = TRUE")
parser.add_option("-e", "--event", default=0, help="which event do we run anti-kt on?", type="int")
parser.add_option("-d", "--debug", action="store_true", default=False, help="debug output")
parser.add_option("-x", "--exclusive", action="store_true", default=False, help="exclusive variant of anti-kt. not recommended")
parser.add_option("-a", "--area", action="store_true", default=False, help="compute active areas")
parser.add_option("-n", "--eta-cutoff", type="float", default=-1, help="If positive, only store jets with |eta| below this cutoff")
#parser.add_option("-g", "--explicit-ghosts", action="store_true", default=False, help="if we compute active areas, include ghosts in output jets and "+
#                                                                                "allow jets consisting entirely of ghosts")
options, args = parser.parse_args()

truth = options.truth
eventNum = options.event
debug = options.debug
exclusive = options.exclusive
inclusive = not exclusive
area = options.area
eta_cutoff = options.eta_cutoff
#explicit_ghosts = options.explicit_ghosts

#if (not area) and explicitGhosts:
#    parser.error("Ghosts are not present if we are not computing jet areas.")

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

if debug: print len(particle_vars)

def modifiedPT(pt):
    return 1.0/pt**2

R = 0.4 if inclusive else 1.0
ptCut = 6.0
dCut = modifiedPT(ptCut)
ghost_area = 0.01 #if this is larger then this program is faster but less accurate. Doesn't seem to work well above 0.01
grid_scatter = 1e-4
pt_scatter = 0.1
mean_ghost_pt = 1e-100
ghost_max_rapidity = 3.0 #maybe increase to 6 if had more computation time

twopi = 2*math.pi

#definitions from FastJet
drap = math.sqrt(ghost_area)
dphi = drap #rectangular grid in eta, phi space. note that eta is basically same as y, since ghosts' masses are small
nphi = int(math.ceil(twopi/dphi))
dphi = twopi/nphi
nrap = int(math.ceil(ghost_max_rapidity/drap))
drap = ghost_max_rapidity / nrap
actual_ghost_area = dphi * drap
n_ghosts = (2*nrap+1)*nphi

if debug: print n_ghosts

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
    return (yI-yJ)**2 + min(phiDiff, twopi - phiDiff)**2

def distIJ(ptI, etaI, phiI, mI, ptJ, etaJ, phiJ, mJ):
    #for now, using eta's instead of y's. maybe not good approximation if m's are non-negligible
#    return min(modifiedPT(ptI), modifiedPT(ptJ)) * deltaIJSquared(etaI, phiI, etaJ, phiJ)/R**2
    return min(modifiedPT(ptI), modifiedPT(ptJ)) * deltaIJSquared(y(etaI, mI, ptI), phiI, y(etaJ, mJ, ptJ), phiJ)/R**2

def antikt_algorithm(event):
    jets_in_event = []
    
    #do anti-kt
    pseudojets = {}
    nextPseudojetIndex = len(event)
    diB = PriorityQueue()
    dIJ = PriorityQueue()
#    diB = {}
#    dIJ = {}
    for i, particle in enumerate(event):
        if truth and not particle['truth']:
            continue
        
        pseudojets[i] = particle
        pseudojets[i]['ghostNumber'] = 0
        if 'index' in particle: pseudojets[i]['particles'] = [particle['index']] #if using a subset of particles, can still match back to original list
        else: pseudojets[i]['particles'] = [i]
        diB.put((modifiedPT(particle['pt']), i))
#        diB[i] = modifiedPT(particle['pt'])
#        dIJ[i] = {}
        for j, otherParticle in enumerate(event):
 #           if j < i:
 #               dIJ[i][j] = dIJ[j][i]
 #           if j == i:
 #               continue            
#            dIJ[i][j] = distIJ(particle['pt'], particle['eta'], particle['phi'], particle['m'],
#                               otherParticle['pt'], otherParticle['eta'], otherParticle['phi'], otherParticle['m'])
            if j < i:
                dIJ.put((distIJ(particle['pt'], particle['eta'], particle['phi'], particle['m'],
                               otherParticle['pt'], otherParticle['eta'], otherParticle['phi'], otherParticle['m']),
                         (i,j)))
    
    if area:
        for irap in range(-nrap, nrap+1):
            for iphi in range(nphi):
                pt = mean_ghost_pt * (1 + random.uniform(0, 1) * pt_scatter)
                phi = (iphi + 0.5) * dphi + dphi * random.uniform(0, 1) * grid_scatter
                if phi < 0:
                    phi += twopi
                if phi > twopi:
                    phi -= twopi
                eta = irap * drap + drap * random.uniform(0, 1) * grid_scatter
                m = 0.0

                diB.put((modifiedPT(pt), nextPseudojetIndex))
    
                for k in pseudojets:
                    pseudojet = pseudojets[k]
                    dIJ.put((distIJ(pt, eta, phi, m,
                                    pseudojet['pt'], pseudojet['eta'], pseudojet['phi'], pseudojet['m']),
                             (k, nextPseudojetIndex)))

                pseudojets[nextPseudojetIndex] = {'ghostNumber': 1, 'particles': [],
                                 'phi': phi,
                                 'eta': eta,
                                 'pt': pt,
                                 'm': m}

                nextPseudojetIndex += 1


    while len(pseudojets) > 0:
        if debug: print len(pseudojets)
#        minDiB, minI = min((diB[i], i) for i in diB)
        
#         minDIJ = None
#         # I < J
#         minIJ = None
        
#        minDIJ, minIJ = min(min((dIJ[i][j], (i,j)) for j in dIJ[i]) for i in dIJ)
        
#         for i in dIJ:
#             for j in dIJ:
#                 if j >= i:
#                     continue
#                 if minDIJ == None or dIJ[i][j] < minDIJ:
#                     minIJ = (i,j)
#                     minDIJ = dIJ[i][j]
        while True:
            if diB.empty():
                minDiB = None
                minI = None
                break
            minDiB, minI = diB.queue[0]
            if minI in pseudojets:
                break
            diB.get()
        while True:
            if dIJ.empty():
                minDIJ = None
                minIJ = None
                break
            minDIJ, minIJ = dIJ.queue[0]
            i, j = minIJ
            if i in pseudojets and j in pseudojets:
                break
            dIJ.get()
        
        #remove pseudojet, form jet if inclusive
        if minDIJ == None or minDiB < minDIJ:
#            if minDIJ != None:
#                dIJ.put((minDIJ, minIJ))

            if (not inclusive) and minDiB > dCut:
                break
            
            diB.get()
            newJet = pseudojets[minI]
            #this second condition makes sure we don't count pseudojets as jets if they're 100% ghosts:
            if inclusive and len(newJet['particles']) > 0 and (eta_cutoff < 0 or abs(newJet['eta']) < eta_cutoff):
                jets_in_event.append({'pt': newJet['pt'], 'eta': newJet['eta'], 'phi': newJet['phi'], 'm': newJet['m'], 
                                      'area': newJet['ghostNumber'] * actual_ghost_area, 'width': None, 
                                      'particle_indices': newJet['particles']})
            del pseudojets[minI]
#            del dIJ[minI]
#            del diB[minI]
#            for j in pseudojets:
#                del dIJ[j][minI]
            
        #merge pseudojets    
        else:
#            if minI != None:
#                diB.put((minDiB, minI))

            if (not inclusive) and minDIJ > dCut:
                break

            dIJ.get()
            i,j = minIJ
            pseudojetI = pseudojets[i]
            pseudojetJ = pseudojets[j]
            ptI = pseudojetI['pt']
            ptJ = pseudojetJ['pt']
            etaI = pseudojetI['eta']
            etaJ = pseudojetJ['eta']
            phiI = pseudojetI['phi']
            phiJ = pseudojetJ['phi']
#            mI = 0
#            mJ = 0
            mI = pseudojetI['m']
            mJ = pseudojetJ['m']

            pI = ptI * math.cosh(etaI)
            pJ = ptJ * math.cosh(etaJ)
            #assume massless or don't
            EI = math.sqrt(pI**2 + mI**2)
            EJ = math.sqrt(pJ**2 + mJ**2)
#            EI = pI
#            EJ = pJ

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
                    newPhi = twopi + math.asin(s)
            else:
                if s >= 0:
                    newPhi = math.acos(c)
                else:
                    newPhi = twopi - math.acos(c)
            
            #different scheme
#            newEtaPrime = (ptI * etaI + ptJ * etaJ) / newPT
#            phiDiff = abs(phiI - phiJ)
#            if phiDiff > math.pi:
#                if phiI < phiJ:
#                    phiI += 2*math.pi
#                else:
#                    phiJ += 2*math.pi
#            newPhiPrime = (ptI * phiI + ptJ * phiJ) / newPT
            
#            newPhiPrime = newPhiPrime if newPhiPrime < 2*math.pi else newPhiPrime - 2*math.pi
#            if debug: print ((newPhi, newPhiPrime), (newEta, newEtaPrime))
            
#            newPT = ptI + ptJ is what I think this merge rule used to be but I forget
#            pIJ = newPT * math.cosh(newEta)            
#            newM = EIJ**2 - pIJ**2
#            newM = math.sqrt(newM) if newM >= 0 else 0

#            diB[nextPseudojetIndex] = modifiedPT(newPT)
            diB.put((modifiedPT(newPT), nextPseudojetIndex))
#            dIJ[nextPseudojetIndex] = {}

            del pseudojets[i]
            del pseudojets[j]
#            del dIJ[i]
#            del dIJ[j]
#            del diB[i]
#            del diB[j]
            for k in pseudojets:
#                del dIJ[k][i]
#                del dIJ[k][j]
                pseudojet = pseudojets[k]
#                dIJ[nextPseudojetIndex][k] = distIJ(newPT, newEta, newPhi, newM,
#                               pseudojet['pt'], pseudojet['eta'], pseudojet['phi'], pseudojet['m'])
#                dIJ[k][nextPseudojetIndex] = dIJ[nextPseudojetIndex][k]
                dIJ.put((distIJ(newPT, newEta, newPhi, newM,
                                pseudojet['pt'], pseudojet['eta'], pseudojet['phi'], pseudojet['m']),
                         (k, nextPseudojetIndex)))
            pseudojets[nextPseudojetIndex] = {'pt': newPT,
                                              'eta': newEta,
                                              'phi': newPhi,
                                              'm': newM,
                                              'ghostNumber': pseudojetI['ghostNumber'] + pseudojetJ['ghostNumber'],
                                              'particles': pseudojetI['particles'] + pseudojetJ['particles']}

            nextPseudojetIndex += 1
            
    if not inclusive:
        for i in pseudojets:
            pseudojet = pseudojets[i]
            if len(pseudojet['particles']) > 0 and (eta_cutoff < 0 or abs(pseudojet['eta']) < eta_cutoff):
                jets_in_event.append(pseudojets[i])
    return jets_in_event

def main():
  #event_vars = numpy.load("../Data/event_vars.npy")
  #jet_vars = numpy.load("../Data/jet_vars.npy")
  particle_vars = numpy.load("../Data/particle_vars.npy")
  #sjet_vars = numpy.load("../Data/sjet_vars.npy")
  #tjet_vars = numpy.load("../Data/tjet_vars.npy")


  num_events = len(particle_vars)
  jets = []
  for i, event in enumerate(particle_vars):
      if i != eventNum:
          continue

      jets_in_event = antikt_algorithm(event)
      jets.append(jets_in_event)

      
  if truth:
      numpy.save("../Data/tjet_vars/our_tjet_vars_"+str(eventNum)+"_new.npy", jets)                
  else:
      numpy.save("../Data/jet_vars/our_jet_vars_"+str(eventNum)+"_new.npy", jets)                

if __name__ == "__main__":
      main()
