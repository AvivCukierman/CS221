import numpy

event = 1

jet_vars = numpy.load("../Data/jet_vars/jet_vars.npy")
our_jet_vars = numpy.load("../Data/jet_vars/our_jet_vars_" + str(event) + ".npy")

tjet_vars = numpy.load("../Data/tjet_vars/tjet_vars.npy")
our_tjet_vars = numpy.load("../Data/tjet_vars/our_tjet_vars_" + str(event) + ".npy")

our_stuffs = sorted([(our_jet_vars[0][j]['pt'], j) for j in range(len(our_jet_vars[0]))])[::-1]
stuffs = sorted([(jet_vars[event][j]['pt'], j) for j in range(len(jet_vars[event]))])[::-1]
our_stuff = [our_jet_vars[0][j] for _,j in our_stuffs]
stuff = [jet_vars[event][j] for _,j in stuffs]

print [thing['pt'] for thing in stuff if abs(thing['eta']) < 2]
print [thing['pt'] for thing in our_stuff if abs(thing['eta']) < 2]

print [(thing['eta'],thing['phi']) for thing in stuff if abs(thing['eta']) < 2]
print [(thing['eta'],thing['phi']) for thing in our_stuff if abs(thing['eta']) < 2]

print [thing['m'] for thing in stuff if abs(thing['eta']) < 2]
print [thing['m'] for thing in our_stuff if abs(thing['eta']) < 2]

print [thing['area'] for thing in stuff if abs(thing['eta'])  < 2]
print [thing['area'] for thing in our_stuff if abs(thing['eta'])  < 2]

print [len(thing['particle_indices']) for thing in our_stuff if abs(thing['eta']) < 2]

print 'Things not matched:'
for thing in stuff:
  matched = False
  for our_thing in our_stuff:
    if abs(thing['pt']-our_thing['pt'])<0.1 and abs(thing['eta']-our_thing['eta'])<0.1 and abs(thing['phi']-our_thing['phi'])<0.1: matched = True
  if not matched:
    print thing['pt'],thing['eta'],thing['phi']

print 'Our things not matched:'
for thing in our_stuff:
  if abs(thing['eta'])>2: continue
  if abs(thing['pt'])<5: continue
  matched = False
  for our_thing in stuff:
    if abs(thing['pt']-our_thing['pt'])<0.1 and abs(thing['eta']-our_thing['eta'])<0.1 and abs(thing['phi']-our_thing['phi'])<0.1: matched = True
  if not matched:
    print thing['pt'],thing['eta'],thing['phi']
