from optparse import OptionParser
import numpy as np
from operator import itemgetter

parser = OptionParser()

# job configuration
parser.add_option("--inputDir", help="dir to store the output", default="../Output")
#parser.add_option("--inputResults", help="input file containing SGD results", default="SGD_results_e200_i10_x")
parser.add_option("--inputResults", help="input file containing SGD results", default="SGD.o688636.")
#parser.add_option("--inputResults", help="input file containing SGD results", default="gda_results_e200")
parser.add_option("--plotDir", help="dir to store the plots", default="../Plots")
parser.add_option("-n","--n_index", help="index of number of hidden features", type=int, default=0)
parser.add_option("-c","--class_weight", help="use classweight for gda", action="store_true", default=False)

(options, args) = parser.parse_args()

ROC_HS = []
ROC_PU = []
ROC_HS_pt = []
ROC_PU_pt = []
if 'SGD' in options.inputResults:
  xs = [0.01,0.02,0.05,0.1,0.2,0.5,1.0,2.0,5.0,10.0]
  ns = [10,20,30,50,75,100,125,150,175,200]
  n_index = options.n_index
if 'gda' in options.inputResults:
  ms = [0.1,0.5,0.75,1.0,2.0,5.0,10.0,20.0,50.0]
  m = ms[options.n_index]
  ls = [0.01,0.1,0.5,0.75,1.0,2.0,5.0,10.0,20.0,50.0]
  xs = ls

for x_index,x in enumerate(xs):
  if 'SGD' in options.inputResults:
    if 'o688636' in options.inputResults: filename = options.inputDir+'/'+options.inputResults+str(n_index*10+x_index+1)
    else: filename = options.inputDir+'/'+options.inputResults+str(x)+'_tjets.txt'
  if 'gda' in options.inputResults:
    if options.class_weight: filename = options.inputDir + '/' + options.inputResults+'_L'+str(x)+'_m'+str(m)+'_c_tjets.txt'
    else: filename = options.inputDir + '/' + options.inputResults+'_L'+str(x)+'_m'+str(m)+'_tjets.txt'

  with open(filename) as f:
    content = f.readlines()

  HS_effs = []
  PU_effs = []

  HS_pt_effs = []
  PU_pt_effs = []
  for line in content:
    if 'testHSError' in line: 
      eff = line.split(':')[1]
      eff = float(eff)
      HS_effs.append(eff)
    if 'testPUError' in line: 
      eff = line.split(':')[1]
      eff = float(eff)
      PU_effs.append(eff)
    if 'HS retained' in line: 
      eff = line.split(':')[1]
      eff = float(eff)
      HS_pt_effs.append(eff)
    if 'PU retained' in line: 
      eff = line.split(':')[1]
      eff = float(eff)
      PU_pt_effs.append(eff)

  ROC_HS.append(np.mean(HS_effs))
  ROC_PU.append(np.mean(PU_effs))
  ROC_HS_pt.append(np.mean(HS_pt_effs))
  ROC_PU_pt.append(np.mean(PU_pt_effs))

import math
import matplotlib.pyplot as plt
if 'SGD' in options.inputResults:
  if 'o688636' in options.inputResults: n=ns[n_index]
  else: n=1
if 'gda' in options.inputResults:
  pass

### Error ROC
fig, ax = plt.subplots()
plt.scatter(ROC_HS,ROC_PU,label='Pileup Classification',marker='o',s=100)
x_offset = []
for i,(x,HS,PU) in enumerate(zip(xs,ROC_HS,ROC_PU)):
  if x==1.0: plt.annotate('Pileup Factor: '+str(x),xy = (HS, PU), xytext = (80,40),textcoords = 'offset points', ha = 'right', va = 'bottom',arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=-0.25'))
  elif x>0.1: plt.annotate(str(x),xy = (HS, PU), xytext = (20, 10),textcoords = 'offset points', ha = 'right', va = 'bottom')
  else:
    if i%2==0:
      plt.annotate(str(x),xy = (HS, PU), xytext = (-20-20*math.log(x), 10),textcoords = 'offset points', ha = 'right', va = 'bottom',arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=-0.25'))
    else:
      plt.annotate(str(x),xy = (HS, PU), xytext = (-20-20*math.log(x), 10),textcoords = 'offset points', ha = 'right', va = 'bottom',arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0.25'))
plt.ylim(0.0,1.0)
plt.xlim(0.0,1.0)
plt.xlabel('Non-pileup classification error rate')
plt.ylabel('Pileup classification error rate')
plt.legend(loc='upper right',frameon=False,numpoints=1)
plotlabel = 'roc_error'
if 'SGD' in options.inputResults:
  if 'o688636' in options.inputResults: outfilename = options.plotDir+'/'+plotlabel+'_n'+str(n)+'.png'
  else: outfilename = options.plotDir+'/'+plotlabel+'.png'
if 'gda' in options.inputResults:
  if options.class_weight: outfilename = options.plotDir+'/'+plotlabel+'_m'+str(m)+'_c.png'
  else: outfilename = options.plotDir+'/'+plotlabel+'_m'+str(m)+'.png'
fig.savefig(outfilename)
###

### Error ROC log
fig, ax = plt.subplots()
plt.scatter(ROC_HS,ROC_PU,label='Pileup Classification',marker='o',s=100)
x_offset = []
for i,(x,HS,PU) in enumerate(zip(xs,ROC_HS,ROC_PU)):
  if x==1.0: plt.annotate('Pileup Factor: '+str(x),xy = (HS, PU), xytext = (80,40),textcoords = 'offset points', ha = 'right', va = 'bottom',arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=-0.25'))
  else: plt.annotate(str(x),xy = (HS, PU), xytext = (20, 10),textcoords = 'offset points', ha = 'right', va = 'bottom')
plt.ylim(0.0001,1.0)
plt.xlim(0.0,1.0)
plt.xlabel('Non-pileup classification error rate')
plt.ylabel('Pileup classification error rate')
plt.legend(loc='upper right',frameon=False,numpoints=1)
ax.set_yscale('log')
plotlabel = 'roc_error_log'
if 'SGD' in options.inputResults:
  if 'o688636' in options.inputResults: outfilename = options.plotDir+'/'+plotlabel+'_n'+str(n)+'.png'
  else: outfilename = options.plotDir+'/'+plotlabel+'.png'
if 'gda' in options.inputResults:
  if options.class_weight: outfilename = options.plotDir+'/'+plotlabel+'_m'+str(m)+'_c.png'
  else: outfilename = options.plotDir+'/'+plotlabel+'_m'+str(m)+'.png'
fig.savefig(outfilename)
###

### pT Error ROC
fig, ax = plt.subplots()
plt.scatter(ROC_HS_pt,ROC_PU_pt,label='Pileup Classification',marker='o',s=100)
x_offset = []
for i,(x,HS,PU) in enumerate(zip(xs,ROC_HS_pt,ROC_PU_pt)):
  if x==1.0: plt.annotate('Pileup Factor: '+str(x),xy = (HS, PU), xytext = (-40,20),textcoords = 'offset points', ha = 'right', va = 'bottom',arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=-0.25'))
  elif x>0.1: plt.annotate(str(x),xy = (HS, PU), xytext = (-10, 10),textcoords = 'offset points', ha = 'right', va = 'bottom')
  else:
    if i%2==0:
      plt.annotate(str(x),xy = (HS, PU), xytext = (20+20*math.log(x), 10),textcoords = 'offset points', ha = 'right', va = 'bottom',arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0.25'))
    else:
      plt.annotate(str(x),xy = (HS, PU), xytext = (20+20*math.log(x), 10),textcoords = 'offset points', ha = 'right', va = 'bottom',arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=-0.25'))
plt.ylim(0.0,1.0)
plt.xlim(0.0,1.0)
plt.xlabel('Non-pileup pT fraction removed from jet')
plt.ylabel('Pileup pT fraction removed from jet')
plt.legend(loc='upper right',frameon=False,numpoints=1)
plotlabel = 'roc_pt_error'
if 'SGD' in options.inputResults:
  if 'o688636' in options.inputResults: outfilename = options.plotDir+'/'+plotlabel+'_n'+str(n)+'.png'
  else: outfilename = options.plotDir+'/'+plotlabel+'.png'
if 'gda' in options.inputResults:
  if options.class_weight: outfilename = options.plotDir+'/'+plotlabel+'_m'+str(m)+'_c.png'
  else: outfilename = options.plotDir+'/'+plotlabel+'_m'+str(m)+'.png'
fig.savefig(outfilename)
###

### pT Error ROC log
fig, ax = plt.subplots()
if 'SGD' in options.inputResults:
  if 'o688636' in options.inputResults: plt.scatter(ROC_HS_pt,ROC_PU_pt,label='Pileup Classification (neural networks)',marker='o',s=100)
  else: plt.scatter(ROC_HS_pt,ROC_PU_pt,label='Pileup Classification (linear classifier)',marker='o',s=100)
if 'gda' in options.inputResults:
  plt.scatter(ROC_HS_pt,ROC_PU_pt,label='Pileup Classification (gaussian discriminant analysis)',marker='o',s=100)
x_offset = []
for i,(x,HS,PU) in enumerate(zip(xs,ROC_HS_pt,ROC_PU_pt)):
  if 'SGD' in options.inputResults: numberlabel = 'Pileup Factor: '
  if 'gda' in options.inputResults: numberlabel = 'Lambda: '
  if x==1.0: plt.annotate(numberlabel+str(x),xy = (HS, PU), xytext = (-40,20),textcoords = 'offset points', ha = 'right', va = 'bottom',arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=-0.25'))
  elif x==10.0: plt.annotate(str(x),xy = (HS, PU), xytext = (20,20),textcoords = 'offset points', ha = 'right', va = 'bottom',arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=-0.25'))
  else: plt.annotate(str(x),xy = (HS, PU), xytext = (-10, 10),textcoords = 'offset points', ha = 'right', va = 'bottom')
if 'SGD' in options.inputResults:
  plt.ylim(0.001,1.0)
  plt.xlim(0.5,1.0)
if 'gda' in options.inputResults:
  plt.ylim(0.1,1.0)
  plt.xlim(0.9,1.0)
plt.xlabel('Fraction of non-pileup pT remaining in jet')
plt.ylabel('Fraction of pileup pT remaining in jet')
plt.legend(loc='upper left',frameon=False,numpoints=1)
#plt.gca().invert_yaxis()
#ax.set_xscale('log')
ax.set_yscale('log')
plotlabel = 'roc_pt_error_log'
if 'SGD' in options.inputResults:
  if 'o688636' in options.inputResults: outfilename = options.plotDir+'/'+plotlabel+'_n'+str(n)+'.png'
  else: outfilename = options.plotDir+'/'+plotlabel+'.png'
if 'gda' in options.inputResults:
  if options.class_weight: outfilename = options.plotDir+'/'+plotlabel+'_m'+str(m)+'_c.png'
  else: outfilename = options.plotDir+'/'+plotlabel+'_m'+str(m)+'.png'
fig.savefig(outfilename)
###
