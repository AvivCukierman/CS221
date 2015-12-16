from optparse import OptionParser
import numpy as np
from operator import itemgetter

parser = OptionParser()

# job configuration
parser.add_option("--inputDir", help="dir to store the output", default="../Output")
parser.add_option("--inputResults", help="input file containing SGD results", default="SGD_results_e200_i10_x")
parser.add_option("--plotDir", help="dir to store the plots", default="../Plots")
#parser.add_option("-x","--x", help="hard scatter factor", type=int, default=2)

(options, args) = parser.parse_args()

xs = [0.01,0.02,0.05,0.1,0.2,0.5,1.0,2.0,5.0,10.0]
categories = ['fracjetpt','jetpt','truthchargept','chargept','sumpt','rpt','pt','drjet','alphapuppi','puppi','sigmarho','rho','eta']


for c in categories:
  ROC_HS = []
  ROC_PU = []
  ROC_HS_pt = []
  ROC_PU_pt = []
  for x in xs:
    filename = options.inputDir+'/'+options.inputResults+str(x)+'_tjets_'+c+'.txt'

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
  fig.savefig(options.plotDir+'/roc_error_'+c+'.png')
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
  fig.savefig(options.plotDir+'/roc_error_log_'+c+'.png')
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
  fig.savefig(options.plotDir+'/roc_pt_error_'+c+'.png')
###

### pT Error ROC log
  fig, ax = plt.subplots()
  plt.scatter(ROC_HS_pt,ROC_PU_pt,label='Pileup Classification',marker='o',s=100)
  x_offset = []
  for i,(x,HS,PU) in enumerate(zip(xs,ROC_HS_pt,ROC_PU_pt)):
    if x==1.0: plt.annotate('Pileup Factor: '+str(x),xy = (HS, PU), xytext = (-40,20),textcoords = 'offset points', ha = 'right', va = 'bottom',arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=-0.25'))
    elif x==10.0: plt.annotate(str(x),xy = (HS, PU), xytext = (20,20),textcoords = 'offset points', ha = 'right', va = 'bottom',arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=-0.25'))
    else: plt.annotate(str(x),xy = (HS, PU), xytext = (-10, 10),textcoords = 'offset points', ha = 'right', va = 'bottom')
  plt.ylim(0.001,1.0)
  plt.xlim(0.5,1.0)
  plt.xlabel('Fraction of non-pileup pT remaining in jet')
  plt.ylabel('Fraction of pileup pT remaining in jet')
  plt.legend(loc='upper left',frameon=False,numpoints=1)
#plt.gca().invert_yaxis()
#ax.set_xscale('log')
  ax.set_yscale('log')
  fig.savefig(options.plotDir+'/roc_pt_error_log_'+c+'.png')
###
