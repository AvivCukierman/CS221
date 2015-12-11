from optparse import OptionParser
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

parser = OptionParser()

# job configuration
parser.add_option("--plotDir", help="dir to store the output", default="../Plots")
parser.add_option("--inputDir", help="dir where input is coming from", default="../Output")
parser.add_option("--jets", help="which jets", default="weighted_classifier_tjet_vars")
parser.add_option("-p","--ptbinsize", help="pT bin size", default=10)
parser.add_option("-r","--response_correction", help="correct for response by dividing", action="store_true",default=False)

(options, args) = parser.parse_args()

def plot_sigma(means,sigmas,ls,l,c,fig,ax):
  x = sigmas.keys()
  x = [int(xx) for xx in x if xx!='10']
  x.sort()
  x = [str(xx) for xx in x]
  if options.response_correction: y = [sigmas[xx]/means[xx] for xx in x]
  else: y = [sigmas[xx] for xx in x]
  x = [str(int(xx)+int(options.ptbinsize)/2) for xx in x]
  plt.plot(x,y,marker='o',linestyle=ls,label=l,color=c)

  return 0 

def plot_mean(means,ls,l,c,fig,ax):
  x = means.keys()
  x = [int(xx) for xx in x if xx!='10']
  x.sort()
  x = [str(xx) for xx in x]
  y = [means[xx] for xx in x]
  x = [str(int(xx)+int(options.ptbinsize)/2) for xx in x]
  plt.plot(x,y,marker='o',linestyle=ls,label=l,color=c)

  return 0 

import json
def make_plot(var,xs,cs,labels,linestyles,jets,extra):
  fig, ax = plt.subplots()
  for x,c,label,linestyle,jet in zip(xs,cs,labels,linestyles,jets):
    filename = options.inputDir+'/'+jet+'_ptbins'+str(options.ptbinsize)+'_x'+str(x)
    with open(filename+'_response_means'+'.json') as f:
      data_means = json.load(f)
    with open(filename+'_'+var+'_sigmas'+'.json') as f:
      data_sigmas = json.load(f)
    plot_sigma(data_means,data_sigmas,linestyle,label,c,fig,ax)

  plt.xlabel('pT of Truth Jets (GeV)')
  if options.response_correction: plt.ylabel('Width of '+var+' / Average response')
  else: plt.ylabel('Width of '+var+' (Gev)')
  plt.ylim(0,10)
  plt.legend(loc='upper right',frameon=False,numpoints=1)
  if options.response_correction: fig.savefig(options.plotDir+'/'+options.jets+'_'+var+'_response_corr'+extra+'.png')
  else: fig.savefig(options.plotDir+'/'+options.jets+'_'+var+extra+'.png')
  plt.close()

  fig, ax = plt.subplots()
  for x,c,label,linestyle,jet in zip(xs,cs,labels,linestyles,jets):
    filename = options.inputDir+'/'+jet+'_ptbins'+str(options.ptbinsize)+'_x'+str(x)
    with open(filename+'_'+var+'_means'+'.json') as f:
      data_means = json.load(f)
    plot_mean(data_means,linestyle,label,c,fig,ax)
  plt.xlabel('pT of Truth Jets (GeV)')
  plt.ylabel('Mean of '+var+' (GeV)')
  plt.ylim(-5,25)
  plt.legend(loc='upper right',frameon=False,numpoints=1)
  fig.savefig(options.plotDir+'/'+options.jets+'_'+var+'_means'+extra+'.png')
  plt.close()


extra = '_all' 
extra+='ptbins'+str(options.ptbinsize)
###
xs = [0.01,0.02,0.05,0.1,0.2,0.5,1.0,2.0,5.0,10.0]
cs = ['r','r','orange','orange','y','y','g','g','b','b']
labels = ['Pileup factor: '+str(x) for x in xs]
linestyles = ['-']*len(xs)
jets = [options.jets]*len(xs)

xs += [0.0]*3
cs += ['black']*3
labels += ['Uncorrected','"Truth" classifier','Rho-area subtracted']
linestyles += ['-','-.','--']
jets += ['binary_classifier_tjet_vars','weighted_classifier_tjet_vars','binary_classifier_tjet_vars_sub']

vars = ['response','offset']
for var in vars:
  make_plot(var,xs,cs,labels,linestyles,jets,extra)
###

extra = '_best' 
extra+='ptbins'+str(options.ptbinsize)
###
xs = [1.0,0.0,0.0,0.0]
cs = ['g','b','black','r']
labels = ['Pileup Classification','Uncorrected','"Truth" classifier','Rho-area subtracted']
linestyles = ['-','-','--','-']
jets = [options.jets,'binary_classifier_tjet_vars','weighted_classifier_tjet_vars','binary_classifier_tjet_vars_sub']

vars = ['response','offset']
for var in vars:
  make_plot(var,xs,cs,labels,linestyles,jets,extra)
###
