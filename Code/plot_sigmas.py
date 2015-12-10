from optparse import OptionParser
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

parser = OptionParser()

# job configuration
parser.add_option("--plotDir", help="dir to store the output", default="../Plots")
parser.add_option("--inputDir", help="dir where input is coming from", default="../Output")
parser.add_option("--jets", help="which jets", default="weighted_classifier_tjet_vars")
parser.add_option("-r","--response_correction", help="correct for response by divigin", action="store_true",default=False)

(options, args) = parser.parse_args()

def plot_sigma(means,sigmas,ls,l,c,fig,ax):
  x = sigmas.keys()
  x = [int(xx) for xx in x if xx!='10']
  x.sort()
  x = [str(xx) for xx in x]
  if options.response_correction: y = [sigmas[xx]/means[xx] for xx in x]
  else: y = [sigmas[xx] for xx in x]
  x = [str(int(xx)+5) for xx in x]
  plt.plot(x,y,marker='o',linestyle=ls,label=l,color=c)

  return 0 

def plot_mean(means,ls,l,c,fig,ax):
  x = means.keys()
  x = [int(xx) for xx in x if xx!='10']
  x.sort()
  x = [str(xx) for xx in x]
  y = [means[xx] for xx in x]
  x = [str(int(xx)+5) for xx in x]
  plt.plot(x,y,marker='o',linestyle=ls,label=l,color=c)

  return 0 

import json
xs = [0.01,0.02,0.05,0.1,0.2,0.5,1.0,2.0,5.0,10.0,0.0]
cs = ['r','r','orange','orange','y','y','g','g','b','b','black']
def make_plot(var,jets):
  fig, ax = plt.subplots()
  for x,c in zip(xs,cs):
    if x>0:
      filename = options.inputDir+'/'+jets+'_x'+str(x)
      with open(filename+'_response_means'+'.json') as f:
        data_means = json.load(f)
      with open(filename+'_'+var+'_sigmas'+'.json') as f:
        data_sigmas = json.load(f)
      label = 'Pileup factor: '+str(x)
      plot_sigma(data_means,data_sigmas,'-',label,c,fig,ax)
    else:
      for b in range(2):
        if b==0: filename = options.inputDir+'/'+'binary_classifier_tjet_vars'+'_x'+str(x)
        if b==1: filename = options.inputDir+'/'+'weighted_classifier_tjet_vars'+'_x'+str(x)
        with open(filename+'_response_means'+'.json') as f:
          data_means = json.load(f)
        with open(filename+'_'+var+'_sigmas'+'.json') as f:
          data_sigmas = json.load(f)
        if b==0: label = 'Uncorrected jets'
        if b==1: label = 'Truth jets'
        if b==0: ls = '-'
        if b==1: ls = '--'
        plot_sigma(data_means,data_sigmas,ls,label,c,fig,ax)

  plt.xlabel('pT of Truth Jets')
  if options.response_correction: plt.ylabel('Width of '+var+' / Average response')
  else: plt.ylabel('Width of '+var)
  plt.legend(loc='upper right',frameon=False,numpoints=1)
  if options.response_correction: fig.savefig(options.plotDir+'/'+jets+'_'+var+'_response_corr.png')
  else: fig.savefig(options.plotDir+'/'+jets+'_'+var+'.png')
  plt.close()

  fig, ax = plt.subplots()
  for x,c in zip(xs,cs):
    if x>0:
      filename = options.inputDir+'/'+jets+'_x'+str(x)
      with open(filename+'_'+var+'_means'+'.json') as f:
        data_means = json.load(f)
      label = 'Pileup factor: '+str(x)
      plot_mean(data_means,'-',label,c,fig,ax)
    else:
      for b in range(2):
        if b==0: filename = options.inputDir+'/'+'binary_classifier_tjet_vars'+'_x'+str(x)
        if b==1: filename = options.inputDir+'/'+'weighted_classifier_tjet_vars'+'_x'+str(x)
        with open(filename+'_'+var+'_means'+'.json') as f:
          data_means = json.load(f)
        if b==0: label = 'Uncorrected jets'
        if b==1: label = 'Truth jets'
        if b==0: ls = '-'
        if b==1: ls = '--'
        plot_mean(data_means,ls,label,c,fig,ax)

  plt.xlabel('pT of Truth Jets')
  plt.ylabel('Mean of '+var)
  plt.legend(loc='upper right',frameon=False,numpoints=1)
  fig.savefig(options.plotDir+'/'+jets+'_'+var+'_means.png')
  plt.close()

vars = ['response','offset']
for var in vars:
  make_plot(var,options.jets)
