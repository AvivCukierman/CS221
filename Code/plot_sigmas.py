from optparse import OptionParser
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

parser = OptionParser()

# job configuration
parser.add_option("--plotDir", help="dir to store the output", default="../Plots")
parser.add_option("--inputDir", help="dir where input is coming from", default="../Output")

(options, args) = parser.parse_args()

def plot_sigma(data,ls,l,c,fig,ax):
  x = data.keys()
  x = [int(xx) for xx in x]
  x.sort()
  x = [str(xx) for xx in x]
  y = [data[xx] for xx in x]
  plt.plot(x,y,marker='o',linestyle=ls,label=l,color=c)

  return 0 

import json
xs = [2,5,10,20,40,80]
cs = ['r','orange','y','g','b','purple']
def make_plot(var):
  fig, ax = plt.subplots()
  for x,c in zip(xs,cs):
    with open(options.inputDir+'/binary_classifier_jet_vars_x'+str(x)+'_'+var+'.json') as f:
      data = json.load(f)
    plot_sigma(data,'-','Pileup factor: '+str(x),c,fig,ax)
  for x,c in zip(xs,cs):
    with open(options.inputDir+'/weighted_classifier_jet_vars_x'+str(x)+'_'+var+'.json') as f:
      data = json.load(f)
    plot_sigma(data,'--','Pileup factor: '+str(x),c,fig,ax)
  with open(options.inputDir+'/jet_vars_x'+str(1)+'_'+var+'.json') as f:
    data = json.load(f)
  plot_sigma(data,'-','Uncorrected Jets','black',fig,ax)
  plt.xlabel('pT of truth particles')
  if var== 'response_sigmas':plt.ylabel('Width of response')
  if var== 'response_means':plt.ylabel('Mean of response')
  if var== 'offset_sigmas':plt.ylabel('Width of offset')
  if var== 'offset_means':plt.ylabel('Mean of offset')
  plt.legend(loc='upper right',frameon=False,numpoints=1)
  fig.savefig(options.plotDir+'/'+var+'.png')
  plt.close()

vars = ['response_sigmas','response_means','offset_sigmas','offset_means']
for var in vars:
  make_plot(var)
