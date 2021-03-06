from optparse import OptionParser
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

parser = OptionParser()

# job configuration
parser.add_option("--plotDir", help="dir to store the output", default="../Plots")
parser.add_option("--inputDir", help="dir where input is coming from", default="../Output")
parser.add_option("--jets", help="jetname", default="binary_classifier_tjet_vars")
parser.add_option("--x", help="HS factor", default=1)
#parser.add_option("--inputEfficiencies", help="input file containing jet reconstruction efficiences", default="efficiencies.npy")
parser.add_option("--inputTJetPts", help="input file containing truth jet pTs", default="tjetpts.npy")
#parser.add_option("--inputMultiples", help="input file containing rates of multiple jet matches", default="multiple_matches.npy")
parser.add_option("--inputOffsets", help="input file containing jet reconstruction offsets", default="offsets.npy")
parser.add_option("--inputResponses", help="input file containing jet reconstruction responses", default="responses.npy")

(options, args) = parser.parse_args()

#efficiencies = np.load(options.inputDir+'/'+options.jets+'_x'+str(options.x)+'_'+options.inputEfficiencies)
#multiples = np.load(options.inputDir+'/'+options.jets+'_x'+str(options.x)+'_'+options.inputMultiples)
offsets = np.load(options.inputDir+'/'+options.jets+'_x'+str(options.x)+'_'+options.inputOffsets)
responses = np.load(options.inputDir+'/'+options.jets+'_x'+str(options.x)+'_'+options.inputResponses)
tjetpts = np.load(options.inputDir+'/'+options.jets+'_x'+str(options.x)+'_'+options.inputTJetPts)

def fitgaus(data,var,label,legend):
    #data = array(data)
    #m = mean(data)
    #s = std(data)
    #l = data.tolist()
    #data = array([ll for ll in l if ll<m+4*s])
    fig, ax = plt.subplots()
    if 'offsets' in var: binsize = 2.0
    else: binsize = 0.05
    n,bins,patches = ax.hist(data,normed=True,bins = (max(data)-min(data))/binsize+1)
    if 'offsets' in var: binsize = 0.2
    else: binsize = 0.01
    n,bins,patches = ax.hist(data,normed=True,bins = (max(data)-min(data))/binsize+1, alpha=0)
    gfunc = norm
    (mu,sigma) = gfunc.fit(data)
    y = gfunc.pdf( bins, mu, sigma)
    l = plt.plot(bins, y, 'r--', linewidth=2) 
    ax.set_xlabel(label)
    ax.set_ylabel('a.u.')
    ax.text(0.4, 0.95, legend, fontsize=16, verticalalignment='top', horizontalalignment='left', transform=fig.gca().transAxes)
    if 'offsets' in var:
      ax.set_xlim(-20,20)
      ax.set_ylim(0,.2)
    else:
      ax.set_xlim(0.5,1.5)
      ax.set_ylim(0,7)
    print sigma,var
    fig.savefig(options.plotDir+'/'+var+'_hist.png')
    plt.close()
    return mu,sigma

ptbinsize = 10
ptbins = range(20,200,ptbinsize)
offset_means = {p:0 for p in ptbins}
offset_sigmas = {p:0 for p in ptbins}
response_means = {p:0 for p in ptbins}
response_sigmas = {p:0 for p in ptbins}
for mintpt in ptbins:
#for mintpt in [20]:
  indices = [i for i,tpt in enumerate(tjetpts) if mintpt<tpt and mintpt+ptbinsize>tpt]
  #indices = [i for i,tpt in enumerate(tjetpts) if mintpt<tpt and mintpt+180>tpt]
  legend = str(mintpt)+' GeV < truth jet pT < '+ str(mintpt+ptbinsize)+ ' GeV'
  if len(indices)>0:
    (offset_means[mintpt],offset_sigmas[mintpt])=fitgaus([offsets[index] for index in indices],options.jets+'_'+str(options.x)+'_'+'offsets'+str(mintpt)+str(mintpt+ptbinsize),'Offsets (reconstructed jet pT - truth jet pT)',legend)
    (response_means[mintpt],response_sigmas[mintpt])=fitgaus([responses[index] for index in indices],options.jets+'_'+str(options.x)+'_'+'responses'+str(mintpt)+str(mintpt+ptbinsize),'Response (reconstructed jet pT / truth jet pT)',legend)
    #fitgaus([offsets[index] for index in indices],options.jets+'_'+str(options.x)+'_'+'offsets'+str(mintpt)+str(mintpt+180),'Offsets (reconstructed jet pT - truth jet pT)',legend)
    #fitgaus([responses[index] for index in indices],options.jets+'_'+str(options.x)+'_'+'responses'+str(mintpt)+str(mintpt+180),'Offsets (reconstructed jet pT - truth jet pT)',legend)

import json
with open(options.inputDir+'/'+options.jets+'_ptbins'+str(ptbinsize)+'_x'+str(options.x)+'_offset_means.json','w') as f:
  json.dump(offset_means,f)
with open(options.inputDir+'/'+options.jets+'_ptbins'+str(ptbinsize)+'_x'+str(options.x)+'_offset_sigmas.json','w') as f:
  json.dump(offset_sigmas,f)
with open(options.inputDir+'/'+options.jets+'_ptbins'+str(ptbinsize)+'_x'+str(options.x)+'_response_means.json','w') as f:
  json.dump(response_means,f)
with open(options.inputDir+'/'+options.jets+'_ptbins'+str(ptbinsize)+'_x'+str(options.x)+'_response_sigmas.json','w') as f:
  json.dump(response_sigmas,f)
