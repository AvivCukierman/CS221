from optparse import OptionParser
import numpy as np
import matplotlib.pyplot as plt

parser = OptionParser()

# job configuration
parser.add_option("--plotDir", help="dir to store the output", default="../Plots")
parser.add_option("--inputDir", help="dir where input is coming from", default="../Output")
parser.add_option("--inputEfficiencies", help="input file containing jet reconstruction efficiences", default="efficiencies.npy")
parser.add_option("--inputTJetPts", help="input file containing truth jet pTs", default="tjetpts.npy")
parser.add_option("--inputMultiples", help="input file containing rates of multiple jet matches", default="multiple_matches.npy")
parser.add_option("--inputOffsets", help="input file containing jet reconstruction offsets", default="offsets.npy")

(options, args) = parser.parse_args()

efficiencies = np.load(options.inputDir+'/'+options.inputEfficiencies)
multiples = np.load(options.inputDir+'/'+options.inputMultiples)
offsets = np.load(options.inputDir+'/'+options.inputOffsets)
tjetpts = np.load(options.inputDir+'/'+options.inputTJetPts)

def fitgaus(data,var,label,legend):
    #data = array(data)
    #m = mean(data)
    #s = std(data)
    #l = data.tolist()
    #data = array([ll for ll in l if ll<m+4*s])
    fig, ax = plt.subplots()
    binsize = 2.0
    n,bins,patches = ax.hist(data,normed=True,bins = (max(data)-min(data))/binsize+1)
    #c = -0.1
    #(a,b,d) = genextreme.fit(data,c)
    #y = genextreme.pdf(bins,a,b,d)
    #plt.plot(bins,y,'r--',linewidth=2)
    ax.set_xlabel(label)
    ax.set_ylabel('a.u.')
    ax.text(0.4, 0.95, legend, fontsize=16, verticalalignment='top', horizontalalignment='left', transform=fig.gca().transAxes)
    ax.set_xlim(-5,40)
    ax.set_ylim(0,.12)
    print var
    fig.savefig(options.plotDir+'/'+var+'_hist.png')
    plt.close()
    return 0

for mintpt in range(10,200,10):
  indices = [i for i,tpt in enumerate(tjetpts) if mintpt<tpt and mintpt+10>tpt]
  legend = str(mintpt)+' GeV < truth jet pT < '+ str(mintpt+10)+ ' GeV'
  if len(indices)>0:
    fitgaus([offsets[index] for index in indices],'offsets'+str(mintpt)+str(mintpt+10),'Offsets (reconstructed jet pT - truth jet pT)',legend)
