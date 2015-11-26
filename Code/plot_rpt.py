from optparse import OptionParser
import numpy as np
import matplotlib.pyplot as plt

parser = OptionParser()

# job configuration
parser.add_option("--plotDir", help="dir to store the output", default="plots")
parser.add_option("--inputDir", help="dir where input is coming from", default="output")
parser.add_option("--inputChargedRpts", help="input file containing charged rpts", default="charged_rpts.npy")
parser.add_option("--inputAllRpts", help="input file containing charged+uncharged rpts", default="all_rpts.npy")

(options, args) = parser.parse_args()

charged_rpts = np.load(options.inputDir+'/'+options.inputChargedRpts)
all_rpts = np.load(options.inputDir+'/'+options.inputAllRpts)

def fitgaus(data,var,label):
  #data = array(data)
  #m = mean(data)
  #s = std(data)
  #l = data.tolist()
  #data = array([ll for ll in l if ll<m+4*s])
  fig, ax = plt.subplots()
  n,bins,patches = ax.hist(data,normed=True,bins=10)
  #c = -0.1
  #(a,b,d) = genextreme.fit(data,c)
  #y = genextreme.pdf(bins,a,b,d)
  #plt.plot(bins,y,'r--',linewidth=2)
  ax.set_xlabel(label)
  ax.set_ylabel('a.u.')
  print var
  fig.savefig(options.plotDir+'/'+var+'_hist.png')
  plt.close()
  return 0

def scatter(x,y,var,xlabel,ylabel):
  fig,ax = plt.subplots()
  ax.scatter(x,y)
  ax.set_xlabel(xlabel)
  ax.set_ylabel(ylabel)
  print var
  fig.savefig(options.plotDir+'/'+var+'_scatter.png')
  plt.close()
  return 0
  

fitgaus(charged_rpts,'charged_rpts','Fraction of pT from HS PV (only charged)')
fitgaus(all_rpts,'all_rpts','Fraction of pT from HS PV (charged and uncharged)')
scatter(charged_rpts,all_rpts,'rpts','Fraction of pT from HS PV (only charged)','Fraction of pT from HS PV (charged and uncharged)')
