from optparse import OptionParser
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

parser = OptionParser()

# job configuration
parser.add_option("--plotDir", help="dir to store the output", default="../Plots")
parser.add_option("--inputDir", help="dir where input is coming from", default="../Output")
parser.add_option("--jets", help="jetname", default="binary_classifier_jet_vars")
parser.add_option("--x", help="HS factor", default=1)
parser.add_option("--inputFakeRate", help="input file containing jet reconstruction fake rates", default="fake_rates.npy")
(options, args) = parser.parse_args()


jets = ['binary_classifier_jet_vars','weighted_classifier_jet_vars']
cs = ['r','b']
labels = ['Binary classifier','Weighted classifier']
fig, ax = plt.subplots()
for jet,c,label in zip(jets,cs,labels):
  plot_fake_rates = []
  xs = [2,5,10,20,40,80]
  for x in xs:
    fake_rates = np.load(options.inputDir+'/'+jet+'_x'+str(x)+'_'+options.inputFakeRate)
    plot_fake_rates.append(np.mean(fake_rates))

  import matplotlib.pyplot as plt
  plt.scatter(xs,plot_fake_rates,marker='o',s=40,c=c,label=label)
  ax.set_xscale('log')
  plt.ylim(0.0,1.0)
  plt.xlim(1,100)
  plt.xlabel('Pileup factor')
  plt.ylabel('Pileup jet rate')
plt.legend(loc='upper right',frameon=False,numpoints=1)
fig.savefig(options.plotDir+'/fake_rate.png')
