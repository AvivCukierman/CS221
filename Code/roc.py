from optparse import OptionParser
import numpy as np
from operator import itemgetter

parser = OptionParser()

# job configuration
parser.add_option("--inputDir", help="dir to store the output", default="../Output")
parser.add_option("--inputResults", help="input file containing SGD results", default="SGD_results_e80_i10_x")
parser.add_option("--plotDir", help="dir to store the plots", default="../Plots")
#parser.add_option("-x","--x", help="hard scatter factor", type=int, default=2)

(options, args) = parser.parse_args()

ROC_HS = []
ROC_PU = []
xs = [2,5,10,20,40,80]
for x in xs:
  filename = options.inputDir+'/'+options.inputResults+str(x)+'.txt'

  with open(filename) as f:
    content = f.readlines()

  HS_effs = []
  PU_effs = []
  for line in content:
    if 'testHSError' in line: 
      eff = line.split(':')[1]
      eff = float(eff)
      HS_effs.append(eff)
    if 'testPUError' in line: 
      eff = line.split(':')[1]
      eff = float(eff)
      PU_effs.append(eff)

  ROC_HS.append(np.mean(HS_effs))
  ROC_PU.append(np.mean(PU_effs))

import matplotlib.pyplot as plt
fig, ax = plt.subplots()
plt.scatter(ROC_HS,ROC_PU,marker='o',s=100)
plt.ylim(0.6,1.0)
plt.xlim(0.0,.12)
plt.xlabel('Non-pileup tagging error rate')
plt.ylabel('Pileup tagging error rate')
fig.savefig(options.plotDir+'/roc.png')
