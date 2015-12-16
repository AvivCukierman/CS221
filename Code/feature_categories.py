import json

with open('feature_normalization_tjets.json') as f:
  norm = json.load(f)

categories = ['fracjetpt','jetpt','truthchargept','chargept','sumpt','rpt','pt','drjet','alphapuppi','puppi','sigmarho','rho','eta']
feature_categories = {c:{} for c in categories}
for k in norm.keys():
    for c in categories:
        if c in k:
            feature_categories[c][k] = norm[k]
            break 
'''
    if 'fracjetpt' in k:
    elif 'jetpt' in k: pass 
    elif 'truthchargept' in k: pass 
    elif 'chargept' in k: pass 
    elif 'sumpt' in k: pass 
    elif 'rpt' in k: pass 
    elif 'pt' in k: pass 
    elif 'drjet' in k: pass 
    elif 'alphapuppi' in k: pass 
    elif 'puppi' in k: pass 
    elif 'sigmarho' in k: pass
    elif 'rho' in k: pass
    elif 'eta' in k or 'phi' in k: pass 
    else: print k #remaining: truth,index,1
'''

with open('feature_normalization_tjets_categories.json','w') as f:
  json.dump(feature_categories,f)
