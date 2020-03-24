# Use python 3

import os, glob

source = 'ZINC000002015152'
F = open(f'MPro_{source}.yaml','r')
dat = F.read()
F.close()

sdf_FNs = glob.glob('../ligands/0-build/*.sdf')
sdf_FNs.remove(f'../ligands/0-build/{source}.sdf')

for sdf_FN in sdf_FNs:
  destination = os.path.basename(sdf_FN).split('.')[0]
  # if not os.path.isfile(f'MPro_{destination}.yaml'):
  print(f'Writing MPro_{destination}.yaml')
  F = open(f'MPro_{destination}.yaml','w')
  F.write(dat.replace(source, destination))
  F.close()
