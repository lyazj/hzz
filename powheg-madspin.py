#!/usr/bin/env python3

import re
import requests
import json

data = json.load(open('powheg-madspin.json'))
results = data['results']
dsid = { }
for result in results:
    name = result['dataset_name']
    dsid[name] = [*dsid.get(name, []), result['prepid']]

# https://cms-pdmv-prod.web.cern.ch/mcm/public/restapi/requests/get_fragment/BTV-RunIII2024Summer24wmLHEGS-00004
gridpacks = set()
for dataset, prepid in sorted(dsid.items()):
    print(dataset, prepid)
    res = requests.get(f'https://cms-pdmv-prod.web.cern.ch/mcm/public/restapi/requests/get_fragment/{prepid[0]}')
    # '/cvmfs/cms.cern.ch/phys_generator/gridpacks/slc6_amd64_gcc481/13TeV/powheg/V2/ST_tch_4f_13TeV_antitop_powheg_madspin/st_tch_4f_ckm_NLO_antitop_powheg_madspin_tarball.tar.xz'
    gridpack = re.search(r"""'/cvmfs/.*\.tar\.xz'""", res.text).group()
    print(gridpack)
    gridpacks.add(gridpack)

with open('powheg-madspin.txt', 'w') as file:
    print(*gridpacks, sep='\n', file=file)
