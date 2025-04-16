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
with open('powheg-madspin-datasets.txt', 'w') as file:
    print(*dsid, sep='\n', file=file)

miss = open('powheg-madspin.miss', 'w')
gridpacks = set()
for dataset, prepid in sorted(dsid.items()):
    print(dataset, prepid)
    # https://cms-pdmv-prod.web.cern.ch/mcm/public/restapi/requests/get_fragment/BTV-RunIII2024Summer24wmLHEGS-00004
    text = requests.get(f'https://cms-pdmv-prod.web.cern.ch/mcm/public/restapi/requests/get_fragment/{prepid[0]}').text
    # '/cvmfs/cms.cern.ch/phys_generator/gridpacks/slc6_amd64_gcc481/13TeV/powheg/V2/ST_tch_4f_13TeV_antitop_powheg_madspin/st_tch_4f_ckm_NLO_antitop_powheg_madspin_tarball.tar.xz'
    try:
        gridpack = re.search(r"""/cvmfs/.*\.(?:tar\.xz|tgz)""", text).group()
        print(gridpack)
        gridpacks.add(gridpack)
    except Exception:
        print('-' * 80 + '\n', text, '-' * 80 + '\n', file=miss)

with open('powheg-madspin-gridpacks.txt', 'w') as file:
    print(*gridpacks, sep='\n', file=file)
