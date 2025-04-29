exec(open('coffea-reweighting.py').read())  # It's hard to import it...

import uproot
import os
import multiprocessing
import shutil

rootfiles = open('example-files.txt').read().strip().split()
shutil.rmtree('reweighted')
os.makedirs('reweighted')

def process(rootfile):
    try:
        print(f"Processing file: {rootfile}")
        output_filename = os.path.join('reweighted', os.path.basename(rootfile))
        events = NanoEventsFactory.from_root(f'root://cms-xrd-global.cern.ch/{rootfile}', treepath='Events').events()
        output = ReweightProcessor().process(events)
        with uproot.recreate(output_filename) as f: f['tree'] = ak.Array(output)
        print(f"Processed {rootfile} -> {output_filename}")
    except Exception:
        process(rootfile)

multiprocessing.Pool().map(process, rootfiles)
os.system('hadd -f coffea-reweighting-test.root reweighted/*')
