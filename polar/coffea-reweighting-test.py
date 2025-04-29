exec(open('coffea-reweighting.py').read())  # It's hard to import it...

import uproot
import os
import multiprocessing

rootfiles = open('example-files.txt').read().strip().split()

def process(rootfile):
    print(f"Processing file: {rootfile}")
    os.makedirs('reweighted', exist_ok=True)
    output_filename = os.path.join('reweighted', os.path.basename(rootfile))
    events = NanoEventsFactory.from_root(f'root://cms-xrd-global.cern.ch/{rootfile}', treepath='Events').events()
    output = ReweightProcessor().process(events)
    with uproot.recreate(output_filename) as f: f['tree'] = ak.Array(output)
    print(f"Processed {rootfile} -> {output_filename}")

multiprocessing.Pool().map(process, rootfiles)
