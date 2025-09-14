#!/usr/bin/env python3

import os
import sys
import importlib

INPUT = sys.argv[1]
OUTPUT = sys.argv[1] + ".patch"
flist = importlib.import_module(INPUT.rstrip(".py")).flist
prefix = flist[0]
while True:
    for f in flist:
        if not f.startswith(prefix):
            prefix = os.path.dirname(prefix)
            break
    else:
        break
with os.popen("xrdfs root://eoscms.cern.ch ls -R " + prefix) as p:
    fset = set(p.read().strip().split())
flist = [f for f in flist if f in fset]
with open(OUTPUT, "w") as patch:
    print(f"process.mixData.input.fileNames = cms.untracked.vstring({repr(flist)})", file=patch)
