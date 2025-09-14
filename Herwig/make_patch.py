#!/usr/bin/env python3

import sys
import importlib
import subprocess

INPUT = sys.argv[1]
OUTPUT = sys.argv[1] + ".patch"
flist = importlib.import_module(INPUT.rstrip(".py")).flist

xrdfs = subprocess.Popen("xrdfs cms-xrd-global.cern.ch", stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
for f in flist:
    print("ls", f, file=xrdfs.stdin, flush=True)
flist = xrdfs.stdout.read().strip().split()
xrdfs.wait()
with open(OUTPUT, "w") as patch:
    print(f"process.mixData.input.fileNames = cms.untracked.vstring({repr(flist)})", file=patch)
