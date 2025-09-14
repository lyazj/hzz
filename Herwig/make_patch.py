import sys
import importlib
import multiprocessing

INPUT = sys.argv[1]
OUTPUT = sys.argv[1] + ".patch"
flist = importlib.import_module(INPUT.rstrip(".py")).flist


def is_valid(fpath):
    import ROOT

    file = ROOT.TFile.Open("root://cms-xrd-global.cern.ch/" + fpath)
    return bool(file)


pool = multiprocessing.Pool()
vlist = pool.map(is_valid, flist)
flist = [f for (f, v) in zip(flist, vlist) if v]
with open(OUTPUT, "w") as patch:
    print(f"process.mixData.input.fileNames = cms.untracked.vstring({repr(flist)})", file=patch)
