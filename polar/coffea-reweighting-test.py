exec(open('coffea-reweighting.py').read())  # It's hard to import it...

rootfiles = open('example-files.txt').read().strip().split()
print(rootfiles)

for rootfile in rootfiles:
    events = NanoEventsFactory.from_root(f'root://cms-xrd-global.cern.ch/{rootfile}', treepath='Events').events()
    output = ReweightProcessor().process(events)
    print(output)
    break
