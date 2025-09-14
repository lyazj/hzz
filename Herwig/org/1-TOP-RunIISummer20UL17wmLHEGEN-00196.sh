#!/bin/bash

# Binds for singularity containers
# Mount /afs, /eos, /cvmfs, /etc/grid-security for xrootd
export APPTAINER_BINDPATH='/afs,/cvmfs,/cvmfs/grid.cern.ch/etc/grid-security:/etc/grid-security,/eos,/etc/pki/ca-trust,/run/user,/var/run/user'

cat <<'EndOfGenScriptFile' > TOP-RunIISummer20UL17wmLHEGEN-00196_gen_script.sh
#!/bin/bash

echo "Running CMS GEN request script using cms-sw containers. Architecture: el9:x86_64"
python3 -m venv cms_gen_venv_TOP-RunIISummer20UL17wmLHEGEN-00196 && source ./cms_gen_venv_TOP-RunIISummer20UL17wmLHEGEN-00196/bin/activate

# Install the PdmV REST client
pip install git+https://github.com/cms-PdmV/mcm_scripts &> /dev/null

echo "Packages installed"
pip freeze
echo ""

# GEN Script begin
rm -f request_fragment_check.py
wget -q https://gitlab.cern.ch/cms-gen/genproductions_scripts/-/raw/master/bin/utils/request_fragment_check.py
chmod +x request_fragment_check.py

./request_fragment_check.py --bypass_status --prepid TOP-RunIISummer20UL17wmLHEGEN-00196

# End of CMS GEN script file: TOP-RunIISummer20UL17wmLHEGEN-00196_gen_script.sh
EndOfGenScriptFile
chmod +x TOP-RunIISummer20UL17wmLHEGEN-00196_gen_script.sh

# Run in singularity container
singularity run --home $PWD:$PWD /cvmfs/unpacked.cern.ch/registry.hub.docker.com/cmssw/el9:x86_64 $(echo $(pwd)/TOP-RunIISummer20UL17wmLHEGEN-00196_gen_script.sh)

GEN_ERR=$?
if [ $GEN_ERR -ne 0 ]; then
  echo "GEN Checking Script returned exit code $GEN_ERR which means there are $GEN_ERR errors"
  echo "Validation WILL NOT RUN"
  echo "Please correct errors in the request and run validation again"
  exit $GEN_ERR
fi
echo "Running VALIDATION. GEN Request Checking Script returned no errors"
# GEN Script end

# Download fragment from McM
curl -s -k https://cms-pdmv-prod.web.cern.ch/mcm/public/restapi/requests/get_fragment/TOP-RunIISummer20UL17wmLHEGEN-00196 --retry 3 --create-dirs -o Configuration/GenProduction/python/TOP-RunIISummer20UL17wmLHEGEN-00196-fragment.py
[ -s Configuration/GenProduction/python/TOP-RunIISummer20UL17wmLHEGEN-00196-fragment.py ] || exit $?;

# Check if fragment contais gridpack path ant that it is in cvmfs
if grep -q "gridpacks" Configuration/GenProduction/python/TOP-RunIISummer20UL17wmLHEGEN-00196-fragment.py; then
  if ! grep -q -e "/cvmfs/cms.cern.ch/phys_generator/gridpacks" -e "/cvmfs/cms-griddata.cern.ch/phys_generator/gridpacks_tarball" Configuration/GenProduction/python/TOP-RunIISummer20UL17wmLHEGEN-00196-fragment.py; then
    echo "Gridpack inside fragment is not in cvmfs."
    exit -1
  fi
fi

# Dump actual test code to a TOP-RunIISummer20UL17wmLHEGEN-00196_test.sh file that can be run in Singularity
cat <<'EndOfTestFile' > TOP-RunIISummer20UL17wmLHEGEN-00196_test.sh
#!/bin/bash

export SCRAM_ARCH=slc7_amd64_gcc700

source /cvmfs/cms.cern.ch/cmsset_default.sh
if [ -r CMSSW_10_6_24/src ] ; then
  echo release CMSSW_10_6_24 already exists
else
  scram p CMSSW CMSSW_10_6_24
fi
cd CMSSW_10_6_24/src
eval `scram runtime -sh`

mv ../../Configuration .
scram b
cd ../..

# Maximum validation runtime: 28800s
# Minimum validation runtime: 600s
# Output events to run for the validation job (from application's setting): 100
# Event efficiency: Computed using the request efficiency and its error.
# Event efficiency: `efficiency - (2 * efficiency_error)`: `1 - (2 * 0)` = 1
# Input events: `int(output_events / event_efficiency)`: `int(100 / 1)` = 100
# Time per event (s): Computed adding all the time_per_event values on every sequence
# Time per event (s): 0.861
# Initial target input events: 100
# Initial target output events: 100
# Validation runtime will not run for long enough than expected, extending the time
# Target input events changed to: `minimum_runtime / time_per_event * number_of_threads`: `600 / 0.861 * 1` = 697
# Target output events changed to: `target_input_events * event_efficiency`: `697 * 1` = 697
# Final target input events: 696
# Final target output events: 696
# This validation will be computed based on the target output events!
EVENTS=696


# cmsDriver command
cmsDriver.py Configuration/GenProduction/python/TOP-RunIISummer20UL17wmLHEGEN-00196-fragment.py --eventcontent RAWSIM,LHE --customise Configuration/DataProcessing/Utils.addMonitoring --datatier GEN,LHE --conditions 106X_mc2017_realistic_v6 --beamspot Realistic25ns13TeVEarly2017Collision --customise_commands process.source.numberEventsInLuminosityBlock="cms.untracked.uint32(100)" --step LHE,GEN --geometry DB:Extended --era Run2_2017 --python_filename TOP-RunIISummer20UL17wmLHEGEN-00196_1_cfg.py --fileout file:TOP-RunIISummer20UL17wmLHEGEN-00196.root --number 696 --number_out 696 --no_exec --mc || exit $? ;

# Run generated config
REPORT_NAME=TOP-RunIISummer20UL17wmLHEGEN-00196_report.xml
# Run the cmsRun
cmsRun -e -j $REPORT_NAME TOP-RunIISummer20UL17wmLHEGEN-00196_1_cfg.py || exit $? ;

# Parse values from TOP-RunIISummer20UL17wmLHEGEN-00196_report.xml report
processedEvents=$(grep -Po "(?<=<Metric Name=\"NumberEvents\" Value=\")(.*)(?=\"/>)" $REPORT_NAME | tail -n 1)
producedEvents=$(grep -Po "(?<=<TotalEvents>)(\d*)(?=</TotalEvents>)" $REPORT_NAME | tail -n 1)
threads=$(grep -Po "(?<=<Metric Name=\"NumberOfThreads\" Value=\")(.*)(?=\"/>)" $REPORT_NAME | tail -n 1)
peakValueRss=$(grep -Po "(?<=<Metric Name=\"PeakValueRss\" Value=\")(.*)(?=\"/>)" $REPORT_NAME | tail -n 1)
peakValueVsize=$(grep -Po "(?<=<Metric Name=\"PeakValueVsize\" Value=\")(.*)(?=\"/>)" $REPORT_NAME | tail -n 1)
totalSize=$(grep -Po "(?<=<Metric Name=\"Timing-tstoragefile-write-totalMegabytes\" Value=\")(.*)(?=\"/>)" $REPORT_NAME | tail -n 1)
totalSizeAlt=$(grep -Po "(?<=<Metric Name=\"Timing-file-write-totalMegabytes\" Value=\")(.*)(?=\"/>)" $REPORT_NAME | tail -n 1)
totalJobTime=$(grep -Po "(?<=<Metric Name=\"TotalJobTime\" Value=\")(.*)(?=\"/>)" $REPORT_NAME | tail -n 1)
totalJobCPU=$(grep -Po "(?<=<Metric Name=\"TotalJobCPU\" Value=\")(.*)(?=\"/>)" $REPORT_NAME | tail -n 1)
eventThroughput=$(grep -Po "(?<=<Metric Name=\"EventThroughput\" Value=\")(.*)(?=\"/>)" $REPORT_NAME | tail -n 1)
avgEventTime=$(grep -Po "(?<=<Metric Name=\"AvgEventTime\" Value=\")(.*)(?=\"/>)" $REPORT_NAME | tail -n 1)
if [ -z "$threads" ]; then
  echo "Could not find NumberOfThreads in report, defaulting to 1"
  threads=1
fi
if [ -z "$eventThroughput" ]; then
  eventThroughput=$(bc -l <<< "scale=4; 1 / ($avgEventTime / $threads)")
fi
if [ -z "$totalSize" ]; then
  totalSize=$totalSizeAlt
fi
if [ -z "$processedEvents" ]; then
  processedEvents=$EVENTS
fi
echo "Validation report of TOP-RunIISummer20UL17wmLHEGEN-00196 sequence 1/1"
echo "Processed events: $processedEvents"
echo "Produced events: $producedEvents"
echo "Threads: $threads"
echo "Peak value RSS: $peakValueRss MB"
echo "Peak value Vsize: $peakValueVsize MB"
echo "Total size: $totalSize MB"
echo "Total job time: $totalJobTime s"
echo "Total CPU time: $totalJobCPU s"
echo "Event throughput: $eventThroughput"
echo "CPU efficiency: "$(bc -l <<< "scale=2; ($totalJobCPU * 100) / ($threads * $totalJobTime)")" %"
echo "Size per event: "$(bc -l <<< "scale=4; ($totalSize * 1024 / $producedEvents)")" kB"
echo "Time per event: "$(bc -l <<< "scale=4; (1 / $eventThroughput)")" s"
echo "Filter efficiency percent: "$(bc -l <<< "scale=8; ($producedEvents * 100) / $processedEvents")" %"
echo "Filter efficiency fraction: "$(bc -l <<< "scale=10; ($producedEvents) / $processedEvents")

# End of TOP-RunIISummer20UL17wmLHEGEN-00196_test.sh file
EndOfTestFile

# Make file executable
chmod +x TOP-RunIISummer20UL17wmLHEGEN-00196_test.sh

if [ -e "/cvmfs/unpacked.cern.ch/registry.hub.docker.com/cmssw/el7:amd64" ]; then
  CONTAINER_NAME="el7:amd64"
elif [ -e "/cvmfs/unpacked.cern.ch/registry.hub.docker.com/cmssw/el7:x86_64" ]; then
  CONTAINER_NAME="el7:x86_64"
else
  echo "Could not find amd64 or x86_64 for el7"
  exit 1
fi
# Run in singularity container
export SINGULARITY_CACHEDIR="/tmp/$(whoami)/singularity"
singularity run --home $PWD:$PWD /cvmfs/unpacked.cern.ch/registry.hub.docker.com/cmssw/$CONTAINER_NAME $(echo $(pwd)/TOP-RunIISummer20UL17wmLHEGEN-00196_test.sh)
