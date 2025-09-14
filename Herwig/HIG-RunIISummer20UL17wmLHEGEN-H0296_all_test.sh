#!/bin/bash

source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc7_amd64_gcc700
[ -z "$1" ] && EVENTS=10 || EVENTS="$1"
[ -z "$2" ] && UPLOAD="" || UPLOAD="/eos/user/l/legao/NanoAODStore/V0/2017/MC/GluGluHToZZTo4L_M125_TuneCH3_13TeV_powheg2_JHUGenV7011_herwig7/RunIISummer20UL17wmLHEGEN-H0296/$2.root"
echo "Events: ${EVENTS}"
echo "Upload: ${UPLOAD}"

voms-proxy-info || exit

enter () {
  [ -r $1/src ] || scram p CMSSW $1
  cd $1/src
  cp -r ../../Configuration .
  eval `scram runtime -sh`
  scram b
  cd ../..
}

enter CMSSW_10_6_28_patch1
cmsDriver.py Configuration/GenProduction/python/HIG-RunIISummer20UL17wmLHEGEN-H0296-fragment.py --eventcontent RAWSIM,LHE --customise Configuration/DataProcessing/Utils.addMonitoring --datatier GEN,LHE --conditions 106X_mc2017_realistic_v6 --beamspot Realistic25ns13TeVEarly2017Collision --customise_commands process.source.numberEventsInLuminosityBlock="cms.untracked.uint32(100)" --step LHE,GEN --geometry DB:Extended --era Run2_2017 --python_filename HIG-RunIISummer20UL17wmLHEGEN-H0296_1_cfg.py --fileout file:HIG-RunIISummer20UL17wmLHEGEN-H0296.root --mc -n $EVENTS || exit $?

enter CMSSW_10_6_17_patch1
cmsDriver.py  --eventcontent RAWSIM --customise Configuration/DataProcessing/Utils.addMonitoring --datatier GEN-SIM --conditions 106X_mc2017_realistic_v6 --beamspot Realistic25ns13TeVEarly2017Collision --step SIM --geometry DB:Extended --era Run2_2017 --python_filename HIG-RunIISummer20UL17SIM-00599_1_cfg.py --fileout file:HIG-RunIISummer20UL17SIM-00599.root --filein file:HIG-RunIISummer20UL17wmLHEGEN-H0296.root --runUnscheduled --mc -n $EVENTS || exit $?

enter CMSSW_10_6_17_patch1
cmsDriver.py  --eventcontent PREMIXRAW --customise Configuration/DataProcessing/Utils.addMonitoring --datatier GEN-SIM-DIGI --conditions 106X_mc2017_realistic_v6 --step DIGI,DATAMIX,L1,DIGI2RAW --procModifiers premix_stage2 --geometry DB:Extended --datamix PreMix --era Run2_2017 --python_filename HIG-RunIISummer20UL17DIGIPremix-00599_1_cfg.py --fileout file:HIG-RunIISummer20UL17DIGIPremix-00599.root --filein file:HIG-RunIISummer20UL17SIM-00599.root --pileup_input "dbs:/Neutrino_E-10_gun/RunIISummer20ULPrePremix-UL17_106X_mc2017_realistic_v6-v3/PREMIX" --runUnscheduled --mc -n $EVENTS --no_exec || exit $?
./patch_premix_inputs.py HIG-RunIISummer20UL17DIGIPremix-00599_1_cfg.py RunIISummer20ULPrePremix-UL17_106X_mc2017_realistic_v6-v3.py.patch
cmsRun HIG-RunIISummer20UL17DIGIPremix-00599_1_cfg.py || exit $?

enter CMSSW_9_4_14_UL_patch1
cmsDriver.py  --eventcontent RAWSIM --customise Configuration/DataProcessing/Utils.addMonitoring --datatier GEN-SIM-RAW --conditions 94X_mc2017_realistic_v15 --customise_commands 'process.source.bypassVersionCheck = cms.untracked.bool(True)' --step HLT:2e34v40 --geometry DB:Extended --era Run2_2017 --python_filename HIG-RunIISummer20UL17HLT-00599_1_cfg.py --fileout file:HIG-RunIISummer20UL17HLT-00599.root --filein file:HIG-RunIISummer20UL17DIGIPremix-00599.root --mc -n $EVENTS || exit $?

enter CMSSW_10_6_17_patch1
cmsDriver.py  --eventcontent AODSIM --customise Configuration/DataProcessing/Utils.addMonitoring --datatier AODSIM --conditions 106X_mc2017_realistic_v6 --step RAW2DIGI,L1Reco,RECO,RECOSIM --geometry DB:Extended --era Run2_2017 --python_filename HIG-RunIISummer20UL17RECO-00599_1_cfg.py --fileout file:HIG-RunIISummer20UL17RECO-00599.root --filein file:HIG-RunIISummer20UL17HLT-00599.root --runUnscheduled --mc -n $EVENTS || exit $?

enter CMSSW_10_6_20
cmsDriver.py  --eventcontent MINIAODSIM --customise Configuration/DataProcessing/Utils.addMonitoring --datatier MINIAODSIM --conditions 106X_mc2017_realistic_v9 --step PAT --procModifiers run2_miniAOD_UL --geometry DB:Extended --era Run2_2017 --python_filename HIG-RunIISummer20UL17MiniAODv2-00675_1_cfg.py --fileout file:HIG-RunIISummer20UL17MiniAODv2-00675.root --filein file:HIG-RunIISummer20UL17RECO-00599.root --runUnscheduled --mc -n $EVENTS || exit $?

enter CMSSW_13_2_2
cmsDriver.py  --eventcontent NANOAODSIM --customise Configuration/DataProcessing/Utils.addMonitoring --datatier NANOAODSIM --conditions 106X_mc2017_realistic_v9 --step NANO --era Run2_2017,run2_nanoAOD_106Xv2 --python_filename HIG-RunIISummer20UL17NanoAODv12-00641_1_cfg.py --fileout file:HIG-RunIISummer20UL17NanoAODv12-00641.root --filein file:HIG-RunIISummer20UL17MiniAODv2-00675.root --mc -n $EVENTS || exit $?

[ ! -z "${UPLOAD}" ] && xrdcp -f HIG-RunIISummer20UL17NanoAODv12-00641.root "root://eosuser.cern.ch/${UPLOAD}"
