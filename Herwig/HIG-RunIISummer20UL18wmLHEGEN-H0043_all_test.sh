#!/bin/bash

source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc7_amd64_gcc700
[ -z "$1" ] && EVENTS=10 || EVENTS="$1"
[ -z "$2" ] && UPLOAD="" || UPLOAD="/eos/user/l/legao/NanoAODStore/V0/2018/MC/GluGluHToZZTo4L_M125_TuneCH3_13TeV_powheg2_JHUGenV7011_herwig7/RunIISummer20UL18wmLHEGEN-H0043/$2.root"
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

#enter CMSSW_10_6_28_patch1
#cmsDriver.py Configuration/GenProduction/python/HIG-RunIISummer20UL18wmLHEGEN-H0043-fragment.py --eventcontent RAWSIM,LHE --customise Configuration/DataProcessing/Utils.addMonitoring --datatier GEN,LHE --conditions 106X_upgrade2018_realistic_v4 --beamspot Realistic25ns13TeVEarly2018Collision --customise_commands process.source.numberEventsInLuminosityBlock="cms.untracked.uint32(100)" --step LHE,GEN --geometry DB:Extended --era Run2_2018 --python_filename HIG-RunIISummer20UL18wmLHEGEN-H0043_1_cfg.py --fileout file:HIG-RunIISummer20UL18wmLHEGEN-H0043.root --mc -n $EVENTS || exit $?
#
#enter CMSSW_10_6_17_patch1
#cmsDriver.py  --eventcontent RAWSIM --customise Configuration/DataProcessing/Utils.addMonitoring --datatier GEN-SIM --conditions 106X_upgrade2018_realistic_v11_L1v1 --beamspot Realistic25ns13TeVEarly2018Collision --step SIM --geometry DB:Extended --era Run2_2018 --python_filename HIG-RunIISummer20UL18SIM-00004_1_cfg.py --fileout file:HIG-RunIISummer20UL18SIM-00004.root --filein file:HIG-RunIISummer20UL18wmLHEGEN-H0043.root --runUnscheduled --mc -n $EVENTS || exit $?

enter CMSSW_10_6_17_patch1
cmsDriver.py  --eventcontent PREMIXRAW --customise Configuration/DataProcessing/Utils.addMonitoring --datatier GEN-SIM-DIGI --conditions 106X_upgrade2018_realistic_v11_L1v1 --step DIGI,DATAMIX,L1,DIGI2RAW --procModifiers premix_stage2 --geometry DB:Extended --datamix PreMix --era Run2_2018 --python_filename HIG-RunIISummer20UL18DIGIPremix-00004_1_cfg.py --fileout file:HIG-RunIISummer20UL18DIGIPremix-00004.root --filein file:HIG-RunIISummer20UL18SIM-00004.root --pileup_input "dbs:/Neutrino_E-10_gun/RunIISummer20ULPrePremix-UL18_106X_upgrade2018_realistic_v11_L1v1-v2/PREMIX" --runUnscheduled --mc -n $EVENTS --no_exec || exit $?
./patch_premix_inputs.py HIG-RunIISummer20UL18DIGIPremix-00004_1_cfg.py RunIISummer20ULPrePremix-UL18_106X_upgrade2018_realistic_v11_L1v1-v2.py.patch
cmsRun HIG-RunIISummer20UL18DIGIPremix-00004_1_cfg.py || exit $?

enter CMSSW_8_0_36_UL_patch1
cmsDriver.py  --eventcontent RAWSIM --customise Configuration/DataProcessing/Utils.addMonitoring --datatier GEN-SIM-RAW --conditions 102X_upgrade2018_realistic_v15 --customise_commands 'process.source.bypassVersionCheck = cms.untracked.bool(True)' --step HLT:2018v32 --geometry DB:Extended --era Run2_2018 --python_filename HIG-RunIISummer20UL18HLT-00004_1_cfg.py --fileout file:HIG-RunIISummer20UL18HLT-00004.root --filein file:HIG-RunIISummer20UL18DIGIPremix-00004.root --mc -n $EVENTS || exit $?

enter CMSSW_10_6_17_patch1
cmsDriver.py  --eventcontent AODSIM --customise Configuration/DataProcessing/Utils.addMonitoring --datatier AODSIM --conditions 106X_upgrade2018_realistic_v11_L1v1 --step RAW2DIGI,L1Reco,RECO,RECOSIM,EI --geometry DB:Extended --era Run2_2018 --python_filename HIG-RunIISummer20UL18RECO-00004_1_cfg.py --fileout file:HIG-RunIISummer20UL18RECO-00004.root --filein file:HIG-RunIISummer20UL18HLT-00004.root --runUnscheduled --mc -n $EVENTS || exit $?

enter CMSSW_10_6_25
cmsDriver.py  --eventcontent MINIAODSIM --customise Configuration/DataProcessing/Utils.addMonitoring --datatier MINIAODSIM --conditions 106X_upgrade2018_realistic_v16_L1v1 --step PAT --procModifiers run2_miniAOD_UL --geometry DB:Extended --era Run2_2018 --python_filename HIG-RunIISummer20UL18MiniAODv2-00166_1_cfg.py --fileout file:HIG-RunIISummer20UL18MiniAODv2-00166.root --filein file:HIG-RunIISummer20UL18RECO-00004.root --runUnscheduled --mc -n $EVENTS || exit $?

enter CMSSW_13_2_2
cmsDriver.py  --eventcontent NANOAODSIM --customise Configuration/DataProcessing/Utils.addMonitoring --datatier NANOAODSIM --conditions 106X_upgrade2018_realistic_v16_L1v1 --step NANO --era Run2_2018,run2_nanoAOD_106Xv2 --python_filename HIG-RunIISummer20UL18NanoAODv9-00162_1_cfg.py --fileout file:HIG-RunIISummer20UL18NanoAODv9-00162.root --filein file:HIG-RunIISummer20UL18MiniAODv2-00166.root --mc -n $EVENTS || exit $?

[ ! -z "${UPLOAD}" ] && xrdcp -f file:HIG-RunIISummer20UL18NanoAODv9-00162.root "root://eosuser.cern.ch/${UPLOAD}"
