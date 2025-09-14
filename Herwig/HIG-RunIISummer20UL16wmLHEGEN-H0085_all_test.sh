#!/bin/bash

source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc7_amd64_gcc700
[ -z "$1" ] && EVENTS=10 || EVENTS="$1"
[ -z "$2" ] && UPLOAD="" || UPLOAD="/eos/user/l/legao/NanoAODStore/V0/2016/MC/GluGluHToZZTo4L_M125_TuneCH3_13TeV_powheg2_JHUGenV7011_herwig7/HIG-RunIISummer20UL16wmLHEGEN-H0085/$2.root"
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
cmsDriver.py Configuration/GenProduction/python/HIG-RunIISummer20UL16wmLHEGEN-H0085-fragment.py --eventcontent RAWSIM,LHE --customise Configuration/DataProcessing/Utils.addMonitoring --datatier GEN,LHE --conditions 106X_mcRun2_asymptotic_v13 --beamspot Realistic25ns13TeV2016Collision --customise_commands process.source.numberEventsInLuminosityBlock="cms.untracked.uint32(100)" --step LHE,GEN --geometry DB:Extended --era Run2_2016 --python_filename HIG-RunIISummer20UL16wmLHEGEN-H0085_1_cfg.py --fileout file:HIG-RunIISummer20UL16wmLHEGEN-H0085.root --mc -n $EVENTS || exit $?

enter CMSSW_10_6_17_patch1
cmsDriver.py  --eventcontent RAWSIM --customise Configuration/DataProcessing/Utils.addMonitoring --datatier GEN-SIM --conditions 106X_mcRun2_asymptotic_v13 --beamspot Realistic25ns13TeV2016Collision --step SIM --geometry DB:Extended --era Run2_2016 --python_filename HIG-RunIISummer20UL16SIM-00498_1_cfg.py --fileout file:HIG-RunIISummer20UL16SIM-00498.root --filein file:HIG-RunIISummer20UL16wmLHEGEN-H0085.root --runUnscheduled --mc -n $EVENTS || exit $?

enter CMSSW_10_6_17_patch1
cmsDriver.py  --eventcontent PREMIXRAW --customise Configuration/DataProcessing/Utils.addMonitoring --datatier GEN-SIM-DIGI --conditions 106X_mcRun2_asymptotic_v13 --step DIGI,DATAMIX,L1,DIGI2RAW --procModifiers premix_stage2 --geometry DB:Extended --datamix PreMix --era Run2_2016 --python_filename HIG-RunIISummer20UL16DIGIPremix-00496_1_cfg.py --fileout file:HIG-RunIISummer20UL16DIGIPremix-00496.root --filein file:HIG-RunIISummer20UL16SIM-00498.root --pileup_input "dbs:/Neutrino_E-10_gun/RunIISummer20ULPrePremix-UL16_106X_mcRun2_asymptotic_v13-v1/PREMIX" --runUnscheduled --mc -n $EVENTS --no_exec || exit $?
cmsRun HIG-RunIISummer20UL16DIGIPremix-00496_1_cfg.py || exit $?

enter CMSSW_8_0_36_UL_patch1
cmsDriver.py  --eventcontent RAWSIM --outputCommand "keep *_mix_*_*,keep *_genPUProtons_*_*" --customise Configuration/DataProcessing/Utils.addMonitoring --datatier GEN-SIM-RAW --inputCommands "keep *","drop *_*_BMTF_*","drop *PixelFEDChannel*_*_*_*" --conditions 80X_mcRun2_asymptotic_2016_TrancheIV_v6 --customise_commands 'process.source.bypassVersionCheck = cms.untracked.bool(True)' --step HLT:25ns15e33_v4 --geometry DB:Extended --era Run2_2016 --python_filename HIG-RunIISummer20UL16HLT-00498_1_cfg.py --fileout file:HIG-RunIISummer20UL16HLT-00498.root --filein file:HIG-RunIISummer20UL16DIGIPremix-00496.root --mc -n $EVENTS || exit $?

enter CMSSW_10_6_17_patch1
cmsDriver.py  --eventcontent AODSIM --customise Configuration/DataProcessing/Utils.addMonitoring --datatier AODSIM --conditions 106X_mcRun2_asymptotic_v13 --step RAW2DIGI,L1Reco,RECO,RECOSIM --geometry DB:Extended --era Run2_2016 --python_filename HIG-RunIISummer20UL16RECO-00498_1_cfg.py --fileout file:HIG-RunIISummer20UL16RECO-00498.root --filein file:HIG-RunIISummer20UL16HLT-00498.root --runUnscheduled --mc -n $EVENTS || exit $?

enter CMSSW_10_6_25
cmsDriver.py  --eventcontent MINIAODSIM --customise Configuration/DataProcessing/Utils.addMonitoring --datatier MINIAODSIM --conditions 106X_mcRun2_asymptotic_v17 --step PAT --procModifiers run2_miniAOD_UL --geometry DB:Extended --era Run2_2016 --python_filename HIG-RunIISummer20UL16MiniAODv2-00490_1_cfg.py --fileout file:HIG-RunIISummer20UL16MiniAODv2-00490.root --filein file:HIG-RunIISummer20UL16RECO-00498.root --runUnscheduled --mc -n $EVENTS || exit $?

enter CMSSW_13_2_2
cmsDriver.py  --eventcontent NANOAODSIM --customise Configuration/DataProcessing/Utils.addMonitoring --datatier NANOAODSIM --conditions 106X_mcRun2_asymptotic_v17 --step NANO --era Run2_2016,run2_nanoAOD_106Xv2 --python_filename HIG-RunIISummer20UL16NanoAODv12-00528_1_cfg.py --fileout file:HIG-RunIISummer20UL16NanoAODv12-00528.root --filein file:HIG-RunIISummer20UL16MiniAODv2-00490.root --mc -n $EVENTS || exit $?

[ ! -z "${UPLOAD}" ] && xrdcp -f HIG-RunIISummer20UL16NanoAODv12-00528.root "root://eosuser.cern.ch/${UPLOAD}"
