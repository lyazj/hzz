import FWCore.ParameterSet.Config as cms

externalLHEProducer = cms.EDProducer("ExternalLHEProducer",
    args = cms.vstring('/cvmfs/cms.cern.ch/phys_generator/gridpacks/2017/13TeV/powheg/V2/gg_H_quark-mass-effects_ZZ_NNPDF31_13TeV/gg_H_quark-mass-effects_NNPDF31_13TeV_M125/v3/gg_H_quark-mass-effects_NNPDF31_13TeV_M125.tgz'),
    nEvents = cms.untracked.uint32(5000),
    numberOfParameters = cms.uint32(1),
    outputFile = cms.string('cmsgrid_final.lhe'),
    scriptName = cms.FileInPath('GeneratorInterface/LHEInterface/data/run_generic_tarball_cvmfs.sh')
)

# Link to cards:
# https://raw.githubusercontent.com/cms-sw/genproductions/fd7d34a91c3160348fd0446ded445fa28f555e09/bin/Powheg/production/2017/13TeV/Higgs/gg_H_ZZ_quark-mass-effects_NNPDF31_13TeV/makecards.py
#    gg_H_ZZ_quark-mass-effects_NNPDF31_13TeV_M125.input
# https://raw.githubusercontent.com/cms-sw/genproductions/fd7d34a91c3160348fd0446ded445fa28f555e09/bin/JHUGen/cards/decay/ZZ4l_withtaus.input

from Configuration.Generator.Herwig7Settings.Herwig7LHECommonSettings_cfi import *
from Configuration.Generator.Herwig7Settings.Herwig7StableParticlesForDetector_cfi import *
from Configuration.Generator.Herwig7Settings.Herwig7CH3TuneSettings_cfi import *
from Configuration.Generator.Herwig7Settings.Herwig7LHEPowhegSettings_cfi import *
from Configuration.Generator.Herwig7Settings.Herwig7PSWeightsSettings_cfi import *
from Configuration.Generator.Herwig7Settings.Herwig7_7p1SettingsFor7p2_cfi import *

generator = cms.EDFilter("Herwig7GeneratorFilter",
    herwig7LHECommonSettingsBlock,
    herwig7LHEPowhegSettingsBlock,
    herwig7p1SettingsFor7p2Block,
    herwig7StableParticlesForDetectorBlock,
    herwig7PSWeightsSettingsBlock,
    herwig7CH3SettingsBlock,
    configFiles = cms.vstring(),
    crossSection = cms.untracked.double(-1),
    dataLocation = cms.string('${HERWIGPATH:-6}'),
    eventHandlers = cms.string('/Herwig/EventHandlers'),
    filterEfficiency = cms.untracked.double(1.0),
    generatorModule = cms.string('/Herwig/Generators/EventGenerator'),
    hw_user_settings = cms.vstring(
        'cd /Herwig/EventHandlers',
        'set EventHandler:LuminosityFunction:Energy 13000*GeV',
        'cd /',
        'set /Herwig/Particles/h0:NominalMass 125.0'
    ),
    parameterSets = cms.vstring(
        'hw_lhe_common_settings',
        'hw_lhe_powheg_settings',
        'hw_7p1SettingsFor7p2',
        'herwig7CH3PDF',
        'herwig7CH3AlphaS',
        'herwig7CH3MPISettings',
        'herwig7StableParticlesForDetector',
        'hw_PSWeights_settings',
        'hw_user_settings'
    ),
    repository = cms.string('${HERWIGPATH}/HerwigDefaults.rpo'),
    run = cms.string('InterfaceMatchboxTest'),
    runModeList = cms.untracked.string('read,run'),
)
