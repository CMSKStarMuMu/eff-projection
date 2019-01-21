#!/bin/bash

cd /afs/cern.ch/work/a/aboletti/private/Kstmumu-Run2/CMSSW_8_0_30/src
source  /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scram runtime -sh`

cd /afs/cern.ch/work/a/aboletti/private/Kstmumu-Run2/eff-projection/
pwd

root -l -b -q 'projEff_spHarm_fromDataset.cc('${1}','${2}','${3}','${4}','${5}','${6}',false)'

root -l -b -q 'plotEff_fromDataset.cc('${1}','${2}','${3}','${4}','${5}','${6}',false)'
