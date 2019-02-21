#!/bin/bash

cd /afs/cern.ch/work/a/aboletti/private/Kstmumu-Run2/CMSSW_8_0_30/src
source  /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scram runtime -sh`

cd /afs/cern.ch/work/a/aboletti/private/Kstmumu-Run2/eff-projection/
pwd

bin=${1}
tag=${2}
ord=${3}
xbin=${4}
ybin=${5}
zbin=${6}

root -l -b -q 'projEff_spHarm_fromDataset.cc('${bin}','${tag}','${ord}','${xbin}','${ybin}','${zbin}',false)' \
    &>logs/log_sub_projEff_singleBin_${bin}_${tag}_${ord}_${xbin}_${ybin}_${zbin}.out
root -l -b -q 'plotEff_fromDataset.cc('${bin}','${tag}','${ord}','${xbin}','${ybin}','${zbin}',false)' \
    &>>logs/log_sub_projEff_singleBin_${bin}_${tag}_${ord}_${xbin}_${ybin}_${zbin}.out
IFS=$' ' command eval 'XYZ=($(tail -n1 logs/log_sub_projEff_singleBin_'${bin}'_'${tag}'_'${ord}'_'${xbin}'_'${ybin}'_'${zbin}'.out))'
echo "${xbin} ${ybin} ${zbin} ${XYZ[6]}" >> list_neg_eff_${bin}_${tag}_${ord}.list 
if [ "${XYZ[6]}" = "0" ]; then
    ./fit_recoMC_projEff ${bin} ${tag} ${ord} ${xbin} ${ybin} ${zbin} 0 \
	&>logs/log_fit_recoMC_projEff_${bin}_${tag}_${ord}_${xbin}_${ybin}_${zbin}.out
fi
echo "done"
