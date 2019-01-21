#!/bin/bash

bin=${1}
tag=${2}

for ord in {5..6..1}; do
    for xbin in {5..45..5}; do
	for ybin in {3..9..1}; do
	    for zbin in {5..45..5}; do
		# background on lxplus
		# source run_projEff_singleBin.sh ${bin} ${tag} ${ord} ${xbin} ${ybin} ${zbin} \
		#     &> logs/log_sub_projEff_singleBin_${bin}_${tag}_${ord}_${xbin}_${ybin}_${zbin}.log &

		# condor jobs
		cat > csub/sub_${bin}_${tag}_${ord}_${xbin}_${ybin}_${zbin}.jdl << EOF
Executable = run_projEff_singleBin.sh
Arguments  = ${bin} ${tag} ${ord} ${xbin} ${ybin} ${zbin}
Log        = logs/log_sub_projEff_singleBin_${bin}_${tag}_${ord}_${xbin}_${ybin}_${zbin}.log
Output     = logs/log_sub_projEff_singleBin_${bin}_${tag}_${ord}_${xbin}_${ybin}_${zbin}.out
Error      = logs/log_sub_projEff_singleBin_${bin}_${tag}_${ord}_${xbin}_${ybin}_${zbin}.out
Queue
EOF
		condor_submit csub/sub_${bin}_${tag}_${ord}_${xbin}_${ybin}_${zbin}.jdl
		
	    done
	done
    done
done
