#!/bin/bash

#parse arguments
if [ $# -ne 5 ]
    then
    echo "Usage: ./generate.sh cfg_name script_name sample_tag queue name_addon "
    exit 0
fi

cfg_name=$1
script_name=$2
tag=$3
queue=$4
name_addon=$5
input="/afs/cern.ch/work/k/ktos/public/CMSSW_8_0_17/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/FILE_TESTS//InputFor${tag}.txt"

while IFS= read -r line
do
  linAr=($line)
  divisions=${linAr[0]}
  dir_name=${linAr[1]}
  file_path=${linAr[2]}
  full_dir_name="${linAr[1]}${name_addon}"
  echo ""
  echo ""
  mkdir -p BSUB/$full_dir_name
  cd BSUB/$full_dir_name
  xsec=`grep "$dir_name" /afs/cern.ch/work/k/ktos/public/CMSSW_8_0_17/src/AnalyzerGeneratorRecoVariousFunctions/Analyzer/FILE_TESTS/CrossSections.py`
  xsec=${xsec##* }
  echo "xsec=$xsec"  
  summedWeights=`grep "$dir_name" /afs/cern.ch/user/k/ktos/NewSkimDir/CMSSW_8_0_30/src/GGHAA2Mu2TauAnalysis/SkimMuMuTauTau/test/SkimSequence/SummedWeightsFiles/SummedWeightsValues.out`
  summedWeights=${summedWeights##* }
  echo "summedWeights=$summedWeights"
  echo "Dir= $full_dir_name    count=$divisions"
  echo "Sub_dir_name = $dir_name"
  echo "File_path= $file_path"
  temp=${file_path##*ktos/} 
  pkl_path=${temp%%/*}
 
  divisions=$((divisions + 5))
  echo "python file= ${cfg_name}_${full_dir_name}.py"
  echo "Script name= ${script_name}_${full_dir_name}.sh"
  sed -e "s|SAMPLE_NAME|$pkl_path|g" -e "s|XSECVALUE|$xsec|g" -e "s|SUMMEDWEIGHTSVALUE|${summedWeights}|g" -e "s|FILE_PATH|${file_path}|g" -e "s|DIRNAME|${full_dir_name}|g" -e "s|ENDFILENUMVALUE|${divisions}|g"  ../../${cfg_name}.py > ${cfg_name}_${full_dir_name}.py
  sed -e "s|ANALYZER|${cfg_name}_${full_dir_name}|g" -e "s|DIRNAME|${full_dir_name}|g"  ../../${script_name}.sh > ${script_name}_${full_dir_name}.sh
  chmod u+x ${script_name}_${full_dir_name}.sh
  bsub -q $queue -J ${cfg_name}_${full_dir_name} < ${script_name}_${full_dir_name}.sh
  cd ../../
done <"$input"
exit 0
