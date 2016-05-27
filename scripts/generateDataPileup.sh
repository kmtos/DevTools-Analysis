#!/bin/bash
version=`python -c "from DevTools.Utilities.utilities import getCMSSWVersion; print getCMSSWVersion()"`
if [ "$version" == "76X" ]; then
    lumimask="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_v2.txt"
    pileupjson="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/PileUp/pileup_latest.txt"
else
    lumimask="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-273450_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt"
    pileupjson="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt"
fi
mkdir pileup
for xsec in 71000; do
    up=$(echo "$xsec*1.05" | bc)
    down=$(echo "$xsec*0.95" | bc)
    echo $xsec
    pileupCalc.py -i $lumimask --inputLumiJSON $pileupjson --calcMode true  --minBiasXsec $xsec --maxPileupBin 80 --numPileupBins 80 pileup/PileUpData.root
    echo $up
    pileupCalc.py -i $lumimask --inputLumiJSON $pileupjson --calcMode true  --minBiasXsec $up --maxPileupBin 80 --numPileupBins 80 pileup/PileUpData_up.root
    echo $down
    pileupCalc.py -i $lumimask --inputLumiJSON $pileupjson --calcMode true  --minBiasXsec $down --maxPileupBin 80 --numPileupBins 80 pileup/PileUpData_down.root
done
for xsec in 65000 66000 67000 68000 69000 70000 71000 72000 73000 74000 75000 76000 77000 78000 79000 80000; do
    echo $xsec
    pileupCalc.py -i $lumimask --inputLumiJSON $pileupjson --calcMode true  --minBiasXsec $xsec --maxPileupBin 80 --numPileupBins 80 pileup/PileUpData_$xsec.root
done
