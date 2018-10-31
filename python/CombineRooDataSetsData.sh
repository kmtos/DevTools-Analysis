for i in BSUB/MiniAOD_SingleMu1*.root; 
do
  echo ""
  echo "$i"
  sed -i "s|SUBSCRIPT|${i##*SingleMu1_}|g" rootMacro_CombineRooDataSets_Data.C
  root -l rootMacro_CombineRooDataSets_Data.C 
  sed -i "s|${i##*SingleMu1_}|SUBSCRIPT|g" rootMacro_CombineRooDataSets_Data.C
done

