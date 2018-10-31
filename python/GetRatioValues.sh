#!/bin/bash

for i in ./BSUB/*; 
do 
  var=`grep -r "Ratio of tau selection to no" $i/L*`
  echo ${i##*/}${var##*=} >> BSUB/RatioValues.out
done
