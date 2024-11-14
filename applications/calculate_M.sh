#!/bin/bash

pathname=MSA_Validation_Set/
filename=effective_validation_set.csv
echo "MSA,Meff" > $filename

for msa in $(ls ${pathname}*.faclean)
do
    echo "$(basename $msa .faclean),$(../../../../Programs/sequeff/build/sequeff -u -f $msa)" >> $filename
done
