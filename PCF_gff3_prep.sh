#!/bin/bash/

ml BEDOPS/2.4.39-foss-2019b
ml Python/3.8.2-GCCcore-8.3.0

wrk_dir=$PWD

for i in yalo yfil acam agam ambo ltul
  do
    for hap in hap1 hap2
      do
        cd "./$i/data/$hap/"
        for j in *.gff3
          do
            echo "Converting $j gff3..."
            gff2bed < $j > ${j%.*}'.bed'
            echo "Renaming $j bed IDs..."
            python $wrk_dir'/'bed_renamer.py ${j%.*}'.bed' > ${j%.*}'_slim_ids.bed'
            echo "Extracting $j gene lines..."
            grep -w "gene" ${j%.*}'_slim_ids.bed' > ${j%.*}'_slim_ids_gene.bed' 
          done
        cd $wrk_dir
      done
  done
