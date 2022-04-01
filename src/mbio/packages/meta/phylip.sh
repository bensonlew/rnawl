#!/bin/bash

#$1 phylip.file
software=$2
bootstrap=$3
#software="/mnt/ilustre/users/sanger-dev/sg-users/houshuang/phylo_tree/phylo_software/phylip-3.697/exe"  #here
#outpath="/mnt/ilustre/users/sanger-dev/sg-users/houshuang/phylo_tree/mp_test/phylip"  #here
record_file="record.out"

echo "bootstrap:" >> $record_file
starttime=`date +'%Y-%m-%d %H:%M:%S'`
echo 'start:'$starttime >> $record_file
#nt, bootstrap
echo -e "$1\nR\n$bootstrap\nY\n5\n" | $software"/seqboot"
mv outfile infile
#D     Distance (F84, Kimura, Jukes-Cantor, LogDet)?  F84
#G     Gamma distributed rates across sites?  No
echo -e "M\nD\n$bootstrap\n5\n1\nY\n" | $software"/dnapars"
mv infile mp.seq
mv outfile dnamp
mv outtree intree
echo -e "Y\n" | $software"/consense"
mv intree dnamp_tree
mv outtree mp.bootstrap
mv outfile mp.consense
endtime=`date +'%Y-%m-%d %H:%M:%S'`
echo 'end:'$endtime >> $record_file	

#nt, branch length
echo "length:" >> $record_file
echo -e "$1\nY\n" | $software"/dnapars"
mv outfile mp.length_out
mv outtree mp.tree
endtime=`date +'%Y-%m-%d %H:%M:%S'`
echo 'end:'$endtime >> $record_file
echo " " >> $record_file
