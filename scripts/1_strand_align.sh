#!/bin/sh

#A script to fix strand, position etc

file=$1
#File name stem for the input files 
instem=$2

#File name stem for the output files 
outstem=$3


echo "Input stem is $instem"
echo "Output stem is $outstem"

#Cut the strand file into a series of Plink slices
rsids_file=$instem.rsids
chr_file=$instem.chr
pos_file=$instem.pos
flip_file=$instem.flip
delete_atcg=$instem.delete
alleles_file=$instem.1kgAlleles

###List of rsids to keep
cat $chr_file | cut -f 1 > $rsids_file


#Because Plink only allows you to update one attribute at a time, we need lots of temp
#Plink files
temp_prefix=TEMP_FILE_$(expr $RANDOM \* $RANDOM)_
temp1=$temp_prefix"1"
temp2=$temp_prefix"2"
temp3=$temp_prefix"3"
temp4=$temp_prefix"4"

## VERSION 4
#1) Only keep the valid rsIDS which are in 1000Genomes
echo "plink --noweb --allow-no-sex --bfile ${file} --extract "${rsids_file}" --make-bed --out $temp1"
#plink --noweb --allow-no-sex --bfile ${file} --extract ${rsids_file} --make-bed --out $temp1
plink2 --allow-no-sex --bfile ${file} --extract ${rsids_file} --make-bed --out $temp1  

#2)  Apply the chr
echo "plink  --noweb --allow-no-sex --bfile $temp1 --update-map $chr_file --update-chr --make-bed --out $temp2"
#plink  --noweb --allow-no-sex --bfile $temp1 --update-map $chr_file --update-chr --make-bed --out $temp2
plink2 --allow-no-sex --bfile $temp1 --update-map $chr_file --update-chr --make-bed --out $temp2 
#3) Apply the pos
echo "plink --noweb --allow-no-sex --bfile $temp2 --update-map $pos_file --make-bed --out $temp3"
#plink --noweb --allow-no-sex --bfile $temp2 --update-map $pos_file --make-bed --out $temp3
plink2 --allow-no-sex --bfile $temp2 --update-map $pos_file --make-bed --out $temp3 
#4) Flip to RevComp
echo "plink  --bfile  $temp3 --flip  $flip_file  --make-bed --out $temp4 --noweb"
#plink --bfile  $temp3 --flip  $flip_file  --make-bed --out $temp4 --noweb
plink2  --bfile  $temp3 --flip  $flip_file  --make-bed --out $temp4 --noweb --a2-allele $alleles_file 3 2 #reference allele is A2 in plink

#5) Exclude A-T & C-G that we cant flip as ambigous MAF
echo "plink  --bfile  $temp4 --exclude $delete_atcg   --make-bed --out $outstem --noweb "
#plink --bfile  $temp4 --exclude $delete_atcg   --make-bed --out $outstem --noweb 
plink2 --bfile  $temp4 --exclude $delete_atcg   --make-bed --out $outstem --a2-allele $alleles_file 3 2

