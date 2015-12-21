## run from base dir
#directory structure
#mkdir data temp strandAligned vcf

infile=$1
#filter samples missing overall >2% or SNPs missing >3%
plink2 --bfile data/${infile} --make-bed --out temp/data --allow-no-sex #--mind 0.02 --geno 0.03

#Split into chromosome files
parallel "plink2 --bfile temp/data --make-bed --out temp/data{}_tmp --allow-no-sex --chr {}" ::: $(seq 1 22)

#filter samples missing >5% per chromosome
parallel "plink2 --bfile temp/data{}_tmp --make-bed --out temp/data{} --allow-no-sex --chr {} #--mind 0.05" ::: $(seq 1 22)

#make sure all rs names are lower case
parallel 'perl -pi -e "s/RS/rs/g" temp/data{}.bim' ::: $(seq 1 22)

#strand align each chromosome
parallel 'instem=temp ; Rscript ~/Murray/Impute_pipeline/Impute.git/scripts/new_v4_ComparePOS1000GWithBimFiles.r {} temp/data{}.bim temp/${instem} ~/Murray \
        && ~/Murray/Impute_pipeline/Impute.git/scripts/1_strand_align.sh temp/data{} temp/${instem}_chr{} strandAligned/data{}_aligned {}' ::: $(seq 1 22)

cat temp/*1kgAlleles > vcf/1kgAlleles.txt

#make vcf
parallel 'plink2 --bfile strandAligned/data{}_aligned --keep-allele-order --recode vcf --out vcf/data{}_strandaligned' ::: $(seq 1 22)
## work out what samples are missing and then remove all missing samples from every chromosome
parallel 'zgrep ^#CHROM vcf/data{}_strandaligned.vcf.gz | cut -f 10- | sed "s/\t/\n/g" > vcf/sample{}.txt ' ::: $(seq 1 22)
parallel 'bcftools view -s {2} --force-samples vcf/data{1}_strandaligned.vcf.gz -o vcf/data{1}_strandaligned.samples_removed.vcf.gz -O z ' ::: $(seq 1 22) ::: ^$(cat vcf/sample* | sort | uniq -c | awk '{if($1 != 22){print $2}}' | tr '\n' ',')
#merge all chromosomes into a single vcf
parallel 'bgzip vcf/data{}_strandaligned.vcf && tabix -p vcf vcf/data{}_strandaligned.vcf.gz' ::: $(seq 1 22)
bcftools concat -O z -o vcf/all_chr_strandaligned.vcf.gz $(for i in $(seq 1 22); do echo vcf/data${i}_strandaligned.vcf.gz ; done | tr '\n' ' ') 


#QC steps:
#remove positions that can't have ref allele set properly
grep "Warning: Impossible A. allele assignment for variant" strandAligned/*log | cut -d' ' -f8 |  sed "s/\.//g" | sort -u > vcf/remove.txt
zgrep -v "#" < vcf/all_chr_strandaligned.vcf.gz | cut -f3 | cat - vcf/remove.txt  | sort | uniq -d |tr "\n" "," > vcf/rs_remove.txt
if [ $(wc -l vcf/rs_remove.txt | cut -c1) -eq 0 ]
then
	echo NO SNPS TO REMOVE!
	plink2 --vcf vcf/all_chr_strandaligned.vcf.gz --recode vcf --out vcf/all_chr_strandaligned.filt_refallele_maf0.001 --a2-allele vcf/1kgAlleles.txt 3 2 --maf 0.001
else
	echo Removing SNPs
	plink2 --vcf vcf/all_chr_strandaligned.vcf.gz --recode vcf --out vcf/all_chr_strandaligned.filt_refallele_maf0.001 --a2-allele vcf/1kgAlleles.txt 3 2 --exclude-snps $(cat vcf/rs_remove.txt) --maf 0.001
fi
#filter on missingness >2% total
#plink2 --vcf vcf/all_chr_strandaligned.filt_refallele_maf0.001.vcf --a2-allele vcf/1kgAlleles.txt 3 2 --mind 0.02 --recode vcf --out vcf/all_chr_strandaligned.filt_refallele_maf0.001_mind0-02

#filter SNPs >3% geno missing
#plink2 --vcf vcf/all_chr_strandaligned.filt_refallele_maf0.001_mind0-02.vcf --a2-allele vcf/1kgAlleles.txt 3 2 --mind 0.02 --recode vcf --out vcf/all_chr_strandaligned.filt_refallele_maf0.001_mind0-02_geno0-03

python ~/Murray/Impute_pipeline/checkVCF.py -r ~/Murray/Bioinformatics/Reference_Files/FASTA/hs37d5/hs37d5.fa -o vcf_checking vcf/all_chr_strandaligned.filt_refallele_maf0.001.vcf
