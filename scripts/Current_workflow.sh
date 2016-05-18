## run from base dir
#directory structure
#mkdir data temp strandAligned vcf

PREFIX='/Volumes/BiochemXsan/staff_groups/merrimanlab/Documents/'


infile=$1 # input bfile
type=$2 #'chip' or 'seq'
phase=$3 #'phase' or 'nophase'
marker=$4 # 'pos' or 'rs'
hwe=$5 #blank or 'no_hwe'
if [ $type == 'chip' ]
then
	echo "chip data - filtering missingness"
	#filter samples missing overall >2% or SNPs missing >3%
	plink2 --bfile data/${infile} --make-bed --out temp/data --allow-no-sex --mind 0.02 --geno 0.03 --missing

	#Split into chromosome files
	parallel "plink2 --bfile temp/data --make-bed --out temp/data{}_tmp --allow-no-sex --chr {}" ::: $(seq 1 22)

	#filter samples missing >5% per chromosome
	parallel "plink2 --bfile temp/data{}_tmp --make-bed --out temp/data{} --allow-no-sex --chr {} --mind 0.05 --missing" ::: $(seq 1 22)
fi

if [ $type == 'seq' ]
then
	echo "seq data - not filtering on missingness"
	#Don't filter on missingness
	plink2 --bfile data/${infile} --make-bed --out temp/data --allow-no-sex --missing

	#Split into chromosome files
	parallel "plink2 --bfile temp/data --make-bed --out temp/data{}_tmp --allow-no-sex --chr {}" ::: $(seq 1 22)

	#filter samples missing >5% per chromosome
	parallel "plink2 --bfile temp/data{}_tmp --make-bed --out temp/data{} --allow-no-sex --chr {} --missing" ::: $(seq 1 22)
fi


#make sure all rs names are lower case
parallel 'perl -pi -e "s/RS/rs/g" temp/data{}.bim' ::: $(seq 1 22)


if [ $marker == 'rs' ]
then
#strand align each chromosome
parallel 'instem=temp ; Rscript ~/Murray/Impute_pipeline/Impute.git/scripts/new_v4_CompareRS100GWithBimFiles.r {} temp/data{}.bim temp/${instem} ~/Murray \
        && ~/Murray/Impute_pipeline/Impute.git/scripts/1_strand_align.sh temp/data{} temp/${instem}_chr{} strandAligned/data{}_aligned {}' ::: $(seq 1 22)
fi

if [ $marker == 'pos' ]
then
#strand align each chromosome
parallel 'instem=temp ; Rscript ~/Murray/Impute_pipeline/Impute.git/scripts/new_v4_ComparePOS1000GWithBimFiles.r {} temp/data{}.bim temp/${instem} ~/Murray \
        && ~/Murray/Impute_pipeline/Impute.git/scripts/1_strand_align.sh temp/data{} temp/${instem}_chr{} strandAligned/data{}_aligned {}' ::: $(seq 1 22)
fi

cat temp/*1kgAlleles > vcf/1kgAlleles.txt

#make vcf
parallel 'plink2 --bfile strandAligned/data{}_aligned --keep-allele-order --recode vcf --out vcf/data{}_strandaligned' ::: $(seq 1 22)


parallel 'bgzip vcf/data{}_strandaligned.vcf && tabix -p vcf vcf/data{}_strandaligned.vcf.gz' ::: $(seq 1 22)
## work out what samples are missing and then remove all missing samples from every chromosome
parallel 'zgrep ^#CHROM vcf/data{}_strandaligned.vcf.gz | cut -f 10- | sed "s/\t/\n/g" > vcf/sample{}.txt ' ::: $(seq 1 22)
parallel 'bcftools view -s {2} --force-samples vcf/data{1}_strandaligned.vcf.gz -o vcf/data{1}_strandaligned.samples_removed.vcf.gz -O z ' ::: $(seq 1 22) ::: ^$(cat vcf/sample* | sort | uniq -c | awk '{if($1 != 22){print $2}}' | tr '\n' ',')




#merge all chromosomes into a single vcf
bcftools concat -O z -o vcf/all_chr_strandaligned.vcf.gz $(for i in $(seq 1 22); do echo vcf/data${i}_strandaligned.samples_removed.vcf.gz ; done | tr '\n' ' ') 





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

#filter hwe
if [ $hwe == 'no_hwe' ]
then
	plink2 --vcf vcf/all_chr_strandaligned.filt_refallele_maf0.001.vcf --recode vcf --out vcf/all_chr_strandaligned.filt_refallele_maf0.01_hwe_checked --a2-allele vcf/1kgAlleles.txt 3 2 --maf 0.01 --missing 
else
	plink2 --vcf vcf/all_chr_strandaligned.filt_refallele_maf0.001.vcf --out vcf/hardy --a2-allele vcf/1kgAlleles.txt 3 2 --maf 0.01 --hardy
	Rscript ~/Murray/Impute_pipeline/hwe_filter.R vcf/hardy.hwe vcf/ bonferroni
	if [ $(wc -l vcf/hardy_exclude.txt | cut -c1) -eq 0 ]
	then
		echo NO HWE SNPS TO REMOVE!
		plink2 --vcf vcf/all_chr_strandaligned.filt_refallele_maf0.001.vcf --recode vcf --out vcf/all_chr_strandaligned.filt_refallele_maf0.01_hwe_checked --a2-allele vcf/1kgAlleles.txt 3 2 --maf 0.01 --missing

	else
		echo Removing HWE SNPs
		plink2 --vcf vcf/all_chr_strandaligned.filt_refallele_maf0.001.vcf --recode vcf --out vcf/all_chr_strandaligned.filt_refallele_maf0.01_hwe_checked --a2-allele vcf/1kgAlleles.txt 3 2 --exclude-snps $(cat vcf/hardy_exclude.txt | sort | uniq) --maf 0.01 --missing
	fi
fi

if [ $phase == 'phase' ]
then
	mkdir -p vcf/Phased
	## can phase with shapeit2 here
	for i in $(seq 1 22); do plink2 --vcf vcf/all_chr_strandaligned.filt_refallele_maf0.01_hwe_checked.vcf --recode vcf --out vcf/Phased/chr${i}_strandaligned.filt_refallele_maf0.01_hwe_checked --a2-allele vcf/1kgAlleles.txt 3 2 --chr ${i} ; done

	#parallel -j 2 '~/Murray/src/shapeit.v2.r837/bin/shapeit -M ~/Murray/Bioinformatics/Reference_Files/Impute_ref_files/1000GP_Phase3/genetic_map_chr{}_combined_b37.txt --input-vcf vcf/Phased/chr{}_strandaligned.filt_refallele_maf0.01_hwe_checked.vcf --output-max vcf/Phased/chr{}_strandaligned.filt_refallele_maf0.01_hwe_checked -T 6 ' ::: $(seq 1 22)

	parallel -j 1 ' ~/Murray/src/shapeit.v2.r837/bin/shapeit -check -M ~/Murray/Bioinformatics/Reference_Files/Impute_ref_files/1000GP_Phase3/genetic_map_chr{}_combined_b37.txt --input-vcf vcf/Phased/chr{}_strandaligned.filt_refallele_maf0.01_hwe_checked.vcf --input-ref ~/Murray/Bioinformatics/Reference_Files/Impute_ref_files/1000GP_Phase3/1000GP_Phase3_chr{}.hap ~/Murray/Bioinformatics/Reference_Files/Impute_ref_files/1000GP_Phase3/1000GP_Phase3_chr{}.legend ~/Murray/Bioinformatics/Reference_Files/Impute_ref_files/1000GP_Phase3/1000GP_Phase3.sample --output-log vcf/Phased/chr{}_strandaligned.filt_refallele_maf0.01_hwe_checked -T 12 ' ::: $(seq 1 22)

	parallel 'touch vcf/Phased/chr{}_strandaligned.filt_refallele_maf0.01_hwe_checked.snp.strand.exclude ' ::: $(seq 1 22)
	parallel -j 1 ' ~/Murray/src/shapeit.v2.r837/bin/shapeit -M ~/Murray/Bioinformatics/Reference_Files/Impute_ref_files/1000GP_Phase3/genetic_map_chr{}_combined_b37.txt --input-vcf vcf/Phased/chr{}_strandaligned.filt_refallele_maf0.01_hwe_checked.vcf --input-ref ~/Murray/Bioinformatics/Reference_Files/Impute_ref_files/1000GP_Phase3/1000GP_Phase3_chr{}.hap ~/Murray/Bioinformatics/Reference_Files/Impute_ref_files/1000GP_Phase3/1000GP_Phase3_chr{}.legend ~/Murray/Bioinformatics/Reference_Files/Impute_ref_files/1000GP_Phase3/1000GP_Phase3.sample --output-max vcf/Phased/chr{}_strandaligned.filt_refallele_maf0.01_hwe_checked --exclude-snp vcf/Phased/chr{}_strandaligned.filt_refallele_maf0.01_hwe_checked.snp.strand.exclude -T 12 ' ::: $(seq 1 22)

	parallel '~/Murray/src/shapeit.v2.r837/bin/shapeit -convert --input-haps vcf/Phased/chr{}_strandaligned.filt_refallele_maf0.01_hwe_checked --output-vcf vcf/Phased/chr{}_phased_strandaligned.filt_refallele_maf0.01_hwe_checked.vcf ' ::: $(seq 1 22)

	parallel 'bgzip {} ;tabix -p vcf {}.gz ' ::: $(ls vcf/Phased/*.vcf)

	bcftools concat -l -O z -o vcf/Phased/all_chr_phased_strandaligned.filt_refallele_maf0.01_hwe_checked.vcf.gz $(for i in $(seq 1 22); do echo vcf/Phased/chr${i}_phased_strandaligned.filt_refallele_maf0.01_hwe_checked.vcf.gz ; done | tr '\n' ' ') 

#python ~/Murray/Impute_pipeline/checkVCF.py -r ~/Murray/Bioinformatics/Reference_Files/FASTA/hs37d5/hs37d5.fa -o vcf_checking vcf/all_chr_strandaligned.filt_refallele_maf0.01_hwe_checked_strand_exclude.vcf
fi

if [ $phase == 'nophase' ]
then
	for i in $(seq 1 22); do plink2 --vcf vcf/all_chr_strandaligned.filt_refallele_maf0.01_hwe_checked.vcf --recode vcf --out vcf/chr${i}_strandaligned.filt_refallele_maf0.01_hwe_checked --a2-allele vcf/1kgAlleles.txt 3 2 --chr ${i} ; done

parallel -j 1 ' ~/Murray/src/shapeit.v2.r837/bin/shapeit -check -M ~/Murray/Bioinformatics/Reference_Files/Impute_ref_files/1000GP_Phase3/genetic_map_chr{}_combined_b37.txt --input-vcf vcf/chr{}_strandaligned.filt_refallele_maf0.01_hwe_checked.vcf --input-ref ~/Murray/Bioinformatics/Reference_Files/Impute_ref_files/1000GP_Phase3/1000GP_Phase3_chr{}.hap ~/Murray/Bioinformatics/Reference_Files/Impute_ref_files/1000GP_Phase3/1000GP_Phase3_chr{}.legend ~/Murray/Bioinformatics/Reference_Files/Impute_ref_files/1000GP_Phase3/1000GP_Phase3.sample --output-log vcf/chr{}_strandaligned.filt_refallele_maf0.01_hwe_checked -T 12 ' ::: $(seq 1 22)

cat vcf/chr*_strandaligned.filt_refallele_maf0.01_hwe_checked*strand | grep -v 'type' | cut -f4 | sort | uniq > vcf/strand_snps_exclude.txt

plink2 --vcf vcf/all_chr_strandaligned.filt_refallele_maf0.01_hwe_checked.vcf --recode vcf --out vcf/all_chr_strandaligned.filt_refallele_maf0.01_hwe_checked_strand_exclude --a2-allele vcf/1kgAlleles.txt 3 2 --exclude-snps $(cat vcf/strand_snps_exclude.txt) --maf 0.01 --missing

python ~/Murray/Impute_pipeline/checkVCF.py -r ~/Murray/Bioinformatics/Reference_Files/FASTA/hs37d5/hs37d5.fa -o vcf_checking vcf/all_chr_strandaligned.filt_refallele_maf0.01_hwe_checked_strand_exclude.vcf

fi

#filter on missingness >2% total
#plink2 --vcf vcf/all_chr_strandaligned.filt_refallele_maf0.001.vcf --a2-allele vcf/1kgAlleles.txt 3 2 --mind 0.02 --recode vcf --out vcf/all_chr_strandaligned.filt_refallele_maf0.001_mind0-02

#filter SNPs >3% geno missing
#plink2 --vcf vcf/all_chr_strandaligned.filt_refallele_maf0.001_mind0-02.vcf --a2-allele vcf/1kgAlleles.txt 3 2 --mind 0.02 --recode vcf --out vcf/all_chr_strandaligned.filt_refallele_maf0.001_mind0-02_geno0-03





