#!/bin/bash
#source conda to activate conda env
source ~/miniconda3/etc/profile.d/conda.sh
conda activate bat-alpha-cov

n=0
while getopts ":n:" option; do
        case "${option}" in
                n) n=$OPTARG;;
        esac
done

if [ $n == 0 ] ; then
        echo "You must specify the folder name with -n."
        exit 1
fi


mkdir ebcov_pipeline_output fastqc_trimmed \
trimmed_fastqc trimmed intermediate results raw_fastqc bamstats
mkdir intermediate/genomecov

# make sample sheet
find $n -name "*_L001_R1_001.fastq.gz"| \
sed "s/_L001_R1_001.fastq.gz//g;s/"$n"//g;s/\///g" \
> ./intermediate/samplesheet.txt

find $n -name "*.fastq.gz" | \
parallel -j 10 \
fastqc -t 2 -o raw_fastqc {}

multiqc ./raw_fastqc/. --outdir raw_fastqc --filename raw_fastqc_multiqc_report

parallel \
'trim_galore --quality 20 --length 36 --cores 2 \
--fastqc --fastqc_args "-o trimmed_fastqc" \
-o trimmed --paired {1}/{2}_L001_R1_001.fastq.gz \
./{1}/{2}_L001_R2_001.fastq.gz' \
::: $n :::: ./intermediate/samplesheet.txt

cat intermediate/samplesheet.txt \
| awk -F '[_]' ' { print $1"_"$2,$1}' OFS=',' \
> ./intermediate/samplesheet_ids.txt

parallel \
'minimap2 -sr -a /mnt/2tbssd/ebcov_ref/ref.mmi \
./trimmed/{}_L001_R1_001_val_1.fq.gz ./trimmed/{}_L001_R2_001_val_2.fq.gz \
> ./intermediate/{}.sam'  :::: ./intermediate/samplesheet.txt

parallel --colsep "," \
'samtools view -bS ./intermediate/{}.sam | samtools sort - -o ./intermediate/{}.bam' :::: ./intermediate/samplesheet.txt

parallel --colsep "," \
'samtools index ./intermediate/{}.bam' :::: ./intermediate/samplesheet.txt

parallel \
'ivar trim -e -i ./intermediate/{}.bam -b /mnt/2tbssd/ebcov_ref/Ebcov.primer.bed -p {}_trimmed' :::: :::: ./intermediate/samplesheet.txt

parallel 'mv {}_trimmed.bam ./intermediate' :::: ./intermediate/samplesheet.txt

parallel \
'samtools sort ./intermediate/{}_trimmed.bam -o ./intermediate/{}_trimmedsorted.bam' :::: ./intermediate/samplesheet.txt

parallel \
'qualimap bamqc -bam ./intermediate/{}.bam -outfile {}_qc.pdf -outdir ./results/{}_bamqc -c' :::: ./intermediate/samplesheet.txt

parallel \
"bedtools genomecov -d -ibam ./intermediate/{}_trimmedsorted.bam|\
sed '1iReference_Name\tPosition\tCoverage' \
> ./intermediate/genomecov/{}_depth.txt" :::: ./intermediate/samplesheet.txt

parallel \
"cat ./results/{}_bamqc/genome_results.txt |\
sed -n '/^>>>>>>> Globals.*/,/std coverage.*/{/^>/d;/There/d;/^$/d;s/^[ \t]*//;s/ = /;/g;p;}' |\
sed 1i'header;value' > ./bamstats/{}-bamstats.txt" :::: ./intermediate/samplesheet_removed.txt

parallel -j 10 \
'bcftools mpileup -f /mnt/2tbssd/ebcov_ref/OL415262.1.fasta ./intermediate/{}_trimmedsorted.bam | \
bcftools call -c | vcfutils.pl vcf2fq > ./intermediate/{}.consensus.fastq' :::: ./intermediate/samplesheet.txt

 
parallel -j 10 \
'seqtk seq -aQ64 -q20 -n N ./intermediate/{}.consensus.fastq > ./results/{}.consensus.fasta'  :::: ./intermediate/samplesheet.txt

for file in ./results/*.fasta;
do
    sed -i "/^>/s/\$/ ${file%%.*}/" "$file"

done
