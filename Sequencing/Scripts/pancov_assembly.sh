#!/bin/bash
#source conda to activate conda env
source ~/miniconda3/etc/profile.d/conda.sh
conda activate mpox

n=0
while getopts ":n:" option; do
        case "${option}" in
                n) n=$OPTARG;;
        esac
done

if [ $n == 0 ] ; then
        echo "You must specify the run name with -n."
        exit 1
fi

rm -rf raw_fastqc trimmed_fastqc trimmed intermediate results bamstats
mkdir raw_fastqc trimmed_fastqc trimmed intermediate results bamstats
mkdir intermediate/hostremoval intermediate/assemblies intermediate/genomecov

find $n -name "*.fastq.gz" | \
parallel -j 10 \
fastqc -t 2 -o raw_fastqc {}

find $n -name "*_L001_R1_001.fastq.gz"|\
sed "s/_L001_R1_001.fastq.gz//g;s/"$n"//g;s/\///g" > ./intermediate/samplesheet.txt

parallel \
'trim_galore --quality 20 --length 36 --cores 2 \
--fastqc --fastqc_args "-o trimmed_fastqc" \
-o ./trimmed --paired {1}/{2}_L001_R1_001.fastq.gz \
{1}/{2}_L001_R2_001.fastq.gz 2>./trimmed_fastqc/{2}_trim_garlore.log' \
::: $n :::: ./intermediate/samplesheet.txt

parallel -j 3 \
'bwa mem /mnt/2tbssd/human_ref/GCA_000001405.15_GRCh38_full_analysis_set.fna -t 8 \
./trimmed/{}_L001_R1_001_val_1.fq.gz ./trimmed/{}_L001_R2_001_val_2.fq.gz \
> ./intermediate/hostremoval/{}_hostmapped.sam' :::: ./intermediate/samplesheet.txt

parallel \
'samtools view -bS ./intermediate/hostremoval/{}_hostmapped.sam > ./intermediate/hostremoval/{}_hostmapped.bam' :::: ./intermediate/samplesheet.txt

parallel \
'samtools view -h -b -q 30 ./intermediate/hostremoval/{}_hostmapped.bam  -U ./intermediate/hostremoval/{}_mapqun30.bam -o ./intermediate/hostremoval/{}_garbage.bam' :::: ./intermediate/samplesheet.txt

parallel \
'samtools sort -n ./intermediate/hostremoval/{}_mapqun30.bam -o ./intermediate/hostremoval/{}_targetsorted.bam' :::: ./intermediate/samplesheet.txt

parallel \
'samtools fixmate ./intermediate/hostremoval/{}_targetsorted.bam ./intermediate/hostremoval/{}_targetsortedfixmate.bam' :::: ./intermediate/samplesheet.txt

parallel \
'bedtools bamtofastq -i ./intermediate/hostremoval/{}_targetsortedfixmate.bam \
-fq ./intermediate/hostremoval/{}_hostremoved_R1.fastq -fq2 ./intermediate/hostremoval/{}_hostremoved_R2.fastq' :::: ./intermediate/samplesheet.txt

gzip ./intermediate/hostremoval/*.fastq

parallel \
'megahit -1 ./intermediate/hostremoval/{}_hostremoved_R1.fastq.gz \
-2 ./intermediate/hostremoval/{}_hostremoved_R2.fastq.gz -o ./intermediate/assemblies/{}_assembly' :::: ./intermediate/samplesheet.txt

parallel \
'minimap2 -sr -a ~/dragen_20230404/dragen_20230404.mmi \
./intermediate/assemblies/{}_assembly/final.contigs.fa \
> ./intermediate/{}.sam'  :::: ./intermediate/samplesheet.txt

parallel \
'samtools view -bS ./intermediate/{}.sam | samtools sort - -o ./intermediate/{}.bam' :::: ./intermediate/samplesheet.txt

parallel \
'samtools index ./intermediate/{}.bam' :::: ./intermediate/samplesheet.txt

parallel \
"bedtools genomecov -d -ibam ./intermediate/{}.bam > ./intermediate/genomecov/{}_genomecov.txt" :::: ./intermediate/samplesheet.txt

parallel \
"cat ./intermediate/genomecov/{}_genomecov.txt |\
awk '{sum3[\$1] += \$3; count3[\$1]++}; END {for (id in sum3) { print id, sum3[id]/count3[id] } }' | \
awk '\$2!=0'| awk '{print \$1}' |  awk 'BEGIN {OFS = \",\"} {\$2 = \"{}\" \$2; print}' \
>> ./intermediate/best_refs.txt" :::: ./intermediate/samplesheet.txt

parallel \
"cat ./intermediate/genomecov/{}_genomecov.txt |\
awk '{sum3[\$1] += \$3; count3[\$1]++}; END {for (id in sum3) { print id, sum3[id]/count3[id] } }' | \
awk '\$2!=0'| awk 'BEGIN {OFS = \" \"} {\$3 = \"{}\" \$3; print}' | sed 's/ /,/g' \
>> ./intermediate/best_refs_hits.csv" :::: ./intermediate/samplesheet.txt

parallel --colsep "," \
'minimap2 -xsr -a ~/dragen_sep_20230404/{1}.mmi \
./intermediate/hostremoval/{2}_hostremoved_R1.fastq.gz  ./intermediate/hostremoval/{2}_hostremoved_R2.fastq.gz\
> ./intermediate/{2}_{1}_bestref.sam'  :::: ./intermediate/best_refs.txt

parallel --colsep "," \
'samtools view -bS ./intermediate/{2}_{1}_bestref.sam | samtools sort - -o ./intermediate/{2}_{1}_bestref.bam' :::: ./intermediate/best_refs.txt

parallel --colsep "," \
'samtools index ./intermediate/{2}_{1}_bestref.bam' :::: ./intermediate/best_refs.txt

parallel -j 10 --colsep "," \
'bcftools mpileup -f ~/dragen_sep_20230404/{1}.fasta ./intermediate/{2}_{1}_bestref.bam | \
bcftools call -c | vcfutils.pl vcf2fq > ./intermediate/{2}_{1}_bestref.consensus.fastq' :::: ./intermediate/best_refs.txt

parallel -j 10 --colsep "," \
'seqtk seq -aQ64 -q20 -n N ./intermediate/{2}_{1}_bestref.consensus.fastq > ./results/{2}_{1}_bestref.consensus.fasta'  :::: ./intermediate/best_refs.txt

parallel --colsep "," \
'qualimap bamqc -bam ./intermediate/{2}_{1}_bestref.bam -outfile {2}_{1}_qc.pdf -outdir ./results/{2}_{1}_bamqc -c' :::: ./intermediate/best_refs.txt

wait

cat ./intermediate/best_refs.txt | sed 's/_S.*//' > ./intermediate/best_refs_cutnames.txt

for filename in results/*;
do
mv "$filename" "$(echo "$filename" | sed -r 's/_S[0-9]{1,2}//')" ;
done

parallel --colsep "," \
"bedtools genomecov -d -ibam ./intermediate/{2}_{1}_bestref.bam |\
sed '1iReference_Name\tPosition\tCoverage' \
> ./intermediate/genomecov/{2}_{1}_depth.txt" :::: ./intermediate/best_refs.txt

for consensus in results/*_bestref.consensus.fasta
do
base=$(basename $consensus "_bestref.consensus.fasta")
sed -i "s/^>.*/>${base}/g" "$consensus"
done

conda deactivate
for consensus in results/*consensus.fasta
do
out=results/$(basename $consensus "_bestref.consensus.fasta" )
       blastn -db nt -taxidlist ~/scripts/11118.txids \
       -query $consensus -outfmt "10 delim=@ \
       qseqid sallseqid stitle pident length qcovs mismatch \
       evalue" -num_threads 20 -out ${out}_blastn.csv
done


head -1 ./results/*_blastn.csv |  grep -v "==>"| sed '/^$/d' | awk 'BEGIN{ OFS=FS="@" }{ sub(/gi\|.*\|gb\|/, "", $2) }1' |\
awk 'BEGIN{ OFS=FS="@" }{ sub(/gi\|.*\|emb\|/, "", $2) }1' | awk 'BEGIN{ OFS=FS="@" }{ sub(/gi\|.*\|dbj\|/, "", $2) }1' |\
awk 'BEGIN{ OFS=FS="@" }{ sub(/\|/, "", $2) }1' |  awk 'BEGIN{ OFS=FS="@" }{ sub(/Severe acute respiratory syndrome coronavirus 2 isolate /, "", $3) }1' |\
awk 'BEGIN{ OFS=FS="@" }{ sub(/Severe acute respiratory syndrome coronavirus 2 genome assembly, ter/, "SARS-CoV2 ", $3) }1' | awk 'BEGIN{FS=OFS="@"} {$3=substr($3,1,30)} 1'\
>> ./results/top_blast_hits.csv

sed -i '1iSampleID_BestRef@BlastMatchID@Blast_Match_Name@ %_ID@Length@%Qcov@Mismatch@E-value' ./results/top_blast_hits.csv


conda activate pancov
parallel --colsep "," \
"cat ./results/{2}_{1}_bamqc/genome_results.txt |\
sed -n '/^>>>>>>> Globals.*/,/std coverage.*/{/^>/d;/There/d;/^$/d;s/^[ \t]*//;s/ = /;/g;p;}' |\
sed 1i'header;value' > ./bamstats/{2}-{1}-bamstats.txt" :::: ./intermediate/best_refs_cutnames.txt 
