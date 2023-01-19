## Resistance gene identification in Metagenomes using CARD/RGI
## Conda environment: rgi
## RGI can be run in two ways:
## 1. * metagenomic contigs are generated (e.g. using MEGAHIT)
##    *
## The first can create a rudimentary heatmap, for presence (perfect/strict) and absence

rgi load --card_json ../../CARD/card.json --local
nohup megahit -1 AFY01_1.fastq.gz -2 AFY01_2.fastq.gz -o Megahit_AFY01 -t 24 > megahit.nhp &
cd Megahit_AFY01
rgi main -n 16 -i final.contigs.fa -o 1FY01_MH --split_prodigal_jobs --input_type contig --local --clean
nohup rgi main -n 16 -i final.contigs.fa -o 1FY01_MH --split_prodigal_jobs --input_type contig --local --clean > rgimain.nhp &

## 2. direct from metagenomic reads, preferable, if 
rgi bwt
sid=11
nohup rgi bwt --read_one AFY${sid}_1.fastq.gz --read_two AFY${sid}_2.fastq.gz --aligner bowtie2 --output_file AFY${sid} --threads 10 --local > rgi_bwt_afy${sid}.nhp &

## TODO: run metagenomic analysis (2.) with WildCARD
dbdir=/mnt/Drive1/ahenschel/SambaShare/rgi
rgi load --wildcard_annotation ${dbdir}wildcard_database_v3.2.2.fasta --wildcard_index ${dbdir}wildcard/index-for-model-sequences.txt --card_annotation ${dbdir}card_database_v3.2.2.fasta

## execute one by one
rgi load --wildcard_annotation ${dbdir}wildcard_database_v3.2.2.fasta --wildcard_index ${dbdir}wildcard/index-for-model-sequences.txt --card_annotation ${dbdir}card_database_v3.2.2.fasta
for sid in `ls AFY??_1.fastq.gz|cut -c 4-5`;
do 
    nohup rgi bwt --read_one AFY${sid}_1.fastq.gz --read_two AFY${sid}_2.fastq.gz --aligner bowtie2 --output_file WildAFY${sid} --threads 10 --include_wildcard --local > rgi_bwt_afy${sid}.nhp 
done
