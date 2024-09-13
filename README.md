# Resistance gene identification in Metagenomes using CARD/RGI
## 0. Minimal demo usage for simple RGI output visualization
* clone this repository
* change into directory
* install jupyter with python libraries pandas, numpy and seaborn. This can be quickly done with anaconda installed and
`conda install -y jupyter pandas numpy scipy seaborn`
or
`conda env create -f rgimini.yaml` (assumunig you have done the RGI jobs already, otherwise use rgi.yaml, which includes the RGI dependencies, however double check with RGI installation instructions)
* start a jupyter notebook
* open rgivisual_demo.ipynb
* follow the instructions in the notebook

RGI can be run in two ways for Metagenomes
## 1. metagenomic contigs are generated (e.g. using MEGAHIT)

`rgi load --card_json ../../CARD/card.json --local`

`nohup megahit -1 AFY01_1.fastq.gz -2 AFY01_2.fastq.gz -o Megahit_AFY01 -t 24 > megahit.nhp &`
`cd Megahit_AFY01`
`rgi main -n 16 -i final.contigs.fa -o 1FY01_MH --split_prodigal_jobs --input_type contig --local --clean`
`nohup rgi main -n 16 -i final.contigs.fa -o 1FY01_MH --split_prodigal_jobs --input_type contig --local --clean > rgimain.nhp &`

Used Conda environment: rgi (see rgi.yaml and conda environment creation from yaml file)

## 2. direct from metagenomic reads
Example run
`rgi bwt`
`sid=11`
`nohup rgi bwt --read_one AFY${sid}_1.fastq.gz --read_two AFY${sid}_2.fastq.gz --aligner bowtie2 --output_file AFY${sid} --threads 10 --local > rgi_bwt_afy${sid}.nhp &`

## TODO: run metagenomic analysis (2.) with WildCARD
`dbdir=<rgi_database_directory>`
`rgi load --wildcard_annotation ${dbdir}wildcard_database_v3.2.2.fasta --wildcard_index ${dbdir}wildcard/index-for-model-sequences.txt --card_annotation ${dbdir}card_database_v3.2.2.fasta`

## execute one by one
`rgi load --wildcard_annotation ${dbdir}wildcard_database_v3.2.2.fasta --wildcard_index ${dbdir}wildcard/index-for-model-sequences.txt --card_annotation ${dbdir}card_database_v3.2.2.fasta`

```for sid in `ls AFY??_1.fastq.gz|cut -c 4-5`;
do 
    nohup rgi bwt --read_one AFY${sid}_1.fastq.gz --read_two AFY${sid}_2.fastq.gz --aligner bowtie2 --output_file WildAFY${sid} --threads 10 --include_wildcard --local > rgi_bwt_afy${sid}.nhp 
done```
