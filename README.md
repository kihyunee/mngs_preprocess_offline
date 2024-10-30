# mngs_preprocess_offline

## Installation in a linux system

```
git clone https://github.com/kihyunee/mngs_preprocess_offline.git
cd mngs_preprocess_offline/
source ./set_up_preprocessing_dependencies.sh
```

## Usage

### Use case: 
- You have basecalled fastq files generated after ONT sequencing and basecalling through MINKNOW.
- There will be a directory such as 'fastq_pass' that contains barcode-by-barcode subdirectories. Inside these barcode subdirectories, we expect to have fastq.gz files.   
- Now you will use this script to (1) pool the reads from each barcode into a single fastq.gz file, (2) apply quality filter and collect the passed reads, (3) remove host genome-aligned reads.

### Input:


### Run:
```
python mngs_ont_read_preprocess_from_basecall_directory.py --barcode_def <barcode definition tsv file> --bd <basecall fastq directory> --threads 4 --host_genome_mmi <your path to T2T-CHM13v2.0_genomic.mmi_mapont> --outdir <your path to pre-processed output files>
```
