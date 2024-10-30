#!/bin/bash

base_dir=$(pwd)

install_dir="${base_dir}/install"
mkdir ${install_dir}

db_dir="${base_dir}/preprocess_db"
mkdir ${db_dir}

# install minimap
cd ${install_dir}
curl -L https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28_x64-linux.tar.bz2 | tar -jxvf -
echo "\n\n### Minimap installed"
echo "\n### Minimap installation check (1) version"
${install_dir}/minimap2-2.28_x64-linux/minimap2 --version
echo "\n### Minimap installation check (2) help"
${install_dir}/minimap2-2.28_x64-linux/minimap2 -h

# install nanoq
cd ${install_dir}
VERSION=0.10.0
RELEASE=nanoq-${VERSION}-x86_64-unknown-linux-musl.tar.gz
wget https://github.com/esteinig/nanoq/releases/download/${VERSION}/${RELEASE}
tar xf nanoq-${VERSION}-x86_64-unknown-linux-musl.tar.gz
echo "\n\n### Nanoq installed"
echo "\n### Nanoq installation check (1) version"
${install_dir}/nanoq-${VERSION}-x86_64-unknown-linux-musl/nanoq --version
echo "\n### Nanoq installation check (2) help"
${install_dir}/nanoq-${VERSION}-x86_64-unknown-linux-musl/nanoq -h

# install pigz
cd ${install_dir}
wget https://zlib.net/pigz/pigz-2.8.tar.gz
tar -xzf pigz-2.8.tar.gz
cd pigz-2.8/
make
echo "\n\n### Pigz installed"
echo "\n### Pigz installation check (1) version"
${install_dir}/pigz-2.8/pigz --version
echo "\n### Pigz installation check (2) help"
${install_dir}/pigz-2.8/pigz -h

# download database
echo "\n\n### Downloading human genome reference"
cd ${base_dir}/preprocess_db
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz
gzip -d GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz
mv GCF_009914755.1_T2T-CHM13v2.0_genomic.fna T2T-CHM13v2.0_genomic.fna

echo "\n\n### Indexing human genome reference"
${install_dir}/minimap2-2.28_x64-linux/minimap2 -t 4 -x map-ont -d T2T-CHM13v2.0_genomic.mmi_mapont T2T-CHM13v2.0_genomic.fna



cd ${base_dir}

