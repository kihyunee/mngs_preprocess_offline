#!/bin/bash

base_dir=$(pwd)

install_dir="${base_dir}/install"
mkdir -p ${install_dir}

db_dir="${base_dir}/preprocess_db"
mkdir -p ${db_dir}

# Install minimap
cd ${install_dir}
curl -L https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28_x64-linux.tar.bz2 | tar -jxvf -
echo -e "\n\n### Minimap installed"
echo -e "\n### Minimap installation check (1) version"
${install_dir}/minimap2-2.28_x64-linux/minimap2 --version
echo -e "\n### Minimap installation check (2) help"
${install_dir}/minimap2-2.28_x64-linux/minimap2 -h

# Install nanoq
cd ${install_dir}
VERSION=0.10.0
RELEASE=nanoq-${VERSION}-x86_64-unknown-linux-musl.tar.gz
wget https://github.com/esteinig/nanoq/releases/download/${VERSION}/${RELEASE}
tar xf ${RELEASE}
echo -e "\n\n### Nanoq installed"
echo -e "\n### Nanoq installation check (1) version"
${install_dir}/nanoq-${VERSION}-x86_64-unknown-linux-musl/nanoq --version
echo -e "\n### Nanoq installation check (2) help"
${install_dir}/nanoq-${VERSION}-x86_64-unknown-linux-musl/nanoq -h

# Install pigz
cd ${install_dir}
wget https://zlib.net/pigz/pigz-2.8.tar.gz
tar -xzf pigz-2.8.tar.gz
cd pigz-2.8/
make
echo -e "\n\n### Pigz installed"
echo -e "\n### Pigz installation check (1) version"
${install_dir}/pigz-2.8/pigz --version
echo -e "\n### Pigz installation check (2) help"
${install_dir}/pigz-2.8/pigz -h

# Download database
echo -e "\n\n### Downloading human genome reference"
cd ${db_dir}
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz
gzip -d GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz
mv GCF_009914755.1_T2T-CHM13v2.0_genomic.fna T2T-CHM13v2.0_genomic.fna

echo -e "\n\n### Indexing human genome reference"
${install_dir}/minimap2-2.28_x64-linux/minimap2 -t 4 -x map-ont -d T2T-CHM13v2.0_genomic.mmi_mapont T2T-CHM13v2.0_genomic.fna

# Add binaries to PATH
cd ${base_dir}
mkdir -p ${install_dir}/bin
cp ${install_dir}/minimap2-2.28_x64-linux/minimap2 ${install_dir}/bin/
cp ${install_dir}/nanoq-${VERSION}-x86_64-unknown-linux-musl/nanoq ${install_dir}/bin/
cp ${install_dir}/pigz-2.8/pigz ${install_dir}/bin/

# Ensure PATH is updated in the current shell and future sessions
export PATH="${install_dir}/bin:$PATH"
if ! grep -q "${install_dir}/bin" ~/.bashrc; then
    echo "export PATH=\"${install_dir}/bin:\$PATH\"" >> ~/.bashrc
fi
if ! grep -q "${install_dir}/bin" ~/.bash_profile; then
    echo "export PATH=\"${install_dir}/bin:\$PATH\"" >> ~/.bash_profile
fi

source ~/.bashrc