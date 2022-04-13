#!/bin/bash

#BSUB -J genome_prep
#BSUB -W 10:00
#BSUB -n 8
#BSUB -R rusage[mem=64GB]
#BSUB -e %J.stderr
#BSUB -o %J.stdout 

module load singularity/3.6.0

SINGULARITY_ENV=/home/schmidth/singularity_images/mmarge.sif
MMARGE="singularity exec --bind /data/leslie/ $SINGULARITY_ENV MMARGE.pl"

DIR=inputs/

mkdir -p $DIR

SNPS=$DIR/129S1_SvImJ.mgp.v5.snps.dbSNP142.vcf
INDELS=$DIR/129S1_SvImJ.mgp.v5.indels.dbSNP142.normed.vcf

if [ ! -f $SNPS ]; then
    echo "Downloading and pre-processing latest SNP files for $(basename $SNPS)."
    echo "======================================================================"
    wget ftp://ftp-mouse.sanger.ac.uk/REL-1505-SNPs_Indels/strain_specific_vcfs/$(basename "$SNPS").gz -O "$SNPS".gz
    gunzip "$SNPS".gz
    grep "^#" $SNPS > "$SNPS".sorted
    grep -v "^#" $SNPS | sort -k1,1 -k2,2n --parallel=24 >> "$SNPS".sorted
fi

if [ ! -f $INDELS ]; then
    echo "Downloading and pre-processing latest INDEL files for $(basename $INDELS)."
    echo "======================================================================"
    wget ftp://ftp-mouse.sanger.ac.uk/REL-1505-SNPs_Indels/strain_specific_vcfs/$(basename "$INDELS").gz -O "$INDELS".gz
    gunzip "$INDELS".gz
    grep "^#" $INDELS > "$INDELS".sorted
    grep -v "^#" $INDELS | sort -k1,1 -k2,2n --parallel=24 >> "$INDELS".sorted
fi

REF_FTP=$(head $INDELS | PERL_BADLANG=0 perl -Xne 'print $1 if /reference=(ftp.*)/')
REF=$DIR/$(basename $REF_FTP)
if [ ! -f $REF ]; then
    echo "Installing reference genome from: $REF_FTP"
    echo "======================================================================"
    wget "$REF_FTP".gz -O $REF.gz
    gunzip $REF.gz
fi

echo "Splitting genome into chromosomes."
echo "======================================================================"

if [ ! -f scripts/faSplit ]; then
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v385/faSplit
    chmod u+x scripts/faSplit
fi

mkdir -p genomes
./scripts/faSplit byname $REF genomes/

for f in genomes/*; do 
    mv $f genomes/chr$(basename $f) 
done

echo "Preparing files with MMARGE"
echo "======================================================================"

$MMARGE prepare_files -core 8 -files "$SNPS".sorted "$INDELS.sorted" -genome genomes -dir . -genome_dir . 
