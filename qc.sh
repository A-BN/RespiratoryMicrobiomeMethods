#!/bin/bash

# antoine.bridier-nahmias@inserm.fr
# 26/11/2018

cd `dirname $0`

fastqc="/home/abn/Documents/tools/FastQC/fastqc"
# Je suis un gentil commentaire.

fastqa=`find ../data/reads/ -name "*fastq.gz"` # Je suis un beau pluriel latin

if [ ! -d ../data/reads/qc ];  then mkdir ../data/reads/qc; fi
for fastqus in ${fastqa}; do
	echo -e "Treating:  `basename $fastqus`"
	$fastqc \
		--outdir "../data/reads/qc" \
		--threads 4 \
		$fastqus

done

multiqc --interactive --outdir "../data/reads/qc/multiqc" "../data/reads/qc"
