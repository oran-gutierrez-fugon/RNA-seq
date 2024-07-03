
#This env was named wrong, for actual crossmap env conda activate /share/lasallelab/Oran/dovetail/luhmes/RNAseq/share/lasallelab/Oran/dovetail/luhmes/RNAseq/crossmap
#For STAR and TrimGalore first part of RNAseq pipeline (Ben Laufer Git hub)
conda activate crossmap
#at /share/lasallelab/Oran/miniconda3/crossmap

#make sure you run the same STAR version for next section too or will error out.

mkdir -p GenomeDir


STAR --runThreadN 30 --runMode genomeGenerate genomeDir -/share/lasallelab/Oran/dovetail/luhmes/RNAseq/star_150 --genomeFastaFiles /share/lasallelab/Oran/dovetail/luhmes/RNAseq/star_150/Homo_sapiens.GRCh37.dna_sm.primary_assembly.fa --sjdbGTFfile /share/lasallelab/Oran/dovetail/luhmes/RNAseq/star_150/Homo_sapiens.GRCh37.87.gtf --sjdbOverhang 150


#STAR --runThreadN 40 --runMode genomeGenerate genomeDir -/share/lasallelab/Oran/dovetail/luhmes/RNAseq/star_150 --genomeFastaFiles /share/lasallelab/Oran/dovetail/refgenomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna --sjdbGTFfile /share/lasallelab/Oran/dovetail/refgenomes/hg38.refGene.gtf --sjdbOverhang 150




#Two ways to use awk to get read length
awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' N3A_S82_L004_R2_001.fastq
awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count}' N3A_S82_L004_R2_001.fastq




#For creating and viewing some of the tracks
/share/lasallelab/Oran/dovetail/luhmes/RNAseq/star_150


track type=bam name=RNA_LUHMES_Neurons_Tagged description="RNA_LUHMES_Neurons_Tagged"bigDataUrl="https://bioshare.bioinformatics.ucdavis.edu/bioshare/download/4jo97m7umuy7a9a/RNA-seq/hg19/Neurons-RNA-cat_REsorted_hg19_tagged.bam"

track name=reads description="LUHMES" type=bam bigDataUrl="https://bioshare.bioinformatics.ucdavis.edu/bioshare/download/hbl76b4wa1i6dvm/phased_possorted_bam.roi.bam"