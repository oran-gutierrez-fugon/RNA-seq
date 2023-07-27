#For forward strand in Undifferentiated sample 1


#First, we subsampled the genome-wide alignments to the locus of interest and only reads in forward strand (paternal):
#(for forward strand change to -F16)

module load samtools

samtools view -f16 -Sb N3A_S82_L004_Aligned.sortedByCoord.out.bam chr15:22351127-27474172 > N3A_S82_L004_Aligned.sortedByCoord.out.ube3a.bam


#Then, we added read groups:
module load picard-tools/2.26.11

picard AddOrReplaceReadGroups I=N3A_S82_L004_Aligned.sortedByCoord.out.ube3a.bam O=N3A_S82_L004_Aligned.sortedByCoord.out.ube3a.rg.bam SO=coordinate RGID=N3A_S82_L004 RGLB=library RGPL=illumina RGPU=N3A_S82_L004 RGSM=N3A_S82_L004

#We marked duplicates:

picard MarkDuplicates I=N3A_S82_L004_Aligned.sortedByCoord.out.ube3a.rg.bam O=N3A_S82_L004_Aligned.sortedByCoord.out.ube3a.rg.dedupped.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics

#We reformatted RNA-seq alignments for variant calling:

module load gatk/4.2.0.0

gatk SplitNCigarReads -R /share/lasallelab/Oran/dovetail/refgenomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -I N3A_S82_L004_Aligned.sortedByCoord.out.ube3a.rg.dedupped.bam -O N3A_S82_L004_Aligned.sortedByCoord.out.ube3a.rg.dedupped.split.bam 

#Calling variants using GATK pipeline:

gatk HaplotypeCaller \
-R /share/lasallelab/Oran/dovetail/refgenomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
-I N3A_S82_L004_Aligned.sortedByCoord.out.ube3a.rg.dedupped.split.bam \
--intervals chr15:22351127-27474172 \
-O Variants_called.vcf

#Download Variant_called.vcf and open in text editor and search for the coordinates in Hg38.

#Zooming in into UBE3A locus:

#module load htslib/1.10.2
#bgzip -c Variants_called.vcf > Variants_called.vcf.gz
#tabix Variants_called.vcf.gz

#module load bcftools/1.15
#bcftools view -r chr15:25333728-25439024 Variants_called.vcf.gz