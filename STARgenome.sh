#make sure you run the same STAR version for next section too or will error out.

mkdir -p GenomeDir


STAR --runThreadN 40 --runMode genomeGenerate genomeDir -/star_150 --genomeFastaFiles /share/lasallelab/Oran/dovetail/refgenomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna --sjdbGTFfile /share/lasallelab/Oran/dovetail/refgenomes/hg38.refGene.gtf --sjdbOverhang 150
#STAR --runThreadN 40 --runMode genomeGenerate genomeDir -/share/lasallelab/Oran/dovetail/luhmes/RNAseq/star_150 --genomeFastaFiles /share/lasallelab/Oran/dovetail/refgenomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna --sjdbGTFfile /share/lasallelab/Oran/dovetail/refgenomes/hg38.refGene.gtf --sjdbOverhang 150




#Two ways to use awk to get read length
awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' N3A_S82_L004_R2_001.fastq
awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count}' N3A_S82_L004_R2_001.fastq