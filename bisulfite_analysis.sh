#--------------------------
# Define folder paths
#--------------------------
sample=$2
ref="ref/$4"    

#---------------------------------------------------------
### Trim raw data ###
echo $(date) Start trimming for $sample >> script.progress.info
#-----------------------------------------  ----------------
mkdir raw/trimmed
$trim_galore \
    --paired \
    --fastqc \
    -q 30 \
    --path_to_cutadapt $cutadapt \
    -j 4 \
    -o raw/trimmed \
    raw/$sample\_1.fq.gz \
    raw/$sample\_2.fq.gz 

mkdir -p results/$sample
cd results/$sample

#---------------------------------------------------------
### Align trimmed reads ###
echo $(date) Start bismark for $sample >> script.progress.info
#---------------------------------------------------------
mkdir $logs
bismark $ref \
    -1 raw/trimmed/$sample\_1_val_1.fq.gz \
    -2 raw/trimmed/$sample\_2_val_2.fq.gz \
    --nucleotide_coverage \
    --temp_dir results/$sample/$sample\_tmp/ \
    --path_to_bowtie2 $bowtie2 \
    --samtools_path $samtools \
    --multicore 6 \
    -N=1 \
    -D=20 \
    -R=3 \
    --score_min 'L,0,-0.2' \
    --non_directional\
    > >(tee $logs/bismark_$sample\.err) 2> >(tee $logs/bismark_$sample\.log >&2)

#### QC alignment maps ###
samtools sort results/$sample/$sample\_1_val_1_bismark_bt2_pe.bam -O bam  -@ 10 > results/$sample/$sample\_1_val_1_bismark_bt2_pe.sorted.bam
qualimap bamqc -bam results/$sample/$sample\_1_val_1_bismark_bt2_pe.sorted.bam

# ------------------------------------------------------------
### Extract methylation levels ###
echo $(date) Start bismark_me for $sample >> script.progress.info
# ------------------------------------------------------------
$bismark_me -p \
    --bedGraph \
    --buffer_size 10G \
    --scaffolds \
    --counts \
    --cytosine_report \
    --multicore 6 \
    --output ./ \
    --genome_folder $ref \
    results/$sample/$sample\_1_val_1_bismark_bt2_pe.bam \
    > >(tee $logs/bismark_$sample\-methyl-extractor.err) 2> >(tee $logs/bismark_$sample\-methyl-extractor.log >&2)

#-------------------------------------------------------------
### Prepare summary report ###
echo $(date) Start of Bismark report preparation for $sample >> script.progress.info
#-------------------------------------------------------------
$bismark2report
$bismark2summary
#-------------------------------------------------------------
echo $(date) Pipeline ends now for $sample >> script.progress.info
#-------------------------------------------------------------
