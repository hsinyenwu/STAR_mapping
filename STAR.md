## STAR_mapping
Read the manual if possible. It will help you to think about what to do for the multimapping reads.    
[Manual for STAR 2.7.9a](https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&ved=2ahUKEwjxpeSS36L_AhViFVkFHYAgApkQFnoECBAQAQ&url=https%3A%2F%2Fraw.githubusercontent.com%2Falexdobin%2FSTAR%2Fd14a0a992f94ba3a64c26dd08ac58e2b4ab134f3%2Fdoc%2FSTARmanual.pdf&usg=AOvVaw0_N6V0i7zVO0D7jpRPqpOX)  
[Google Group Help Page](https://groups.google.com/g/rna-star)  

Here is an example for mapping one sample (i.e., sample D123567):  
### Step1 generating index:
```
#!/bin/bash -l
#SBATCH -n 12 -N 1 --time=00:30:00 --mem=30gb
#SBATCH -J STAR_index
#SBATCH -o STAR_RNA_Ribo_index.out
#SBATCH -e STAR_RNA_Ribo_index.err

############################################
#xxx is the path
FASTA=/xxx/TAIR10_chr_all_2.fas
GTF=/xxx/Araport11_20181206.gtf

##########
module purge
module load GCC/8.3.0
module load STAR/2.7.9a

###############
## RNA index ##
###############
newINDEX=/xxx/STAR1_Index/RNA

###########################################
# Generate index with Araport11

mkdir -p $newINDEX
cd $newINDEX

echo "Generate index with  for STAR"
STAR --runThreadN 12 \
--runMode genomeGenerate \
--genomeDir $newINDEX \
--genomeFastaFiles $FASTA \
--sjdbGTFfile $GTF \
--sjdbOverhang 99 \

#--sjdbOverhang 99 beacuse of the reads are paired end 100 (max length reads), this value is 100-1=99

###############
## Ribo index ##
###############
# Ribo-seq is shorter. I estimate the maximum read lenth is ~35
newINDEX=/xxx/STAR1_Index/Ribo

###########################################
# Generate index with Araport11
mkdir -p $newINDEX
cd $newINDEX

echo "Generate index with  for STAR"
STAR --runThreadN 12 \
--runMode genomeGenerate \
--genomeDir $newINDEX \
--genomeFastaFiles $FASTA \
--sjdbGTFfile $GTF \
--sjdbOverhang 34 \ 

############################################
# output details
echo "details of the job and queue"
echo ---------

#####################
scontrol show job ${SLURM_JOBID}
```

### Step2A mapping RNA-seq:
```
#!/bin/bash -l
#SBATCH -n 10 --time=2:00:00 --mem=50gb
#SBATCH -J RNA_mapping
#SBATCH --output=STAR1-RNA-mapping.out

############################################
INPUT=/xxx/RNA_preprocessed
OUTPUT=/xxx/RNA_CTRL_merged
FASTA=/xxx/TAIR10_chr_all_2.fas
starIndex1=/xxx/STAR1_Index/RNA
GTF=/xxx/Araport11_20181206.gtf
##########

mkdir -p $OUTPUT
cd $OUTPUT

##########
module purge
module load GCC/8.3.0
module load STAR/2.7.9a

#1st pass STAR alignment        
STAR --runThreadN 10 \
--genomeDir $starIndex1 \
--readFilesCommand zcat \
--readFilesIn $INPUT/D123567.r1.fastq.gz $INPUT/D123567.r2.fastq.gz \
--alignIntronMax 5000 \
--alignIntronMin 15 \
--outFilterMismatchNmax 2 \
--outFilterMultimapNmax 20 \
--outFilterType BySJout \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 2 \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM \
--outSAMmultNmax 1 \
--outMultimapperOrder Random \
--outFileNamePrefix "star_RNA_" \

############################################
# output details
echo "details of the job and queue"
echo ---------
scontrol show job ${SLURM_JOBID}
```
### Step2B mapping Ribo-seq:
```
#!/bin/bash -l
#SBATCH -n 10 --time=2:00:00 --mem=50gb
#SBATCH -J Ribo_mapping
#SBATCH --output=STAR1-Ribo-mapping.out

############################################
INPUT=/mnt/home/larrywu/CTRL_arabidopsis/data/ribo_preprocessed
OUTPUT=/mnt/home/larrywu/Kyle_CTRL/analysis_Araport11_1st+2nd_Riboseq/STAR1/Ribo_CTRL_merged
FASTA=/mnt/research/riboplant/Reference/TAIR10_chr_all_2.fas
starIndex1=/mnt/home/larrywu/Kyle_CTRL/analysis_Araport11_1st+2nd_Riboseq/STAR1_Index/Ribo
GTF=/mnt/home/larrywu/ABA/20181206/reference/Araport11_20181206.gtf
##########

mkdir -p $OUTPUT
cd $OUTPUT

##########
module purge
module load GCC/8.3.0
module load STAR/2.7.9a

STAR --runThreadN 10 \
--genomeDir $starIndex1 \
--readFilesCommand zcat \
--readFilesIn $INPUT/D123567.noContam4.fastq.gz \
--alignIntronMax 5000 \
--alignIntronMin 15 \
--outFilterMismatchNmax 1 \
--outFilterMultimapNmax 20 \
--outFilterType BySJout \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 2 \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM \
--outSAMmultNmax 1 \
--outMultimapperOrder Random \
--outFileNamePrefix "star_ribo" \

###############################################
# output details
echo "details of the job and queue"
echo ---------

scontrol show job ${SLURM_JOBID}
```
### Convert bam to fastq (for Kallisto)
```
module purge
module load GCC/10.2.0
module load BEDTools/2.30.0

bamToFastq -i /xxx/$Input.bam -fq /xxx/$Input.r1.fastq -fq2 /xxx/$Input.r2.fastq
```
