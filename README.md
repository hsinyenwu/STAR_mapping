### STAR_mapping
(Manual)[https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&ved=2ahUKEwjxpeSS36L_AhViFVkFHYAgApkQFnoECBAQAQ&url=https%3A%2F%2Fraw.githubusercontent.com%2Falexdobin%2FSTAR%2Fd14a0a992f94ba3a64c26dd08ac58e2b4ab134f3%2Fdoc%2FSTARmanual.pdf&usg=AOvVaw0_N6V0i7zVO0D7jpRPqpOX]

Here is an example for mapping one sample:
Step1 generating index:
```
#!/bin/bash -l
#SBATCH -n 12 --time=00:30:00 --mem=30gb
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
