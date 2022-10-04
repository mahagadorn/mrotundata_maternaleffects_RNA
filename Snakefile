"""
Author: Karen M. Kapheim
Date: 04 April 2021
Aim: Process mRNA sequences from Megachile rotundata maternal RNA project
Run: From Frisco node, `./run_workflow.sh`; Make sure it is executable
"""

#### Set working directory
workdir: "/uufs/chpc.utah.edu/common/home/kapheim-group1/mrot_maternalRNA/alignment_quantification_mRNAs/"

#### Get sample names

#This os library let's you read out filenames of directories

import os

#we will populate this list with the unique prefixes of our paired end reads

seq_names = []

# Get .fastq file names. I'm only going to get read one names, because read two have the same prefix

for f_name in os.listdir('/uufs/chpc.utah.edu/common/home/kapheim-group1/mrot_maternalRNA/alignment_quantification_mRNAs/reads'):
        if f_name.endswith('_R1_001.fastq.gz'):
                #just getting the unique prefixes of the fastq files by stripping the righthand side of the filenames
                f_name = f_name.rstrip('_R1_001.fastq.gz')
                #building up the list
                seq_names.append(f_name)
seqs=seq_names
print(len(seqs))

# Alternative method:
#
# sampleIDs=[
#    "20129_ATATCTCG-ACTAAGAT_L00M",
#    "202230_AATGCCTC-TGGATCGA_L00M"
#]


#### Get to work

rule all:
    input:
        expand("/uufs/chpc.utah.edu/common/home/kapheim-group1/mrot_maternalRNA/alignment_quantification_mRNAs/trimmed/{sample}_R1.fastq.gz", sample=seq_names),
        expand("/uufs/chpc.utah.edu/common/home/kapheim-group1/mrot_maternalRNA/alignment_quantification_mRNAs/trimmed/{sample}_R2.fastq.gz", sample=seq_names),
        expand("/uufs/chpc.utah.edu/common/home/kapheim-group1/mrot_maternalRNA/alignment_quantification_mRNAs/trimmed/{sample}_R1.unpaired.fastq.gz", sample=seq_names),
        expand("/uufs/chpc.utah.edu/common/home/kapheim-group1/mrot_maternalRNA/alignment_quantification_mRNAs/trimmed/{sample}_R2.unpaired.fastq.gz", sample=seq_names),
        "/uufs/chpc.utah.edu/common/home/kapheim-group1/mrot_maternalRNA/alignment_quantification_mRNAs/genome_index",
        expand("/uufs/chpc.utah.edu/common/home/kapheim-group1/mrot_maternalRNA/alignment_quantification_mRNAs/star_aligned/{sample}_Aligned.out.sam", sample=seq_names),
        expand("/uufs/chpc.utah.edu/common/home/kapheim-group1/mrot_maternalRNA/alignment_quantification_mRNAs/htseq/counts_{sample}.txt", sample=seq_names)

rule trimmomatic_pe:
    input:
        trimm_path = "/uufs/chpc.utah.edu/sys/installdir/trimmomatic/Trimmomatic-0.39/",
        r1="/uufs/chpc.utah.edu/common/home/kapheim-group1/mrot_maternalRNA/alignment_quantification_mRNAs/reads/{sample}_R1_001.fastq.gz",
        r2="/uufs/chpc.utah.edu/common/home/kapheim-group1/mrot_maternalRNA/alignment_quantification_mRNAs/reads/{sample}_R2_001.fastq.gz"
    output:
        r1_p="/uufs/chpc.utah.edu/common/home/kapheim-group1/mrot_maternalRNA/alignment_quantification_mRNAs/trimmed/{sample}_R1.fastq.gz",
        r2_p="/uufs/chpc.utah.edu/common/home/kapheim-group1/mrot_maternalRNA/alignment_quantification_mRNAs/trimmed/{sample}_R2.fastq.gz",
        # reads where trimming entirely removed the mate
        r1_up="/uufs/chpc.utah.edu/common/home/kapheim-group1/mrot_maternalRNA/alignment_quantification_mRNAs/trimmed/{sample}_R1.unpaired.fastq.gz",
        r2_up="/uufs/chpc.utah.edu/common/home/kapheim-group1/mrot_maternalRNA/alignment_quantification_mRNAs/trimmed/{sample}_R2.unpaired.fastq.gz"
    threads:
        8
    shell:"""
        module load trimmomatic/0.39
        java -jar {input.trimm_path}trimmomatic-0.39.jar PE \
        {input.r1} {input.r2} \
        {output.r1_p} {output.r1_up} \
        {output.r2_p} {output.r2_up} \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads \
        SLIDINGWINDOW:4:20 \
        MINLEN:51
        """

rule star_index:
    input:
        path_assembly = "/uufs/chpc.utah.edu/common/home/kapheim-group1/mrot_maternalRNA/alignment_quantification_mRNAs/genome/Mrot.scaf.fa",
        path_GFF = "/uufs/chpc.utah.edu/common/home/kapheim-group1/mrot_maternalRNA/alignment_quantification_mRNAs/genome/Mrot_v1.1.gff"
    output:
        genome_dir = directory("/uufs/chpc.utah.edu/common/home/kapheim-group1/mrot_maternalRNA/alignment_quantification_mRNAs/genome_index")
    threads:
        8
    shell:"""
        module load star/2.7.8a
        STAR --runThreadN {threads} \
        --runMode genomeGenerate \
        --genomeDir {output.genome_dir} \
        --genomeFastaFiles {input.path_assembly} \
        --sjdbGTFfile {input.path_GFF} \
        --sjdbOverhang 99 \
        --sjdbGTFtagExonParentTranscript Parent \
        --sjdbGTFfeatureExon CDS \
        --genomeChrBinNbits 15  \
        --genomeSAindexNbases 13
        """
# For some reason, this rule does not run well in snakemake. 
# Seems to be deleting the output directory upon starting and then not being able to find it.
# Since this is only run once, I decided just to enter the code above manually through interactive job.
# It will probably be skipped when running in the future, because the output already exists.
# Probably is likely related to output being a directory?

rule star_pe_multi:
    input:
        genome_dir = "/uufs/chpc.utah.edu/common/home/kapheim-group1/mrot_maternalRNA/alignment_quantification_mRNAs/genome_index",
        seqs_R1 = "/uufs/chpc.utah.edu/common/home/kapheim-group1/mrot_maternalRNA/alignment_quantification_mRNAs/trimmed/{sample}_R1.fastq.gz",
        seqs_R2 = "/uufs/chpc.utah.edu/common/home/kapheim-group1/mrot_maternalRNA/alignment_quantification_mRNAs/trimmed/{sample}_R2.fastq.gz"
    params:
        prefix = "/uufs/chpc.utah.edu/common/home/kapheim-group1/mrot_maternalRNA/alignment_quantification_mRNAs/star_aligned/{sample}_",
    output:
        sam = "/uufs/chpc.utah.edu/common/home/kapheim-group1/mrot_maternalRNA/alignment_quantification_mRNAs/star_aligned/{sample}_Aligned.out.sam"
    threads:
        8
    shell:"""
        module load star/2.7.8a
        STAR --runThreadN {threads} \
        --genomeDir {input.genome_dir} \
        --readFilesIn {input.seqs_R1} {input.seqs_R2} \
        --readFilesCommand zcat \
        --outFileNamePrefix {params.prefix} \
        --clip3pNbases 0 \
        --clip5pNbases 0
        """

# From STAR manual: "output unsorted Aligned.out.bam file. The paired ends of an alignment are always adjacent, 
# and multiple alignments of a read are adjacent as well. This ”unsorted” file can be directly 
# used with downstream software such as HTseq, without the need of name sorting. The order 
# of the reads will match that of the input FASTQ(A) files only if one thread is used 
# --runThread 1, and --outFilterType --BySJout is not used"

rule counting:
    input: 
        gff = "/uufs/chpc.utah.edu/common/home/kapheim-group1/mrot_maternalRNA/alignment_quantification_mRNAs/genome/Mrot_v1.1.gff",
        sam = "/uufs/chpc.utah.edu/common/home/kapheim-group1/mrot_maternalRNA/alignment_quantification_mRNAs/star_aligned/{sample}_Aligned.out.sam"
    output:
       "/uufs/chpc.utah.edu/common/home/kapheim-group1/mrot_maternalRNA/alignment_quantification_mRNAs/htseq/counts_{sample}.txt"
    threads:
        8
    shell:"""
        module load python/3.6.3
        python -m HTSeq.scripts.count \
        -m union \
        -s reverse \
        -t CDS \
        -i Parent \
        {input.sam} {input.gff} > {output}
        """