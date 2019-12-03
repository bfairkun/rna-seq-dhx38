rule STAR_alignment_FirstPass_PE:
    input:
        index_human = config["Human_ref"]["genome_dir"] + "chrLength.txt",
        R1 = "Sample_Fastq_PE/{sample}.R1.fastq.gz",
        R2 = "Sample_Fastq_PE/{sample}.R2.fastq.gz"
    log:
        "logs/STAR_FirstPass/{sample}.log"
    threads: 8
    params:
        index = Human_genomeDir
    output:
        SJout = "Alignments/FirstPass_PE/{sample}/SJ.out.tab"
    shell:
        """
        ulimit -v 31538428657
        STAR --runThreadN {threads} --genomeDir {params.index} --readFilesIn {input.R1} {input.R2} --readFilesCommand zcat --outSAMmultNmax 1 --outSAMtype None --outFileNamePrefix Alignments/FirstPass_PE/{wildcards.sample}/ &> {log}
        """

rule STAR_alignment_FirstPass:
    input:
        index_human = config["Human_ref"]["genome_dir"] + "chrLength.txt",
        R1 = "Sample_Fastq/{sample}.fastq.gz"
    log:
        "logs/STAR_FirstPass/{sample}.log"
    threads: 8
    params:
        index = Human_genomeDir
    output:
        SJout = "Alignments/FirstPass/{sample}/SJ.out.tab"
    shell:
        """
        ulimit -v 31538428657
        STAR --runThreadN {threads} --genomeDir {params.index} --readFilesIn {input.R1} --readFilesCommand zcat --outSAMmultNmax 1 --outSAMtype None --outFileNamePrefix Alignments/FirstPass/{wildcards.sample}/ &> {log}
        """

rule STAR_make_index_human_AfterFirstPass:
    input:
        fasta=config["Human_ref"]["genome_fasta"],
        gtf=config["Human_ref"]["genome_gtf"],
        FirstPass_sjdb=expand( "Alignments/FirstPass_PE/{sample}/SJ.out.tab", sample=SampleList),
    output:
        index = "Alignments/FirstPass_sjdb_index/chrLength.txt"
    log:
        "logs/STAR_make_index_human_AfterFirstPass.log"
    shell:
        """
        STAR --runMode genomeGenerate --genomeSAsparseD 2 --runThreadN 4 --genomeDir Alignments/FirstPass_sjdb_index/ --sjdbGTFfile {input.gtf} --genomeFastaFiles {input.fasta} --sjdbFileChrStartEnd {input.FirstPass_sjdb} &> {log}
        """

rule STAR_alignment_SecondPass_PE:
    input:
        index = "Alignments/FirstPass_sjdb_index/chrLength.txt",
        R1 = "Sample_Fastq_PE/{sample}.R1.fastq.gz",
        R2 = "Sample_Fastq_PE/{sample}.R2.fastq.gz"
    log:
        "logs/STAR_SecondPass_PE/{sample}.log"
    threads: 8
    output:
        SJout = "Alignments/SecondPass/{sample}/SJ.out.tab",
        bam = "Alignments/SecondPass/{sample}/Aligned.sortedByCoord.out.bam",
        bai = "Alignments/SecondPass/{sample}/Aligned.sortedByCoord.out.bam.bai",
        # R1Unmapped = "Alignments/SecondPass_PE/{sample}/Unmapped.out.mate1.fastq.gz",
        # R2Unmapped = "Alignments/SecondPass_PE/{sample}/Unmapped.out.mate2.fastq.gz"
    shell:
        """
        ulimit -v 31538428657
        STAR --runThreadN {threads} --genomeDir Alignments/FirstPass_sjdb_index --readFilesIn {input.R1} {input.R2} --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix Alignments/SecondPass/{wildcards.sample}/ --outReadsUnmapped Fastx --alignEndsType EndToEnd &> {log}
        samtools index {output.bam}
        """

rule sortBam:
    input:
        "{Filepath}.bam"
    output:
        sortedbam = "{Filepath}.sorted.bam",
        index = "{Filepath}.sorted.bam.bai"
    shell:
        """
        samtools sort {input} > {output.sortedbam}
        samtools index {output.sortedbam}
        """

rule STAR_alignment_SecondPass:
    input:
        index = "Alignments/FirstPass_sjdb_index/chrLength.txt",
        R1 = "Sample_Fastq/{sample}.fastq.gz"
    log:
        "logs/STAR_SecondPass/{sample}.log"
    threads: 8
    output:
        SJout = "Alignments/SecondPass_SE/{sample}/SJ.out.tab",
        bam = "Alignments/SecondPass_SE/{sample}/Aligned.sortedByCoord.out.bam",
        bai =  "Alignments/SecondPass_SE/{sample}/Aligned.sortedByCoord.out.bam.bai"
    shell:
        """
        ulimit -v 31538428657
        STAR --runThreadN {threads} --genomeDir Alignments/FirstPass_sjdb_index --readFilesIn <(zcat {input.R1} | fastx_trimmer -l 46) --outSAMtype BAM SortedByCoordinate --outWigType wiggle --outFileNamePrefix Alignments/SecondPass/{wildcards.sample}/ &> {log}
        samtools index {output.bam}
        """
