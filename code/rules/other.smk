rule gtf2IntronsBed:
    input:
        gtf=config["Human_ref"]["genome_gtf"],
    output:
        bedgz = "Misc/AnnotatedIntronBeds/Annotated_all_introns.bed.gz",
        bed = "Misc/AnnotatedIntronBeds/Annotated_all_introns.bed",
        exonsbedgz = "Misc/AnnotatedIntronBeds/Annotated_all_exons.txt.gz"
    params:
        leafcutter_path=config["Path_to_leafcutter_repo"]
    shell:
        """
        {params.leafcutter_path}leafviz/gtf2leafcutter.pl -o Misc/AnnotatedIntronBeds/Annotated {input.gtf}
        zcat {output.bedgz} > {output.bed}
        """

rule AnnotateSplicingTypes:
    input:
        SJout = "Alignments/SecondPass/{sample}/SJ.out.tab",
        AnnotatedIntrons = "Misc/AnnotatedIntronBeds/Annotated_all_introns.bed"
    output:
        SJOutASBed = "Alignments/JunctionBeds/{sample}.junc",
        SJout_AS_annotated = "Alignments/SecondPass/{sample}/SJ.out.annotated.tab"
    shell:
        """
        awk -F'\\t' -v OFS='\\t' '$4=="2" {{print $1,$2-1,$3,".",$7,"-"}} $4=="1" {{print $1,$2-1,$3,".",$7,"+"}}' {input.SJout} | tee {output.SJOutASBed} | scripts/AnnotateSplicingType.py -I - -A {input.AnnotatedIntrons} -O {output.SJout_AS_annotated}
        """

rule faidx:
    input: config["Human_ref"]["genome_fasta"]
    output: config["Human_ref"]["genome_fasta"] + ".fai"
    shell: "samtools faidx {input}"

rule BedgraphAndBigwigs:
    input:
        bam = "Alignments/SecondPass/{sample}/Aligned.sortedByCoord.out.bam",
        fai =config["Human_ref"]["genome_fasta"] + ".fai"
    output:
        PlusBg = "Alignments/SecondPass/{sample}/Plus.bg",
        MinusBg = "Alignments/SecondPass/{sample}/Minus.bg",
        UnstrandBg = "Alignments/SecondPass/{sample}/Unstranded.bg",
        PlusBw = "Alignments/SecondPass/{sample}/Plus.bw",
        MinusBw = "Alignments/SecondPass/{sample}/Minus.bw",
        UnstrandBw = "Alignments/SecondPass/{sample}/Unstranded.bw",
    shell:
        """
        samtools view -bh -F256 {input.bam} | bedtools genomecov -max 10000 -ibam - -bg -strand + -split -scale $(bc <<< "scale=3;1000000000/$(samtools idxstats {input.bam} | awk '{{sum+=$2}} END{{printf sum}}')") | sort -k1,1 -k2,2n > {output.PlusBg}
        samtools view -bh -F256 {input.bam} | bedtools genomecov -max 10000 -ibam - -bg -strand - -split -scale $(bc <<< "scale=3;1000000000/$(samtools idxstats {input.bam} | awk '{{sum+=$2}} END{{printf sum}}')") | sort -k1,1 -k2,2n > {output.MinusBg}
        samtools view -bh -F256 {input.bam} | bedtools genomecov -max 10000 -ibam - -bg -split -scale $(bc <<< "scale=3;1000000000/$(samtools idxstats {input.bam} | awk '{{sum+=$2}} END{{printf sum}}')") | sort -k1,1 -k2,2n > {output.UnstrandBg}
        bedGraphToBigWig {output.PlusBg} {input.fai} {output.PlusBw}
        bedGraphToBigWig {output.MinusBg} {input.fai} {output.MinusBw}
        bedGraphToBigWig {output.UnstrandBg} {input.fai} {output.UnstrandBw}
        """

rule copyBigwigs:
    input:
        UnstrandBw = "Alignments/SecondPass/{sample}/Unstranded.bw",
    output:
        UnstrandedBigwigs="UnstrandedBigwigs/{sample}.bw"
    shell:
        "cp {input} {output}"

rule MakeTopHatJuncBedFiles:
    input:
        SJ = "Alignments/SecondPass/{sample}/SJ.out.tab",
    output:
        bed = "Junctions/{sample}.junctions.bed.gz",
    log:
        "logs/MakeTopHatJuncBedFiles/{sample}.log"
    shell:
        """
        bash scripts/SJ_to_junctions.sh {input.SJ} Junctions/{wildcards.sample}.junctions.bed
        gzip Junctions/{wildcards.sample}.junctions.bed
        """

rule CountsPerBam:
    input:
        expand("Alignments/SecondPass/{sample}/Aligned.sortedByCoord.out.bam", sample=SampleList),
    output:
        config["gitinclude_output"] + "CountsPerBam.txt"
    log:
        "logs/CountsPerBam.log"
    shell:
        """
        bash scripts/CountReadsPerBam.sh {output} {input} 2> {log}
        """
