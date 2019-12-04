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

rule MergedIntronList:
    input:
        SJout_AS_annotated = expand( "Alignments/SecondPass/{sample}/SJ.out.annotated.tab", sample=SampleList),
        fai =config["Human_ref"]["genome_fasta"] + ".fai"
    output:
        "Misc/FullIntronList.bed"
    shell:
        """
        cat {input.SJout_AS_annotated} | awk -F'\\t' -v OFS='\\t' '{{ print $1,$2,$3,$1"_"$2"_"$3"_"$6, "." ,$6,$7 }}' | sort | uniq | bedtools sort -i - -faidx {input.fai} > {output}
        """

rule GetMergedIntronDonorAndAcceptorSequences:
    """
    Will be used later to score donor, acceptor, and branch motifs
    """
    input:
        FullIntronList = "Misc/FullIntronList.bed",
        fasta=config["Human_ref"]["genome_fasta"],
    output:
        DonorFasta = "Misc/FullIntronList.Donors.fasta",
        BranchRegionFasta = "Misc/FullIntronList.BranchRegion.fasta",
        AcceptorFasta = "Misc/FullIntronList.Acceptors.fasta"
    shell:
        """
        awk -F'\\t' -v OFS='\\t' '$6=="+" {{print $1,$2-3,$2+6,$4,$5,$6}} $6=="-" {{print $1, $3-6, $3+3, $4, $5, $6}}' {input.FullIntronList} | bedtools getfasta -fi {input.fasta} -s -name -bed - > {output.DonorFasta}
        awk -F'\\t' -v OFS='\\t' '$6=="+" {{print $1,$3-100,$3,$4,$5,$6}} $6=="-" {{print $1, $2, $2+100, $4, $5, $6}}' {input.FullIntronList} | bedtools getfasta -fi {input.fasta} -s -name -bed - > {output.BranchRegionFasta}
        awk -F'\\t' -v OFS='\\t' '$6=="+" {{print $1,$3-20,$3+3,$4,$5,$6}} $6=="-" {{print $1, $2-3, $2+20, $4, $5, $6}}' {input.FullIntronList} | bedtools getfasta -fi {input.fasta} -s -name -bed - > {output.AcceptorFasta}
        """

rule ScoreBranchpoints:
    input:
        BranchRegionFasta = "Misc/FullIntronList.BranchRegion.fasta",
    output:
        Branchpoints = "Misc/MotifScores/Branchpoints.txt"
    params:
        bp_svm_script = config["Path_to_svm_bpfinder"]
    shell:
        """
        {params.bp_svm_script} -i {input.BranchRegionFasta} -s Hsap > {output}
        """

rule CalculateBestBpPerIntron:
    input:
        Branchpoints = "Misc/MotifScores/Branchpoints.txt"
    output:
        Branchpoints = "Misc/MotifScores/Branchpoints.bestPerIntron.txt"
    shell:
        """
        perl scripts/calculate_best_BP_per_intron.pl < {input.Branchpoints} > {output.Branchpoints}
        """

rule GetAnnotated_Donor_And_Acceptors:
    input:
        bed = "Misc/AnnotatedIntronBeds/Annotated_all_introns.bed",
        fasta=config["Human_ref"]["genome_fasta"],
        fai =config["Human_ref"]["genome_fasta"] + ".fai"
    output:
        AnnotatedDonorFa = "Misc/MotifScores/PWMs/AnnotatedDonors.fa",
        AnnotatedAccepterFa = "Misc/MotifScores/PWMs/AnnotatedAcceptors.fa"
    shell:
        """
        awk -F'\\t' -v OFS='\\t' '$6=="+" {{print $1,$2-3,$2+6,".",".",$6}} $6=="-" {{print $1, $3-6, $3+3, $1"_"$2"_"$3"_"$6, ".", $6}}' {input.bed} | sort | uniq | bedtools sort -i - -faidx {input.fai} | bedtools getfasta -fi {input.fasta} -s -name+ -bed - > {output.AnnotatedDonorFa}
        awk -F'\\t' -v OFS='\\t' '$6=="+" {{print $1,$3-20,$3+3,".",".",$6}} $6=="-" {{print $1, $2-3, $2+20, ".", ".", $6}}' {input.bed} | sort | uniq | bedtools sort -i - -faidx {input.fai} | bedtools getfasta -fi {input.fasta} -s -name+ -bed - > {output.AnnotatedAccepterFa}
        """

rule ScoreDonorAndAcceptorPSSM:
    input:
        AnnotatedDonorFa = "Misc/MotifScores/PWMs/AnnotatedDonors.fa",
        AnnotatedAccepterFa = "Misc/MotifScores/PWMs/AnnotatedAcceptors.fa",
        DonorFasta = "Misc/FullIntronList.Donors.fasta",
        AcceptorFasta = "Misc/FullIntronList.Acceptors.fasta"
    output:
        DonorsScored = "Misc/MotifScores/Donors.txt",
        AcceptorsScored = "Misc/MotifScores/Acceptors.txt"
    shell:
        """
        python3 scripts/CreateDonorAndAcceptorPWM.py {input.AnnotatedDonorFa} {input.AnnotatedAccepterFa} {input.DonorFasta} {input.AcceptorFasta} {output.DonorsScored} {output.AcceptorsScored}
        """

rule MergeIntronFeatureScores:
    input:
        DonorsScored = "Misc/MotifScores/Donors.txt",
        AcceptorsScored = "Misc/MotifScores/Acceptors.txt",
        Branchpoints = "Misc/MotifScores/Branchpoints.bestPerIntron.txt",
        IntronList = "Misc/FullIntronList.bed"
    output:
        "../output/IntronFeatures.txt.gz"
    shell:
        """
        Rscript scripts/MergeIntronFeatures.R
        gzip ../output/IntronFeatures.txt 
        """

