rule SubreadCount:
    input:
        expand("Alignments/SecondPass/{sample}/Aligned.sortedByCoord.out.bam", sample=SampleList),
    output:
        Output = "GeneExpression/GeneExpressionCountTable.subread.txt.gz",
        OutputGit = config["gitinclude_output"] + "GeneExpressionCountTable.subread.txt.gz",
        OutputSummary = "GeneExpression/GeneExpressionCountTable.subread.txt.summary",
    log:
        "logs/SubreadCount.txt"
    params:
        Gtf = config["Human_ref"]["genome_gtf"],
        featureCounts = config["Path_to_subread_feature_counts_executable"]
    shell:
        """
        {params.featureCounts} -a {params.Gtf} -o GeneExpression/GeneExpressionCountTable.temp.txt {input} &> {log}
        # Format the header for the count table to just contain sample name
        # instead of bam filepath
        cat GeneExpression/GeneExpressionCountTable.temp.txt | sed -e '2s/Alignments\/SecondPass\///g' | sed -e '2s/\/Aligned.sortedByCoord.out.bam//g' | gzip - > {output.Output}
        cat GeneExpression/GeneExpressionCountTable.temp.txt.summary |  sed -e '1s/Alignments\/SecondPass\///g' | sed -e '1s/\/Aligned.sortedByCoord.out.bam//g' > {output.OutputSummary}
        rm GeneExpression/GeneExpressionCountTable.temp.txt.summary GeneExpression/GeneExpressionCountTable.temp.txt 
        cp {output.Output} {output.OutputGit}
        """
