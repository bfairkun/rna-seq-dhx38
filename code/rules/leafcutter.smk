rule gtf2IntronsBed:
    """
    this rule uses a  useful script in the leafcutter repo parses out all
    intron regions from a gtf. However, it is not correct bed format in that it
    seems to be not half-closed interval format. In other words, the third
    column is supposed to be non-inclusive so I had to manually subtract 1 from
    the third column output of the script from the third column.
    """
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
        zcat Misc/AnnotatedIntronBeds/Annotated_all_introns.bed.gz | awk -F'\\t' -v OFS='\\t' '{{print $1,$2,$3-1,$4,$5,$6}}' > {output.bed}
        cat {output.bed} | gzip - > {output.bedgz}
        """

rule AnnotateSplicingTypesAndCreateJuncFiles:
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

rule leacutter_cluster:
    input:
        expand("Alignments/JunctionBeds/{sample}.junc", sample=SampleList)
    output:
        juncfilelist = "leafcutter/junctionfiles.txt",
        counts = "leafcutter/leafcutter_perind.counts.gz",
        numers = "leafcutter/leafcutter_perind_numers.counts.gz"
    log:
        "logs/leafcutter_cluster.log"
    params:
        leafcutter_path=config["Path_to_leafcutter_repo"]
    shell:
        """
        echo {input} | tr " " '\\n' > {output.juncfilelist} 2> {log}
        {params.leafcutter_path}clustering/leafcutter_cluster.py -s True -j {output.juncfilelist} -r leafcutter/ &>> {log}
        """

rule leafcutter_ds:
    input:
        numers = "leafcutter/leafcutter_perind_numers.counts.gz",
        groupfile = config["leafcutter_groupfile"],
        exonsbedgz = "Misc/AnnotatedIntronBeds/Annotated_all_exons.txt.gz"
    output:
        "leafcutter/differential_splicing/leafcutter_effect_sizes.txt",
        "leafcutter/differential_splicing/leafcutter_cluster_significance.txt"
    threads: 4
    params:
        leafcutter_path=config["Path_to_leafcutter_repo"],
        outputprefix = "-o leafcutter/differential_splicing/leafcutter",
        extra_params = config["leafcutter_ds_optional_params"]
    log:
        "logs/leafcutter_ds.log"
    shell:
        """
        mkdir -p leafcutter/differential_splicing/
        /software/R-3.4.3-el7-x86_64/bin/Rscript {params.leafcutter_path}scripts/leafcutter_ds.R -p {threads} -e {input.exonsbedgz} {params.outputprefix} {params.extra_params} {input.numers} {input.groupfile} &> {log}
        """

rule move_leafcutter_results:
    input:
        effects = "leafcutter/differential_splicing/leafcutter_effect_sizes.txt",
        sig = "leafcutter/differential_splicing/leafcutter_cluster_significance.txt",
    output:
        effects = "../output/differential_splicing_effect_sizes.txt.gz",
        sig = "../output/differential_splicing_significance.txt.gz"
    shell:
        """
        cat {input.effects} | gzip - > {output.effects}
        cat {input.sig} | gzip - > {output.sig}
        """

rule move_leafcutter_count_table:
    input:
        "leafcutter/leafcutter_perind.counts.gz"
    output:
        numers = config["gitinclude_output"] + "leafcutter.perind.counts.numers.gz",
        denoms = config["gitinclude_output"] + "leafcutter.perind.counts.denoms.gz",
    shell:
        """
        zcat {input} | perl -lne 'if ($.==1) {{print}} else {{$_ =~ s/\d+\///g; print}}' | gzip - > {output.denoms}
        zcat {input} | perl -lne 'if ($.==1) {{print}} else {{$_ =~ s/\/\d+//g; print}}' | gzip - > {output.numers}
        """

rule leafcutter_ds_from_list:
    input:
        numers = "leafcutter/leafcutter_perind_numers.counts.gz",
        groupfile = lambda wildcards: Leafcutter_ds_analyses_params.at[wildcards.analysis_name, "GroupsFile"],
        exonsbedgz = "Misc/AnnotatedIntronBeds/Annotated_all_exons.txt.gz"
    output:
        "leafcutter/differential_splicing/Analyses/{analysis_name}_effect_sizes.txt",
        "leafcutter/differential_splicing/Analyses/{analysis_name}_cluster_significance.txt"
    threads: 4
    params:
        leafcutter_path=config["Path_to_leafcutter_repo"],
        outputprefix = lambda wildcards: "-o leafcutter/differential_splicing/Analyses/{analysis_name}".format(analysis_name=wildcards.analysis_name),
        extra_params = lambda wildcards: Leafcutter_ds_analyses_params.at[wildcards.analysis_name, "ExtraLeafcutterParams"]
    log:
        "logs/leafcutter_ds_analyses/{analysis_name}.log"
    shell:
        """
        /software/R-3.4.3-el7-x86_64/bin/Rscript {params.leafcutter_path}scripts/leafcutter_ds.R -p {threads} -e {input.exonsbedgz} {params.outputprefix} {params.extra_params} {input.numers} {input.groupfile} &> {log}
        """

rule move_leafcutter_ds_from_list_results:
    input:
        effects = "leafcutter/differential_splicing/Analyses/{analysis_name}_effect_sizes.txt",
        sig = "leafcutter/differential_splicing/Analyses/{analysis_name}_cluster_significance.txt"
    output:
        effects = "../output/leafcutter_ds_analyses/{analysis_name}_effect_sizes.txt.gz",
        sig = "../output/leafcutter_ds_analyses/{analysis_name}_cluster_significance.txt.gz"
    shell:
        """
        cat {input.effects} | gzip - > {output.effects}
        cat {input.sig} | gzip - > {output.sig}
        """
