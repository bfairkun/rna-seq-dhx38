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
        groupfile = config["leafcutter_groupfile"]
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
        {params.leafcutter_path}scripts/leafcutter_ds.R -p {threads} {params.outputprefix} {params.extra_params} {input.numers} {input.groupfile} &> {log}
        """
