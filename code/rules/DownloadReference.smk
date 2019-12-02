rule Download_hg38:
    output:
        fa = config["Human_ref"]["genome_fasta"],
        gtf = config["Human_ref"]
    shell:
        """
        # get chromosome fasta and concatentate
        wget "ftp://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.*.fa.gz"
        zcat Homo_sapiens.GRCh38.dna.chromosome.*.fa.gz > {output.fa}
        rm Homo_sapiens.GRCh38.dna.chromosome.*.fa.gz

        # get chromosome gtf
        wget -O {output.gtf} "ftp://ftp.ensembl.org/pub/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh38.98.chr.gtf.gz"
        """

Human_genomeDir = config["Human_ref"]["genome_dir"][:-1]

rule STAR_make_index_human:
    input:
        fasta=config["Human_ref"]["genome_fasta"],
        gtf=config["Human_ref"]["genome_gtf"]
    output:
        index = config["Human_ref"]["genome_dir"] + "chrLength.txt",
    log:
        "logs/STAR_make_index_human.log"
    params:
        genomeDir = Human_genomeDir
    shell:
        """
        STAR --runMode genomeGenerate --genomeSAsparseD 2 --runThreadN 4 --genomeDir {params.genomeDir} --sjdbGTFfile {input.gtf} --genomeFastaFiles {input.fasta} &> {log}
        """
