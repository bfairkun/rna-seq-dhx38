import pandas as pd
import os

###### Config file and sample sheets #####
configfile: "config.yaml"

# samples = pd.read_csv(config["samples"],sep='\t', index_col=0)
FastQList = pd.read_csv(config["FastQList"],sep='\t')

SRA_fastq_for_download = expand("SRA_Fastq/{SRA_Accension}.fastq.gz", SRA_Accension=[x for x in FastQList["SRA"].tolist() if str(x) != 'nan'])

# print( expand("SRA_Fastq/{SRA_Accension}.fastq.gz", SRA_Accension=[x for x in FastQList.loc[FastQList['Sample'] == "19201_cheRNA_1"]["R1"].tolist() if str(x) != 'nan'] ))

SampleList = set(FastQList['Sample'])
SampleList_PE = set(FastQList[FastQList.R2.notnull()]['Sample'])

# # How to access values in samples.tsv
# print(samples)
# print( expand("Hello {sample}", sample=samples.index) )
# print( samples.at["A", "R1"] )

Leafcutter_ds_analyses_params = pd.read_csv(config["leafcutter_ds_analysis_params"], sep='\t', index_col=0)
Leafcutter_ds_analyses_params.fillna('', inplace=True)
Leafcutter_ds_analyses_list = Leafcutter_ds_analyses_params.index
