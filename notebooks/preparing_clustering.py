import pandas as pd
import polars as pl
from Bio import SeqIO
import polars.selectors as cs
import os

os.system("mkdir -p /home/onyxia/work/KEGG_Pipeline/data")
os.system("wget -nc https://minio.lab.sspcloud.fr/gamer35/KEGG_db/Prok_proteins.parquet -P /home/onyxia/work/KEGG_Pipeline/data")


data_path= "/home/onyxia/work/KEGG_Pipeline/data/"

df = pd.read_parquet(data_path + "Prok_proteins.parquet")

df_polars = pl.read_parquet(data_path + "Prok_proteins.parquet")
df_polars = df_polars.unique()\
            .with_columns(pl.col('gene_id').str.replace(':','_'))


df_testin = df_polars.with_columns(pl.when(cs.string().str.lengths() >= 700)
                  .then(cs.string().str.slice(0, 700))
                  .otherwise(cs.string())
                  .keep_name()
                )
df_testin = df_testin.to_pandas()

# Function to write a FASTA file from DataFrame
def write_fasta_from_dataframe(dataframe, output_file):
    with open(output_file, 'w') as fasta_file:
        for index, row in dataframe.iterrows():
            fasta_file.write(f'>{row["gene_id"]}\n{row["AA_seq"]}\n')

# Specify the output file name
output_filename = data_path + 'truncated_proteins.fasta'

# Call the function to write the FASTA file
write_fasta_from_dataframe(df_testin, output_filename)

print(f'FASTA file "{output_filename}" has been created.')
#writes a 5go fasta file in 11 minutes