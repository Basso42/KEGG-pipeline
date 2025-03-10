cd /home/onyxia/work/KEGG_Pipeline/data

#conda install -c conda-forge -c bioconda mmseqs2 -y

echo "Preparing the dataset"
python /home/onyxia/work/KEGG_Pipeline/notebooks/preparing_clustering.py

echo "Creating the database"
mmseqs createdb truncated_proteins.fasta DB

sleep 3

echo "Beginning to cluster"

mmseqs easy-cluster truncated_proteins.fasta clusterRes tmp --min-seq-id 0.4 --cov-mode 0
#mmseqs easy-cluster DB DB_clu tmp --min-seq-id 0.4 --cov-mode 0
#mmseqs createtsv DB DB DB_clu DB_clu.tsv