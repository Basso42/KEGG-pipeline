#download BLAST from NCBI
mkdir blast
cd $PWD/blast
wget -nc https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.15.0+-x64-linux.tar.gz       

#will create a subdirectory containing the binaries
tar zxvpf ncbi-blast-2.15.0+-x64-linux.tar.gz

#access program binaries
export PATH=$PATH:$PWD/ncbi-blast-2.15.0+/bin 
#export PATH=$PATH:$/home/onyxia/work/KEGG_Pipeline/ncbi-blast-2.15.0+/bin to write in /home/onyxia/.profile



#a small test 
blastn -help

#Examples 



#wget ftp://ftp.ncbi.nih.gov/refseq/B_taurus/mRNA_Prot/cow.1.protein.faa.gz
#wget ftp://ftp.ncbi.nih.gov/refseq/H_sapiens/mRNA_Prot/human.1.protein.faa.gz
#gunzip *gz

#head cow.1.protein.faa
#grep -c '^>' cow.1.protein.faa #number of sequences in your faasta file

#head -6 cow.1.protein.faa > cow.small.faa

#makeblastdb -in human.1.protein.faa -dbtype prot #creates protein database from human.1.faa file 
## to see how nicolas works with BLAST : https://gitlab.inria.fr/search?search=blast&nav_source=navbar&project_id=26553&search_code=true&repository_ref=EnzBert