"""
This script allows to run NCBI BLAST+ in local to perform annotations of a test set, given a train set by extracting the information given in the alignemnet of the sequences.
"""
import os
import pandas as pd
from pathlib import Path
import shutil
from Bio import SearchIO
from utils import generate_fasta_file
import os
import polars as pl


class SequenceKNN:
    def __init__(
        self,
        path_train_parquet,
        path_test_parquet,
        nb_thread,
        path_output_pred,
        tmp_folder=Path("data/tmp_blastp/"),
        tool="BLASTp",
    ):
        """
        Creates, compares and annotates protein sequences having a certain enzyme

        Args:
            - path_train_parquet: path of the training dataset (parquet)
            - path_test_parquet: path test dataset (parquet)
            - nb_thread: number of processes used in the blast annot
            - path_output_pred: path output predictiom
            - tmp_folder: path of the temporary folder
            - tool: command passed to align the sequences with BLAST+
        """
        
        self.tool = tool
        self.df_train = pd.read_parquet(path_train_parquet)
        self.df_test = pd.read_parquet(path_test_parquet)
        self.tmp_folder = tmp_folder
        self.nb_thread = nb_thread
        if not os.path.exists(tmp_folder):
            os.mkdir(tmp_folder)
        self.name_fasta_train = "train.fasta"
        self.name_fasta_test = "test.fasta"
        self.path_fasta_train = tmp_folder / Path(self.name_fasta_train)
        self.path_fasta_test = tmp_folder / Path(self.name_fasta_test)
        self.output_query = "res_blastp.txt"
        self.path_output_pred = path_output_pred


    def launch_pipeline(self):
        """
        Generate fasta files for train and test
        Args:
            - df_train: pandas dataframe of the train sample
            - path_fasta_train: path to generate fasta file
        """

        generate_fasta_file(df=self.df_train, fasta_path=self.path_fasta_train) #create fasta from train
        self.create_blast_db() #create database
        generate_fasta_file(df=self.df_test, fasta_path=self.name_fasta_test) #create fasta for test
        self.query_seq() #do the alignments
        self.parse_res_and_get_pred()
        shutil.rmtree(self.tmp_folder)

    def create_blast_db(self): 
        """
        Create blast database
        """
        os.chdir(self.tmp_folder)
        if self.tool == "BLASTp":
            command = f"makeblastdb -in {str(self.name_fasta_train)} -dbtype prot"
        elif self.tool == "DIAMOND":
            command = f"diamond makedb  --in {str(self.name_fasta_train)} --dbtype prot"
        else:
            raise ValueError("tool unkwown")
        print("command:", command)
        os.system(command)

    def query_seq(self):
        """
        Perform alignment to compare train and test sequences
        """
        if not os.path.exists(self.output_query):
            if self.tool == "BLASTp":
                command = (
                    "blastp -query "
                    + self.name_fasta_test
                    + " -db "
                    + self.name_fasta_train
                    + " -out "
                    + self.output_query
                    + " -outfmt 6"  # XML output or tabular
                    + " -num_threads "
                    + str(self.nb_thread)
                    # + " -mt_mode 1" # Deprecated options? not working anymore
                )
            elif self.tool == "DIAMOND":
                command = (
                    "diamond blastp -q "
                    + self.name_fasta_test
                    + " -d "
                    + self.name_fasta_train
                    + " -o "
                    + self.output_query
                    + " --outfmt 5"  # XML output
                    + " --threads "
                    + str(self.nb_thread)
                    # + " -mt_mode 1" # Deprecated options? not working anymore
                )
            else:
                raise ValueError("tool unkwown")
            print("command:", command)
            os.system(command)
        else:
            print("File already computed")

    def create_dico_res_blastp(self):
        """
        Parses results of alignment
        """
        qresults = SearchIO.parse(path_data+"tmp_blastp/res_blastp.txt", "blast-tab")
        print(qresults) #new
        dico = {}
        dico_seq_identity = {}
        for qresult in qresults:
            if len(qresult.hits) != 0:
                #print(dir(qresult.hits[0][0]))
                #print("qresult.hits:", qresult.hits[0][0])
                #print("qresult.hits:", qresult.hits[0][0].ident_pct) #ident_num this is the nuñber of same residues
                #print("qresult.hits:", qresult.hits[0][0].aln_span)
                #print("qresult.hits:", qresult.hits[0][0].pos_num)
                #print("qresult.hits:", len(qresult.hits[0][0].query))
                #print("qresult.hits:", len(qresult.hits[0][0].hit))

                dico[qresult.id] = [res.id for res in qresult.hits]
                
                dico_seq_identity[qresult.id] = [qresult.hits[i][0].ident_pct for i in range(len(qresult.hits)) #map le nom d'une requête à une liste de partage d'identité entre 0 et !
                    #res[0].ident_num / qresult.seq_len for res in qresult.hits
                ] # number of identical residues/sequence length of the hit
        return dico, dico_seq_identity

    def create_dico_ec_train(self):
        return {
            row["gene_id"]: row["K0"]
            for index, row in self.df_train.iterrows()
        }

    def parse_res_and_get_pred(self):
        """
        Parse output of blastp file and generates prediction from it
        """
        dico_query_to_target, dico_seq_identity = self.create_dico_res_blastp() #fetch dicts mapping from query to target and % shared with target
        dico_train_id_to_ec = self.create_dico_ec_train() #dict mapping training gene to its function
        os.chdir("../../")
        print("self.path_output_pred:", self.path_output_pred)
        with open(self.path_output_pred, "w") as output_pred:
            output_pred.write("AA_seq,pred_K0,seq_identity,train_id\n") #writes file with 
            for _, row in self.df_test.iterrows():
                sequence = row["AA_seq"]
                kegg_id = row["gene_id"]

#if a query was hit by a target gene, than note the KO predicted and % of sequence identity shared (should not be higher than 40%, MMSeqs2)
                if kegg_id in dico_query_to_target.keys(): 
                    list_ids = dico_query_to_target[kegg_id]
                    list_seq_seq_identity = dico_seq_identity[kegg_id]
                    #if list_ids[0][:3]=="amb": 
                    #    print(list_ids[0])
                    #   list_ids[0]='mag:'+list_ids[0]
                    #if list_ids[0][:3]=="Aci":
                    #    print(list_ids[0])
                    #    list_ids[0]='tsa:'+list_ids[0]
                    pred = dico_train_id_to_ec[list_ids[0]]
                    seq_identity = list_seq_seq_identity[0]
                    

                #if no hit, there is no prediction
                else: 
                    pred = "No"
                    seq_identity = 0.0
                output_pred.write(
                    sequence + "," + pred + "," + str(seq_identity) + "," + str(list_ids[0]) + "\n" 
                )


if __name__=="__main__":
    
    #from utils import *
    #from blast_pred import *
    path_data="/home/onyxia/work/KEGG_Pipeline/data/"

    os.system("wget -nc https://minio.lab.sspcloud.fr/gamer35/KEGG_db/train-validation-test/train_dataset.parquet -P data")
    os.system("wget -nc https://minio.lab.sspcloud.fr/gamer35/KEGG_db/train-validation-test/test_dataset.parquet -P data")
    os.system(f"mkdir -P {path_data}tmp_blastp/")

    train = pl.read_parquet(path_data + "train_dataset.parquet")
    test = pl.read_parquet(path_data + "test_dataset.parquet")

    train = train.with_columns(pl.col('gene_id').str.replace(':','_')).write_parquet(path_data + "train_dataset.parquet")
    test = test.with_columns(pl.col('gene_id').str.replace(':','_')).write_parquet(path_data + "test_dataset.parquet")


    test=SequenceKNN(path_train_parquet= path_data+"train_dataset.parquet",
        path_test_parquet= path_data+"test_dataset.parquet",
        nb_thread=100,
        path_output_pred=path_data + "pred",
        tmp_folder=Path(path_data+"tmp_blastp/"),
        tool="BLASTp")

    test.launch_pipeline()