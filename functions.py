from Bio import SeqIO
from Bio.KEGG import REST

## This notebook aims at building a dataframe composed of eukaryotists annotated genes, including the K0, the DNA sequence and the organisms name.

from bs4 import BeautifulSoup
import requests
import re
import pandas as pd
from itertools import zip_longest
import s3fs
import os
import glob



def get_prokaryotic_organisms():
    # Define the KEGG API URL for all organisms
    kegg_api_url = "http://rest.kegg.jp/list/organism"
    
    try:
        # Send a GET request to the KEGG API
        response = requests.get(kegg_api_url)
        
        # Check if the request was successful (HTTP status code 200)
        if response.status_code == 200:
            # Split the response text into lines
            lines = response.text.strip().split('\n')
            # Extract the organism codes and names from the lines
            organisms = [line.split('\t') for line in lines]
            # Filter for prokaryotic organisms (Bacteria and Archaea)
            prokaryotic_organisms = [org for org in organisms if "Bacteria" in org[3] or "Archaea" in org[3]]
            return prokaryotic_organisms
        else:
            from bs4 import BeautifulSoup
            print(f"Error: Failed to retrieve data from KEGG API. Status code: {response.status_code}")
    except Exception as e:
        print(f"Error: {e}")


def get_eukaryotic_organisms():
    # Define the KEGG API URL for all organisms
    kegg_api_url = "http://rest.kegg.jp/list/organism"
    
    try:
        # Send a GET request to the KEGG API
        response = requests.get(kegg_api_url)
        
        # Check if the request was successful (HTTP status code 200)
        if response.status_code == 200:
            # Split the response text into lines
            lines = response.text.strip().split('\n')
            # Extract the organism codes and names from the lines
            organisms = [line.split('\t') for line in lines]
            # Filter for prokaryotic organisms (Bacteria and Archaea)
            euka_organisms = [org for org in organisms if "Eukaryotes" in org[3]]

            return euka_organisms
        else:
            from bs4 import BeautifulSoup
            print(f"Error: Failed to retrieve data from KEGG API. Status code: {response.status_code}")
    except Exception as e:
        print(f"Error: {e}")


def get_genes_for_organism(organism_code):
    # Define the KEGG API URL for genes associated with the organism
    kegg_api_url = f"http://rest.kegg.jp/list/{organism_code}"
        
    # Make a request to the KEGG API to retrieve genes for the specified organism
    response = requests.get(kegg_api_url)  # Use kegg_api_url here, not kegg_api_base_url
    
    # Check if the request was successful (HTTP status code 200)
    if response.status_code == 200:
        # Split the response content into lines, where each line represents a gene entry
        gene_lines = response.text.strip().split("\n")
    
        # Extract gene IDs from each line
        gene_ids = [line.split("\t")[0] for line in gene_lines]

        return gene_ids  # Return the gene IDs

    return []  # Return an empty list if the request was not successful or no genes found


def grouper(iterable, n, fillvalue=None):
    """
    To iterate on non-overlapping chunks of n elements 
    """
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def check_list_length(lst):
    if len(lst)==0 or (len(lst) == 2 and lst[0]=='['):
        return True
    else:
        return False


def fetch_info(gene_ids):
    gene_group_url = '+'.join(gene_ids)
    full_url = f"http://www.kegg.jp/dbget-bin/www_bget?{gene_group_url}"
    gene_page_response = requests.get(full_url)
    
    # KO
    match = re.findall(r'/entry/(K\d{5})', gene_page_response.text)
    pattern = r'>Entry</span>(.*?)>Organism</span>'
    extracts = re.findall(pattern, gene_page_response.text, re.DOTALL)
    indic_k0 = [True if "KO" in extract else False for extract in extracts]
    index_A = 0
    k0_value = []
    for element_B in indic_k0:
        if element_B and index_A < len(match):
            k0_value.append(match[index_A])
            index_A += 1
        else:
            k0_value.append('')    
            
    #Nt_seq        
    Nt_seq = re.findall(r">nt</code><br>\n(.*?)</td></tr>", gene_page_response.text, re.DOTALL)
    Nt_seq = [seq.replace('<br>', '').replace('\n', '') for seq in Nt_seq]

    #Prot
    pattern = r'>Position</span>(.*?)>NT seq</span>'
    extracts = re.findall(pattern, gene_page_response.text, re.DOTALL)
    indic_AA = [True if "AA seq" in extract else False for extract in extracts]

    return k0_value, Nt_seq, indic_AA


def upload_fold_s3_csv(s3_folder_path, local_path):
    """
    Args: -s3_path: folder to upload .csv file
          -local_csv_file: csv you want to upload
    """
    S3_ENDPOINT_URL = "https://" + os.environ["AWS_S3_ENDPOINT"]
    fs = s3fs.S3FileSystem(client_kwargs={'endpoint_url': S3_ENDPOINT_URL})
    
    BUCKET='gamer35/'
    s3_folder_path = s3_folder_path +'/'
    for fname in glob.glob(local_path + '*.csv'):
        
        # Path to the local CSV file that you want to upload
        LOCAL_CSV_FILE = fname  # Replace with the actual local CSV file path
        file_name_S3 = fname.split('/')[-1]
        
    
        # Upload the CSV file to S3
        with open(LOCAL_CSV_FILE, 'rb') as local_file:
            with fs.open(BUCKET + s3_folder_path + file_name_S3, 'wb') as s3_file:
                s3_file.write(local_file.read())
        
        print(f"Uploaded {BUCKET + s3_folder_path + file_name_S3}")
    print("File uploaded")


def upload_s3_csv(s3_folder_path, file_path):
    """
    Args: -s3_path: folder to upload .csv file
          -local_csv_file: csv you want to upload
    """
    S3_ENDPOINT_URL = "https://" + os.environ["AWS_S3_ENDPOINT"]
    fs = s3fs.S3FileSystem(client_kwargs={'endpoint_url': S3_ENDPOINT_URL})
    
    BUCKET='gamer35/'
    s3_folder_path = s3_folder_path +'/'

        
        # Path to the local CSV file that you want to upload
    LOCAL_CSV_FILE = file_path  # Replace with the actual local CSV file path
    file_name_S3 = file_path.split('/')[-1]
        
    
        # Upload the CSV file to S3
    with open(LOCAL_CSV_FILE, 'rb') as local_file:
        with fs.open(BUCKET + s3_folder_path + file_name_S3, 'wb') as s3_file:
            s3_file.write(local_file.read())
        
    print(f"Uploaded {BUCKET + s3_folder_path + file_name_S3}")




def del_s3_csv(s3_folder_path, file_to_delete):
    """
    inputs: - s3_folder_path: the s3_folder path inside the bucket
            - file_to_delete: the name of the .csv file to be deleted

    fc: delete csv file from a S3 folder
    """
    S3_ENDPOINT_URL = "https://" + os.environ["AWS_S3_ENDPOINT"]
    fs = s3fs.S3FileSystem(client_kwargs={'endpoint_url': S3_ENDPOINT_URL})
    BUCKET='gamer35/'
    file_list = fs.ls(f"{BUCKET + s3_folder_path}")
        
    if f'{BUCKET + s3_folder_path}/{file_to_delete}' in file_list:
        fs.rm(f'{BUCKET + s3_folder_path}/{file_to_delete}')
        print(f"Deleted {file_to_delete} from {s3_folder_path}")
    else:
        print(f"{file_to_delete} not found in {s3_folder_path}")

def check_if_prot_list(page_response_text):

    pattern = r'>Position</span>(.*?)>NT seq</span>'
    extracts = re.findall(pattern, page_response_text, re.DOTALL)
    indic_AA = [True if "AA seq" in extract else False for extract in extracts]

    return indic_AA


def import_s3_folder(s3_folder, local_path, ext = 'csv'):
    """
    Imports all the .{format} file to {local_path} from a given S3 folder
    """
    S3_ENDPOINT_URL = "https://" + os.environ["AWS_S3_ENDPOINT"]
    fs = s3fs.S3FileSystem(client_kwargs={'endpoint_url': S3_ENDPOINT_URL})
    BUCKET='gamer35/'

    file_list = fs.ls(BUCKET + s3_folder)
    
    # Loop through the CSV files and download them
    for file_path in file_list:
        if file_path.endswith(f'.{ext}'):
            # Construct the local file path where the CSV will be saved
            local_file_path = os.path.join(local_path, os.path.basename(file_path))
            
            # Download the CSV file from S3 and save it locally
            with fs.open(file_path, 'rb') as s3_file, open(local_file_path, 'wb') as local_file:
                local_file.write(s3_file.read())
    
    print(f"Downloaded CSV files to {local_path}")

