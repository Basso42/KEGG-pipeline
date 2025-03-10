# The goal of this script is to count the number of new organisms have been added to KEGG database since we have scrapped the ones we have in our dataset. (which have their name contained in organisms_{DATE}.txt)

from functions import get_prokaryotic_organisms
from random import randint

prok = get_prokaryotic_organisms()
prok = [organism[1] for organism in prok]

name_file = 'organisms_15th_nov'
kegg_organism_url = "https://www.genome.jp/kegg-bin/show_organism?org="
genes_organism_url = 'http://rest.kegg.jp/list/'


# reading file
my_file = open(f"/home/onyxia/work/KEGG_Pipeline/{name_file}.txt", "r") 
data = my_file.read()
org_yet_scrapped = data.split("\n") 

# fetching unscrapepd genomes
unscrapped_org = [org for org in prok if org not in org_yet_scrapped]
rand_org = unscrapped_org[randint(0,len(unscrapped_org))]


# File name
output_file = 'bonjour.txt'

# Writing to the file
with open(output_file, 'w') as f:
    # Writing information about the number of unscrapped organisms
    unscrapped_count = len(prok) - len(org_yet_scrapped)
    f.write(f"Relative to {name_file}, the KEGG dataset currently has {unscrapped_count} prokaryotic organisms that haven't been scraped yet.\n \n")

    # Writing information about a random organism to look at and its list of genes
    random_organism_url = kegg_organism_url + rand_org
    genes_organism_url = genes_organism_url + rand_org
    f.write(f"Random organism to look at: {random_organism_url} and its list of genes: {genes_organism_url}\n \n")

    # Writing the list of unscrapped organisms
    if unscrapped_org:
        f.write("List of unscrapped organisms:\n")
        for org in unscrapped_org:
            f.write(f"{org}\n")
    else:
        f.write("All organisms have been successfully scraped!\n")

# File has been automatically closed after exiting the 'with' block



print(f"Relatively to {name_file}, the KEGG dataset has now {len(prok) - len(org_yet_scrapped)} prokaryotes organisms that we haven't scrapped yet. ")

