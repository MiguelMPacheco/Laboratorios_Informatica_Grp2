# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 20:33:07 2020

@author: Utilizador
"""

from Bio import Entrez
from Bio import Medline
from Bio import SeqIO
    
    
Entrez.email = "pg@uminho.pt"

handle_seq_search_WDR73 = Entrez.esearch(db="nucleotide", term="Homo sapiens[Organism] AND WDR73[Gene Name]", idtype="acc")
record_seq_WDR73 = Entrez.read(handle_seq_search_WDR73)
print("\n number of WDR73 ids associated with chosen term: \n",record_seq_WDR73["Count"])
print("\n list of WDR73 ids associated with chosen term: \n",record_seq_WDR73["IdList"])

handle_allseqs_WDR73 =Entrez.efetch(db="nucleotide",id=','.join(record_seq_WDR73["IdList"]),rettype = "gb",retmode="text")
save_allseqs_WDR73=open("allseqs_WDR73.txt","w")
save_allseqs_WDR73.write(handle_allseqs_WDR73.read())
save_allseqs_WDR73.close()

handle_allseqs_WDR73 =Entrez.efetch(db="nucleotide",id=','.join(record_seq_WDR73["IdList"]),rettype = "gb",retmode="text")
for index, record in enumerate(SeqIO.parse(handle_allseqs_WDR73,"gb")):
    print("index %i, ID = %s, length %i, with %i features and %i annotations"
          % (index, record.id, len(record.seq), len(record.features), len(record.annotations)))

records_allseqs_WDR73=SeqIO.parse("allseqs_WDR73.txt","gb")

for record in records_allseqs_WDR73:
    print("\n description: \n",record.description)
    print("\n record id: \n",record.id)
    print("\n number of features: \n", len(record.features))
    print("\n number of annotations: \n",len(record.annotations))
    for l in record.annotations:
        print("\n *** annotations ***")
        print("\n molecule_type = ",record.annotations["molecule_type"])
        print("\n accessions = ",record.annotations["accessions"])
        print("\n organism = ",record.annotations["organism"])

records_allseqs_WDR73=list(SeqIO.parse("allseqs_WDR73.txt","gb"))

feature_CDS_gene=[]
list_NP=[]
for record in records_allseqs_WDR73:
    if record.id[:2] == "NG":
        for feat in record.features:
            print("\n *** features ***")
            print(" \n print record.features: \n",feat)
            if feat.type == "CDS" and ("WDR73" in feat.qualifiers["gene"]):
                feature_CDS_gene.append(feat.qualifiers["translation"])
                print("\n gene: ", feat.qualifiers["gene"])
                print("\n product: ", feat.qualifiers["product"])
                print("\n protein_id: ", feat.qualifiers["protein_id"])
                print("\n db_xref: ", feat.qualifiers["db_xref"])
        print("\n list of feature_CDS: \n", feature_CDS_gene)
                
protein_WDR73=feature_CDS_gene[0][0]
print(protein_WDR73)

protein_WDR73_seq=open("protein_WDR73_seq.txt","w")
protein_WDR73_seq.write(protein_WDR73)
protein_WDR73_seq.close()