# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 20:38:09 2020

@author: Utilizador
"""

from Bio import Entrez
from Bio import Medline
from Bio import SeqIO
    
    
Entrez.email = "pg@uminho.pt"

handle_seq_search_PARD6B = Entrez.esearch(db="nucleotide", term="Homo sapiens[Organism] AND PARD6B[Gene Name]", idtype="acc")
record_seq_PARD6B = Entrez.read(handle_seq_search_PARD6B)
print("\n number of PARD6B ids associated with chosen term: \n",record_seq_PARD6B["Count"])
print("\n list of PARD6B ids associated with chosen term: \n",record_seq_PARD6B["IdList"])

handle_allseqs_PARD6B =Entrez.efetch(db="nucleotide",id=','.join(record_seq_PARD6B["IdList"]),rettype = "gb",retmode="text")
save_allseqs_PARD6B=open("allseqs_PARD6B.txt","w")
save_allseqs_PARD6B.write(handle_allseqs_PARD6B.read())
save_allseqs_PARD6B.close()

handle_allseqs_PARD6B =Entrez.efetch(db="nucleotide",id=','.join(record_seq_PARD6B["IdList"]),rettype = "gb",retmode="text")
for index, record in enumerate(SeqIO.parse(handle_allseqs_PARD6B,"gb")):
    print("index %i, ID = %s, length %i, with %i features and %i annotations"
          % (index, record.id, len(record.seq), len(record.features), len(record.annotations)))

records_allseqs_PARD6B=SeqIO.parse("allseqs_PARD6B.txt","gb")



#falta dict={id:seq}

for record in records_allseqs_PARD6B:
    print("\n description: \n",record.description)
    print("\n record id: \n",record.id)
    print("\n number of features: \n", len(record.features))
    print("\n number of annotations: \n",len(record.annotations))
    for l in record.annotations:
        print("\n *** annotations ***")
        print("\n molecule_type = ",record.annotations["molecule_type"])
        print("\n accessions = ",record.annotations["accessions"])
        print("\n organism = ",record.annotations["organism"])
    #        print("\n comment = ",record.annotations["comment"])
    #        print("\n %s = ",record.annotations["sequence_version"])

records_allseqs_PARD6B=list(SeqIO.parse("allseqs_PARD6B.txt","gb"))

feature_CDS_gene=[]
list_NP=[]
for record in records_allseqs_PARD6B:
    for feat in record.features:
        print("\n *** features ***")
        print(" \n print record.features: \n",feat)
        if feat.type == "CDS":
            feature_CDS_gene.append(feat.qualifiers["translation"])
            print("\n gene: ", feat.qualifiers["gene"])
            print("\n product: ", feat.qualifiers["product"])
            print("\n protein_id: ", feat.qualifiers["protein_id"])
print("\n list of feature_CDS: \n", feature_CDS_gene)
                
protein_PARD6B=feature_CDS_gene[0][0]
print(protein_PARD6B)

protein_PARD6B_seq=open("protein_PARD6B_seq.txt","w")
protein_PARD6B_seq.write(protein_PARD6B)
protein_PARD6B_seq.close()