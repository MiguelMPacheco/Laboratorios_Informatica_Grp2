# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 20:36:19 2020

@author: Utilizador
"""

from Bio import Entrez
from Bio import Medline
from Bio import SeqIO
    
    
Entrez.email = "pg@uminho.pt"

handle_seq_search_KNTC1 = Entrez.esearch(db="nucleotide", term="Homo sapiens[Organism] AND KNTC1[Gene Name]", idtype="acc")
record_seq_KNTC1 = Entrez.read(handle_seq_search_KNTC1)
print("\n number of KNTC1 ids associated with chosen term: \n",record_seq_KNTC1["Count"])
print("\n list of KNTC1 ids associated with chosen term: \n",record_seq_KNTC1["IdList"])

handle_allseqs_KNTC1 =Entrez.efetch(db="nucleotide",id=','.join(record_seq_KNTC1["IdList"]),rettype = "gb",retmode="text")
save_allseqs_KNTC1=open("allseqs_KNTC1.txt","w")
save_allseqs_KNTC1.write(handle_allseqs_KNTC1.read())
save_allseqs_KNTC1.close()

handle_allseqs_KNTC1 =Entrez.efetch(db="nucleotide",id=','.join(record_seq_KNTC1["IdList"]),rettype = "gb",retmode="text")
for index, record in enumerate(SeqIO.parse(handle_allseqs_KNTC1,"gb")):
    print("index %i, ID = %s, length %i, with %i features and %i annotations"
          % (index, record.id, len(record.seq), len(record.features), len(record.annotations)))

records_allseqs_KNTC1=SeqIO.parse("allseqs_KNTC1.txt","gb")



#falta dict={id:seq}

for record in records_allseqs_KNTC1:
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

records_allseqs_KNTC1=list(SeqIO.parse("allseqs_KNTC1.txt","gb"))

feature_CDS_gene=[]
list_NP=[]
for record in records_allseqs_KNTC1:
    if record.id[:2] == "NG":
        for feat in record.features:
#            print("\n *** features ***")
#            print(" \n print record.features: \n",feat)
            if feat.type == "CDS" and ("KNTC1" in feat.qualifiers["gene"]):
                feature_CDS_gene.append(feat.qualifiers["translation"])
                print("\n gene: ", feat.qualifiers["gene"])
                print("\n product: ", feat.qualifiers["product"])
                print("\n protein_id: ", feat.qualifiers["protein_id"])
                print("\n db_xref: ", feat.qualifiers["db_xref"])
        print("\n list of feature_CDS: \n", feature_CDS_gene)

                
protein_KNTC1=feature_CDS_gene[0][0]
print(protein_KNTC1)

protein_KNTC1_seq=open("protein_KNTC1_seq.txt","w")
protein_KNTC1_seq.write(protein_KNTC1)
protein_KNTC1_seq.close()