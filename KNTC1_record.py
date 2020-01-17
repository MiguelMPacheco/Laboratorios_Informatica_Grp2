# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 11:33:23 2020

@author: Utilizador
"""


from Bio import Entrez
from Bio import Medline
from Bio import SeqIO
    

handle_pKNTC1 =Entrez.efetch(db="protein",id="NP_055523.1",rettype = "gb",retmode="text")
save_pKNTC1=open("pKNTC1.txt","w")
save_pKNTC1.write(handle_pKNTC1.read())
save_pKNTC1.close()

records_pKNTC1=SeqIO.parse("pKNTC1.txt","gb")


for record in records_pKNTC1:
    print("\n description: \n",record.description)
    print("\n record id: \n",record.id)
    print("\n number of features: \n", len(record.features))
    print("\n number of annotations: \n",len(record.annotations))
    for l in record.annotations:
        print("\n *** annotations ***")
        print("\n molecule_type = ",record.annotations["topology"])
        print("\n accessions = ",record.annotations["accessions"])
        print("\n organism = ",record.annotations["source"])
        print("\n comment = ",record.annotations["comment"])
        print("\n %s = ",record.annotations["structured_comment"])

records_pKNTC1=list(SeqIO.parse("pKNTC1.txt","gb"))

feature_CDS_gene=[]
for record in records_pKNTC1:
        for feat in record.features:
            print("\n *** features ***")
            print(" \n print record.features: \n",feat)
            print("\n qualifiers: \n", feat.qualifiers.items())
