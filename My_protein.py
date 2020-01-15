# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 09:33:54 2020

@author: Utilizador
"""

class My_protein:
    
    from Bio import Entrez
    from Bio import SeqIO
    from Bio import ExPASy
    from Bio import SwissProt
    
    def __init__(self,

###code:

###a entrada dos ids tem de vir automatizada...

#retrieve data from Uniprot
handle = ExPASy.get_sprot_raw("Q6P4I2")
record = SwissProt.read(handle)
print(record)
print(record.description)
print(record.references)
for i in references:
    print(i)
for ref in record.references:
    print("authors:", ref.authors)
    print("title:", ref.title)
for record in record.annotations:
    print(record)
print(record.annotations["date"])
#dir(record)


handle = ExPASy.get_prosite_raw('PS00001') 
record = Prosite.read(handle) 
handle.close() 
print(record.accession) 

##retrive data from NCBI CDD
handle=Entrez.efetch(db="protein", rettype="gb", retmode="text", id="NP_116245")
protein_record=SeqIO.read(handle, "genbank")
SeqIO.write(protein_record,"file_name.gbk","genbank")
handle.close()
for i in protein_record.features:
    print(i)
    print(i.qualifiers)