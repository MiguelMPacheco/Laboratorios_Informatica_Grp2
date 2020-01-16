# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 09:30:16 2020

@author: Utilizador
"""




from Bio import Entrez
from Bio import Medline
from Bio import SeqIO

Entrez.email = "pg@uminho.pt"

keyword = ["human AND WDR73 AND function", "human AND PARD6B AND function","human AND KNTC1 AND function"]
for i in keyword:
    x = open("file_"+str(i)+".txt","w+")
    handle_search = Entrez.esearch(db="pubmed", term = i)
    record_search = Entrez.read(handle_search)
    print("\n number of records found from Pubmed: \n", i, record_search["Count"])
    print("\n ids records found from Pubmed: \n ", i, record_search["IdList"])
    id_list = record_search["IdList"]
    print(len(id_list))
    handle_ids = Entrez.efetch(db="pubmed", id = id_list , rettype = "medline", idtype="text")
    records_pubmed= Medline.parse(handle_ids)
    for record in records_pubmed:
        print("\n title: \n",record["TI"], file = x)
        print("\n abstract: \n", record["AB"], file = x)
        print("\n abstract: \n", record["CRDT"], file = x)
        print("\n title: \n",record["TI"])
        print("\n abstract: \n", record["AB"])
        print("\n abstract: \n", record["CRDT"])
    x.close()