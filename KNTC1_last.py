# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 20:36:19 2020

@author: rjoana
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
        print("\n comment = ",record.annotations["comment"])
        print("\n %s = ",record.annotations["sequence_version"])

records_allseqs_KNTC1=list(SeqIO.parse("allseqs_KNTC1.txt","gb"))

feature_CDS_gene=[]

for record in records_allseqs_KNTC1:
    if record.id[:2] == "NG":
        for feat in record.features:
            print("\n *** features ***")
            print(" \n print record.features: \n",feat)
            if feat.type == "source":
                print(feat.qualifiers.keys())
                print("",feat.qualifiers["chromosome"])
                print("",feat.qualifiers["map"])            
            if feat.type == "CDS" and ("KNTC1" in feat.qualifiers["gene"]):
                feature_CDS_gene.append(feat.qualifiers["translation"])
                print("\n gene: ", feat.qualifiers["gene"])
                print("\n location and strand: \n", feat.location)
                print("\n product: ", feat.qualifiers["product"])
                print("\n protein_id: ", feat.qualifiers["protein_id"])
                print("\n db_xref: ", feat.qualifiers["db_xref"])
        print("\n list of feature_CDS: \n", feature_CDS_gene)

                
protein_KNTC1=feature_CDS_gene[0][0]
print(protein_KNTC1)

protein_KNTC1_seq=open("protein_KNTC1_seq.txt","w")
protein_KNTC1_seq.write(protein_KNTC1)
protein_KNTC1_seq.close()


from Bio.Blast import NCBIXML
from Bio.Blast import NCBIWWW

pKNTC1o = open("protein_KNTC1_seq.txt")
pKNTC1r = pKNTC1o.read()
print(pKNTC1r)


##nr
blast_results_KNTC1=NCBIWWW.qblast("blastp","nr", pKNTC1r, hitlist_size= 100)
save_blast_results_KNTC1=open("KNTC1_blast_records.xml", "w")
save_blast_results_KNTC1.write(blast_results_KNTC1.read())
save_blast_results_KNTC1.close()

blast_results_KNTC1=open("KNTC1_blast_records.xml")
blast_records_KNTC1= NCBIXML.read(blast_results_KNTC1)

E_VALUE_THRESH = 0.001
for alignment in blast_records_KNTC1.alignments:
    for hsps in alignment.hsps:
        if hsps.expect < E_VALUE_THRESH:
            print("\n *KNTC1 nr Alignment*")
            print("sequence:", alignment.title)
            print("accession:", alignment.accession)            
            print("length:", alignment.length)
            print("e value:", hsps.expect)
            print("score:", hsps.score)
            print("identities:", hsps.identities)
            print(hsps.query[0:100] + "...")
            print(hsps.match[0:100] + "...")
            print(hsps.sbjct[0:100] + "...")


my_entrez_query ="Mus musculus[ORGN]"
mblast_results_KNTC1=NCBIWWW.qblast("blastp","nr", pKNTC1r.format("fasta"),entrez_query = my_entrez_query, hitlist_size= 100)
save_mblast_results_KNTC1=open("KNTC1_mblast_records.xml", "w")
save_mblast_results_KNTC1.write(mblast_results_KNTC1.read())
save_mblast_results_KNTC1.close()

mblast_results_KNTC1=open("KNTC1_mblast_records.xml")
mblast_records_KNTC1= NCBIXML.read(mblast_results_KNTC1)

E_VALUE_THRESH = 0.001
for alignment in blast_records_KNTC1.alignments:
    for hsps in alignment.hsps:
        if hsps.expect < E_VALUE_THRESH:
            print("\n *KNTC1 mouse Alignment*")
            print("sequence:", alignment.title)
            print("accession:", alignment.accession)            
            print("length:", alignment.length)
            print("e value:", hsps.expect)
            print("score:", hsps.score)
            print("identities:", hsps.identities)
            print(hsps.query[0:100] + "...")
            print(hsps.match[0:100] + "...")
            print(hsps.sbjct[0:100] + "...")
       
