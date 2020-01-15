# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 09:32:24 2020

@author: Utilizador
"""

from Bio import Entrez, SeqIO, SeqRecord
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Align.Applications import MuscleCommandline
from Bio.Phylo.Applications import PhymlCommandline
from Bio.Phylo import PhyloXML
from Bio.SeqIO import UniprotIO as BSU#Parse XML Uniprot



###code:
Entrez.email = "pg@uminho.pt"

p1 = open("protein_PARD6B_seq.txt")
p1s = p1.read()
print(p1s)
#atenção, demora mt tempo!!!
#n estou a conseguir usar o file p blast, só handle...
##nr
blast_results_PARD6B=NCBIWWW.qblast("blastp","nr", p1s, hitlist_size= 100)
save_blast_results_PARD6B=open("PARD6B_blast_records.out", "w")
save_blast_results_PARD6B.write(blast_results_PARD6B.read())
save_blast_results_PARD6B.close()

##refseq
blast_ref_results_PARD6B=NCBIWWW.qblast("blastp","refseq_protein", p1s, hitlist_size= 100)
save_blast_ref_results_PARD6B=open("PARD6B_blast_records_refs.out", "w")
save_blast_ref_results_PARD6B.write(blast_ref_results_PARD6B.read())
save_blast_ref_results_PARD6B.close()


###not human
#my_entrez_query1 ="NOT Homo sapiens[ORGN]"
#blast_not_results_PARD6B=NCBIWWW.qblast("blastp","nr", p1s, entrez_query = my_entrez_query1, hitlist_size= 100)
#save_blast_not_results_PARD6B=open("PARD6B_blast_not_records.out", "w")
#save_blast_not_results_PARD6B.write(blast_not_results_PARD6B.read())
#save_blast_not_results_PARD6B.close()

###mouse bd
#my_entrez_query2 ="Mus musculus[ORGN]"
#mblast_results_PARD6B=NCBIWWW.qblast("blastp","nr", p1s,entrez_query = my_entrez_query2, hitlist_size= 100)
#save_mblast_results_PARD6B=open("PARD6B_mblast_records.out", "w")
#save_mblast_results_PARD6B.write(mblast_results_PARD6B.read())
#save_mblast_results_PARD6B.close()

#opens
blast_results_PARD6B=open("PARD6B_blast_records.out")
blast_ref_results_PARD6B=open("PARD6B_blast_records_refs.out")
blast_not_results_PARD6B=open("PARD6B_blast_not_records.out")
mblast_results_PARD6B=open("PARD6B_mblast_records.out")

#reads
blast_records_PARD6B= NCBIXML.read(blast_results_PARD6B)
blast_ref_records_PARD6B= NCBIXML.read(blast_ref_results_PARD6B)
blast_not_records_PARD6B= NCBIXML.read(blast_not_results_PARD6B)
mblast_records_PARD6B= NCBIXML.read(mblast_results_PARD6B)

E_VALUE_THRESH = 0.001
for alignment in blast_records_PARD6B.alignments:
    for hsps in alignment.hsps:
        if hsps.expect < E_VALUE_THRESH:
            print("\n *Alignment*")
            print("sequence:", alignment.title)
            print("accession:", alignment.accession)            
            print("length:", alignment.length)
            print("e value:", hsps.expect)
            print("score:", hsps.score)
            print("identities:", hsps.identities)
            print(hsps.query[0:100] + "...")
            print(hsps.match[0:100] + "...")
            print(hsps.sbjct[0:100] + "...")
            
###read with SearchIO experimental module
from Bio import SearchIO
blast_qresult = SearchIO.read("PARD6B_blast_records.out", "blast-xml")
print(blast_qresult)
blast_hsp = blast_qresult[0][0]    # first hit, first hsp
print(blast_hsp)
###


def get_seqrecs(alignments, threshold):
    for aln in alignments:
        for hsp in aln.hsps:
            if hsp.expect < threshold:
                yield SeqRecord(p1s(hsp.sbjct), id=aln.accession)
                break

best_seqs = get_seqrecs(blast_hsp.alignments, 1e-90)
SeqIO.write(best_seqs, 'PARD6B_family.fasta', 'fasta')

E_VALUE_THRESH = 0.001
for alignment in mblast_records_PARD6B.alignments:
    for hsps in alignment.hsps:
        if hsps.expect < E_VALUE_THRESH:
            print("\n *Alignment*")
            print("sequence:", alignment.title)
            print("accession:", alignment.accession)            
            print("length:", alignment.length)
            print("e value:", hsps.expect)
            print("score:", hsps.score)
            print("identities:", hsps.identities)
            print(hsps.query[0:100] + "...")
            print(hsps.match[0:100] + "...")
            print(hsps.sbjct[0:100] + "...")
            
cmdline=ClustalwCommandline(r'C:\Program Files (x86)\ClustalW2\clustalw2',infile='PARD6B_blast_records_refs.out.txt')
cmdline()

cmdline = PhymlCommandline(cmd='phyml',input='alignprot.phy', datatype='aa', model='WAG', alpha='e', bootstrap=100)
out_log, err_log = cmdline()

muscle_exe = r"C:\Users\Bela\Desktop\Mestrado\Python\Python\Workspace\Labs\muscle3.8.31_i86win32.exe"
cmdline = MuscleCommandline(muscle_exe, input=r"alignprot.fasta", out=r"alignprot.aln",clw=True)
cmdline()