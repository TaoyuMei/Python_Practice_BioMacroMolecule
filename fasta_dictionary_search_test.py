#!/usr/bin/python3

import handin5
fasta_dict=handin5.read_fasta("Ecoli.prot.fasta")
print(len(fasta_dict.keys()))

yhcn=handin5.find_prot(fasta_dict,"YHCN_ECOLI")

boom=handin5.find_prot(fasta_dict,"BOOM_ECOLI")

matches=handin5.find_prot2(fasta_dict,"^.{3}_ECOLI")
print(len(matches))