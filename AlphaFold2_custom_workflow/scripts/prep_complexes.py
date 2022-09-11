#!/usr/bin/env python


import sys
import glob
import re

import os
import shutil

from Bio import SeqIO

# shutil.move("path/to/current/file.foo", "path/to/new/destination/for/file.foo")

sA = sys.argv[1]
sB = sys.argv[2]
outdir_path = sys.argv[3]
fasta_dir_path = sys.argv[4]
msaA = sys.argv[5]
msaB = sys.argv[6]


name_A = re.search(".+/(.+)\.fasta", sA).group(1)
name_B = re.search(".+/(.+)\.fasta", sB).group(1)


os.makedirs("{}/complexes/{}__{}/msas".format(outdir_path, name_A, name_B), exist_ok=True)


multi_fasta = open("{}/{}__{}.fasta".format(fasta_dir_path, name_A, name_B), 'w')

chain_id = open("{}/complexes/{}__{}/msas/chain_id_map.json".format(outdir_path, name_A, name_B), 'w')
chain_id.write("{\n\t\"A\": {\n")

A = SeqIO.to_dict(SeqIO.parse(sA, "fasta"))
for i in A:

    multi_fasta.write(">{}\n{}\n".format(A[i].id, A[i].seq))

    chain_id.write("\t\"description\": \"{}\",\n\t\"sequence\": \"{}\"\n".format(A[i].description, A[i].seq))
    chain_id.write("\t},\n\t\"B\": {\n")

B = SeqIO.to_dict(SeqIO.parse(sB, "fasta"))
for i in B:

    multi_fasta.write(">{}\n{}\n".format(B[i].id, B[i].seq))

    chain_id.write("\t\"description\": \"{}\",\n\t\"sequence\": \"{}\"\n".format(B[i].description, B[i].seq))
    chain_id.write("\t}\n}\n")
