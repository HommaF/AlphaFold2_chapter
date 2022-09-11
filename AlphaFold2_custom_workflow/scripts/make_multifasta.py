#!/usr/bin/env python

import sys
import glob
import re

import os
import shutil

from Bio import SeqIO

sA = sys.argv[1]
sB = sys.argv[2]
fasta_dir_path = sys.argv[3]
outdir_path = sys.argv[4]
repo_msas = sys.argv[5]

name_A = re.search(".+/(.+)\.fasta", sA).group(1)
name_B = re.search(".+/(.+)\.fasta", sB).group(1)


def link_to_af_monomers(path_out, path_repo, name, abbrv):
    if os.path.exists("{}/{}/{}".format(path_out, abbrv, name)) == False:
            os.makedirs("{}/{}".format(path_out, abbrv), exist_ok = True)
            os.symlink("{}/{}".format(path_repo, name), "{}/{}/{}".format(path_out, abbrv, name))

link_to_af_monomers(outdir_path, repo_msas, name_A, "A")
link_to_af_monomers(outdir_path, repo_msas, name_B, "B")


if os.path.isdir("{}".format(fasta_dir_path)) == False:
    os.makedirs("{}".format(fasta_dir_path), exist_ok = True)

multi_fasta = open("{}/{}_{}.fasta".format(fasta_dir_path, name_A, name_B), 'w')

A = SeqIO.to_dict(SeqIO.parse(sA, "fasta"))
for i in A:
    multi_fasta.write(">{}\n{}\n".format(A[i].id, A[i].seq))

B = SeqIO.to_dict(SeqIO.parse(sB, "fasta"))
for i in B:
    multi_fasta.write(">{}\n{}\n".format(B[i].id, B[i].seq))
