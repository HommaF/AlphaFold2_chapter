#!/usr/bin/env python

## Database search as in AlphaFold v2.1.1, using the reduced dataset and including all databases needed to model multimers

import os
import subprocess
import shutil
import sys
import re


data_dir = sys.argv[1]
fasta_loc= sys.argv[2]
out_dir = sys.argv[3]
n_cpu = sys.argv[4]

#name = re.search('.+/([A-Za-z0-9\.]+)\.fasta', fasta_loc).group(1)
#print(name)

#msa_dir = os.path.join(out_dir, name, 'msas')
#os.makedirs(msa_dir)
#msa_file = os.path.join(msa_dir, db_name)

hhblits_binary = shutil.which('hhblits')
n_iter = 3
e_value = 0.001
maxseq = 1000000
realign_max = 100000
maxfilt = 100000
min_prefilter_hits = 1000
p = 20
z = 500
# a3m_path??


os.system('{} -i {} -cpu {} -oa3m {} -o /dev/null -e {} -maxseq {} -realign_max {} -maxfilt {} -min_prefilter_hits {}'.format(hhblits_binary, fasta_loc, n_cpu, a3m_path, e_value, maxseq, realign_max, maxfilt, min_prefilter_hits))
