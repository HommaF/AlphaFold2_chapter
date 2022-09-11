#!/usr/bin/env python

import pandas as pd
import numpy as np
import itertools
import json
import glob
import sys
import re
import os

outdir = sys.argv[1]
summary_path = sys.argv[2]

def scoring(outdir):
    scoring = []

    for files in glob.glob("{}/*/ranking_debug.json".format(outdir)):
        name = re.search(".+\/(.+)\/ranking_debug.json", files).group(1)
        names = re.search("(.+)_(.+)", name)
        name_A = names.group(1)
        name_B = names.group(2)

        f = open(files)
        data = json.load(f)

        [scoring.append([name, name_A, name_B, x, data['iptm+ptm'].get(x)]) for x in data['iptm+ptm']]

    scoring = pd.DataFrame(scoring, columns=['prot_id', 'A', 'B', 'model', 'iptm_ptm'])
    ids = np.unique(scoring[scoring.iptm_ptm > 0.65].prot_id.values)

    return(scoring, ids)


scoring, ids = scoring(outdir)

scoring.to_csv("{}/multimers_scoring.tsv".format(summary_path), sep='\t')
