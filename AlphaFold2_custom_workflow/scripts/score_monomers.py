#!/usr/bin/env python

import pandas as pd
import numpy as np
import itertools
import glob
import json
import glob
import sys
import re
import os

outdir_A = sys.argv[1]
outdir_B = sys.argv[2]
summary_path = sys.argv[3]
fasta_dir = sys.argv[4]
repo_msas = sys.argv[5]
repo_multimer = sys.argv[6]

def scoring(outdir):
    scoring = []

    for files in glob.glob("{}/*/ranking_debug.json".format(outdir)):
        name = re.search(".+\/(.+)\/ranking_debug.json", files).group(1)

        f = open(files)
        data = json.load(f)

        [scoring.append([name, x, data['plddts'].get(x)]) for x in data['plddts']]

    scoring = pd.DataFrame(scoring, columns=['prot_id', 'model', 'plddts'])
#    ids = np.unique(scoring[scoring.plddts > 65].prot_id.values)
    ids = np.unique(scoring[scoring.plddts > 1].prot_id.values)

    return(scoring, ids)

def remove_low(df, ids):
    scoring = df.copy()

    for name in np.unique(scoring.prot_id.values):
        if name not in ids:
            file_list = glob.glob("{}/*{}*.fasta".format(fasta_dir, name))

            for file_path in file_list:
                os.remove(file_path)

def link_dbs(msas_path, multimer_path, comb_list, comb_dir, abbrv):
    os.makedirs("{}/{}_{}/msas/{}".format(multimer_path, comb_list[0], comb_list[1], abbrv, exist_ok=True))
    link_files = glob.glob("{}/{}/msas/*".format(msas_path, comb_dir))
    for files in link_files:
        db_name = re.search(".*/([a-zA-Z_\.0-9]+)", files).group(1)
        os.symlink(files, "{}/{}_{}/msas/{}/{}".format(multimer_path, comb_list[0], comb_list[1], abbrv, db_name))


scoring_A, ids_A = scoring(outdir_A)
scoring_B, ids_B = scoring(outdir_B)

remove_low(scoring_A, ids_A)
remove_low(scoring_B, ids_B)

comb = list(itertools.product(ids_A, ids_B))

out_json = {}
out_json["samples"] = []

for i in comb:
    out_json["samples"].append("{}_{}".format(i[0], i[1]))

    os.makedirs("{}/{}_{}/msas".format(repo_multimer, i[0], i[1]), exist_ok=True)

    if os.path.isdir("{}/{}_{}/msas/A".format(repo_multimer, i[0], i[1])) == False:
        link_dbs(repo_msas, repo_multimer, i, i[0], "A")

    if os.path.exists("{}/{}_{}/msas/B".format(repo_multimer, i[0], i[1])) == False:
        link_dbs(repo_msas, repo_multimer, i, i[1], "B")


with open("config_multimer_complexes.json", "w") as outfile:
    json.dump(out_json, outfile, indent=4)

scoring_A.to_csv("{}/monomers_scoring_A.tsv".format(summary_path), sep="\t")
scoring_B.to_csv("{}/monomers_scoring_B.tsv".format(summary_path), sep="\t")
