#!/usr/bin/env python

import os
import shutil
import re
import sys
import glob

parent_repo = sys.argv[1]
pkl_repo = sys.argv[2]

sample = re.search('.+\/([a-zA-Z0-9_\.]+)', parent_repo).group(1)
print(sample)

os.makedirs(os.path.join(pkl_repo, sample), exist_ok=True)

for file_path in glob.glob(os.path.join(parent_repo, "result*pkl*")):
    shutil.move(file_path, os.path.join(pkl_repo, sample))
