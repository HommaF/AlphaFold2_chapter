#!/usr/bin/env python

import os
import re
import sys
import glob

target_dir = sys.argv[1]

dont_process = ['msas', 'ranking_debug.json']

for file_path in glob.glob(os.path.join(target_dir, "*")):
    name = re.search(".*\/([a-zA-Z0-9_\.]+)", file_path).group(1)

    if name not in dont_process:
        os.system("gzip -9 {}".format(file_path))
