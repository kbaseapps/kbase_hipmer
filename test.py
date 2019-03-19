#!/usr/bin/env python
import glob
for name in glob.glob('dir/*'):
    print name

m = glob.glob('dir/*.txt')[0]
print("m {}".format(m))
