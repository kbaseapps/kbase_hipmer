#!/usr/bin/env python
import sys,os
output_contigs=''
name="final_assembly.fa"
for root, dirs, files in os.walk("."):
#    print("FILES={}".format(files))
    if name in files:
        output_contigs = os.path.join(root, name)

if output_contigs:
    print("OUTPUT CONTIGS ARE HERE {}".format(output_contigs))
else:
    print("output contigs not found")
    sys.exit()
