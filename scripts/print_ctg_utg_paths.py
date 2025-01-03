# -*- coding: utf-8 -*-
#================================================
#
#         Author: Haoyu Cheng
#    Description: -
#
#================================================
from __future__ import print_function
import sys
import os

def take_ref_id(rname, rid, idx):
    zs = idx.split('@')
    for z in zs:
        s = z.split(':')
        if s[1] == rid:
            return s[0]
    return zs[0].split(':')[0]

nums = len(sys.argv)
if nums < 3:
    print("No. of parameter must >= 2 ...")
    print("python3 print_ctg_utg_paths.py query.p_ctg.gfa ref.p_utg.gfa")
    os._exit(0)

query_contig = sys.argv[1]
ref_utg = sys.argv[2]

# First pass: collect orientations from L lines
orientations = {}
with open(ref_utg, 'r') as f1:
    for line in f1:
        if line[0] == 'L':
            s = line.strip().split('\t')
            orientations[s[1]] = s[2]
            orientations[s[3]] = s[4]

# Second pass: collect utg mappings
dict_ref = {}
with open(ref_utg, 'r') as f1:
    for line in f1:
        if line[0] == 'A':
            s = line.strip().split('\t')
            utg_info = s[1] + ':' + s[7].split(':')[2]
            if s[4] not in dict_ref:
                dict_ref[s[4]] = utg_info
            else:
                dict_ref[s[4]] = dict_ref[s[4]] + '@' + utg_info

# Process query contig
list_contig = {}
with open(query_contig, 'r') as f1:
    for line in f1:
        if line[0] == 'A':
            s = line.strip().split('\t')
            if s[4] in dict_ref:
                ref_name = take_ref_id(s[4], s[7].split(':')[2], dict_ref[s[4]])
                orientation = '>' if ref_name in orientations and orientations[ref_name] == '+' else '<'
                if s[1] not in list_contig:
                    list_contig[s[1]] = orientation + ref_name
                else:
                    strl = len(ref_name)
                    if strl <= len(list_contig[s[1]]):
                        if list_contig[s[1]][-strl:] == ref_name:
                            continue
                    list_contig[s[1]] = list_contig[s[1]] + orientation + ref_name

for key in list_contig:
    print(f"{key}\t{list_contig[key]}")
