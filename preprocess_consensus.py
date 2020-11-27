import collections
import pysam
import sys
import argparse
import os
import re

parser = argparse.ArgumentParser(prog='preprocess_consensus.py',description="preprocess consensus sequences for covid phylogenetic analysis")
parser.add_argument('-d','--directory',dest="dir",help='directory of fasta consensus files',required=True)
parser.add_argument('-c','--completeness',dest='threshold',help='completeness threshold',type=float,default=0.75)
parser.add_argument('-t','--table',dest='table',help="translation table. First column is current name, Second column is new name")
parser.add_argument('-x','--exclude',dest='exclude',help="file with list of ids to exclude")
parser.add_argument('-i','--include',dest='include',help="file with list of ids to include")
parser.add_argument('-o','--output',help='output directory', dest="out",required=True)



args=parser.parse_args()
print(args.out)
print(args.dir)

log=open(args.out + "/preprocess.log","w")
cons=open(args.out + "/consensus.fasta","w")

tdict={}
if args.table:
    with open(args.table) as f:
        for line in f:
            key,value=line.strip().split("\t")
            tdict[key]=value
    
xlist=[]
if args.exclude:
    with open(args.exclude) as f:
        for line in f:
            key=line.strip()
            xlist.append(key)

ilist=[]
if args.include:
    with open(args.include) as f:
        for line in f:
            key=line.strip()
            ilist.append(key)


for fn in os.listdir(args.dir):
    fpath=args.dir + "/" + fn
    print(fpath)
    fa = pysam.FastxFile(fpath)

    for record in fa:
        ### if there is a translation dictionary, then replace the name.
        ### otherwise, cleanup the name
        if(record.name in tdict):
            record.name=tdict[record.name]
        else:
            record.name = record.name.replace("Consensus_","")
            record.name = re.sub('.primertrimmed.*','',record.name)
        
        # if there is an include list, then skipif not in the include list
        if( (len(ilist)>0) and (record.name not in ilist)):
            continue
        # if there is an exclude list, then skip if in the exclude llist
        if(record.name in xlist):
            continue
            
        base_counter = collections.Counter()
        for b in record.sequence:
            base_counter.update(b.upper())

        total = 0
        for base, count in base_counter.items():
            total += count
        completeness = 0

        if total > 0:
            completeness = 1 - (float(base_counter['N']) / total)
        if completeness >= args.threshold:
            cons.write(">" + record.name + "\n")
            cons.write(record.sequence + "\n")

        log.write(fn + "\t" + record.name + "\t" + str(completeness) + "\n")

