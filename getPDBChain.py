#!/bin/env python
# vim: set fileencoding=utf-8 :
#
# Author:   toshi123
# License:  MIT License
# Created:  2014-05-29
#

import sys
import argparse
from prody import *

if __name__ == '__main__':
    confProDy(verbosity='none')

    p = argparse.ArgumentParser(description="Download and Write PDB file")
    p.add_argument('id')
    p.add_argument('-c','--chain',default='all',help='chain name')
    p.add_argument('-f','--fasta',default=False,action="store_true",help='output fasta file format')
    p.add_argument('--ca',default=False,action="store_true",help='output only CA atoms')
    args = p.parse_args()

    pdbfile = fetchPDB(args.id)
    if pdbfile is None:
        sys.exit("Fatal: "+args.id+" dosen't exist")

    if args.chain == "all":
        pdb = parsePDB(pdbfile)
        name = args.id
    else:
        pdb = parsePDB(pdbfile, chain=args.chain)
        name = args.id+args.chain
        if pdb is None:
            sys.exit("Fatal: "+args.chain+" chain dosen't exist in "+args.id)

    if args.fasta:
        fastafile = open(name+".faa","w")
        for chain in pdb.iterChains():
            title = args.id+chain.getChid()
            seq = chain.getSequence()
            if seq is not None:
                fastafile.write(">"+title+"\n"+seq+"\n")
            else:
                print "chain "+chain.getChid()+" is not protein"
        fastafile.close()
    elif args.ca:
        writePDB(name+"_ca",pdb.select('ca'))
    else:
        writePDB(name,pdb)

