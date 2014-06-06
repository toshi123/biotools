#!/bin/env python
# vim: set fileencoding=utf-8 :
#
# Author:   toshi123
# License:  MIT License
# Created:  2014-05-09
#

import sys
from Bio.PDB.DSSP import *
from Bio.PDB.Polypeptide import *
from prody import *
import tempfile

def calcASA(pdbfile,dssp="dssp"):
    dsspresult = typeselecter(pdbfile,dssp)
    asa = {}
    for key in dsspresult[1]:
        num = int(key[1][1])
        asa[num] = int(dsspresult[0][key][2])
    return asa

def calcNomalizedASA(pdbfile,dssp="dssp"):
    dsspresult = typeselecter(pdbfile,dssp)
    asa = {}
    for key in dsspresult[1]:
        num = int(key[1][1])
        three = one_to_three(dsspresult[0][key][0].upper())
        asa[num] = float(dsspresult[0][key][2]) / MAX_ACC[three]
    return asa


def runDSSP(pdbfile,dssp="dssp"):
    dsspresult = typeselecter(pdbfile,dssp)
    result = {}
    for key in dsspresult[1]:
        num = int(key[1][1])
        result[num] = {\
            "residue":dsspresult[0][key][0],\
            "secondary":dsspresult[0][key][1],\
            "asa":int(dsspresult[0][key][2]),\
            "phi":float(dsspresult[0][key][3]),\
            "psi":float(dsspresult[0][key][4])\
        }
    return result

def typeselecter(pdbfile,dssp="dssp"):
    if isinstance(pdbfile,(prody.atomic.atomgroup.AtomGroup, prody.atomic.selection.Selection)):
        tmp = tempfile.NamedTemporaryFile(suffix='.pdb')
        writePDB(tmp.name,pdbfile)
        dsspresult = dssp_dict_from_pdb_file(tmp.name,DSSP=dssp)
        tmp.close()
    elif isinstance(pdbfile,basestring):
        dsspresult = dssp_dict_from_pdb_file(pdbfile,DSSP=dssp)
    else:
        sys.exit("dssp.calcASA requires PDB object or file name")

    return dsspresult

