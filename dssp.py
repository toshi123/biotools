#!/bin/env python
# vim: set fileencoding=utf-8 :
#
# Author:   toshi123
# License:  MIT License
# Created:  2014-05-09
#

import sys
from Bio.PDB.DSSP import *
from prody import *
import tempfile

one2three = {\
        "W":"TRP",\
        "Y":"TYR",\
        "F":"PHE",\
        "I":"ILE",\
        "L":"LEU",\
        "M":"MET",\
        "V":"VAL",\
        "H":"HIS",\
        "R":"ARG",\
        "Q":"GLN",\
        "E":"GLU",\
        "T":"THR",\
        "N":"ASN",\
        "K":"LYS",\
        "A":"ALA",\
        "P":"PRO",\
        "D":"ASP",\
        "C":"CYS",\
        "S":"SER",\
        "G":"GLY"\
    }

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
        asa[num] = float(dsspresult[0][key][2]) / MAX_ACC[one2three[dsspresult[0][key][0].upper()]]
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

