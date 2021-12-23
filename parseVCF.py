#!/usr/bin/python

###########################################################################################
#
# Description: Genreic script to process the VCFAnno output. Currently incorporates these annotation sources
#               -- GNOMAD v2.1
#               -- NGC
#               -- ENSEMBL
#               -- CLINVAR (nstd45,nstd51,nstd102)
#
# Usage: python parseVCF.py <vcfanno-annotated-input> <output-file>
#
#
###########################################################################################


import re,sys,gzip,os,getopt,time
import argparse

start = time.time()
script_path = os.path.dirname(os.path.abspath( __file__ ))

sys.path.append(script_path+"/CONFIG/")
sys.path.append(script_path+"/VCFANNO/")
sys.path.append(script_path+"/SLURM/")

from CONFIG import *
from VCFANNO import *
from SLURM import *

if __name__=="__main__":

    objC = CONFIG()
    inp_file, ngc_id, line_index, manifest_file,xml_file,\
    db_type, task, debug_flag, ovFrac = objC.parseVCFCommandArgs()
    print "\n", inp_file, ngc_id, line_index, manifest_file,xml_file, db_type,\
                                                 task, debug_flag, ovFrac,"\n"

    objC = CONFIG()
    configDict = objC.getConfigDict(xml_file)
    objS = SLURM()
    famDict = objS.getManifestDict(manifest_file)

    objV = VCFANNO()
    #header,wh = objV.writeDBSpecificHeader("gnomad",wh)
    objV.processAnnotatedVCF(configDict,inp_file,ngc_id,db_type,famDict,line_index,
                                                debug_flag,float(ovFrac),"overlap")

    end = time.time()
    print "The script tool approx: ",end-start," seconds"
    sys.exit()

