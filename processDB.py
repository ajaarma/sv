#!/usr/bin/python


############################################################################################################
#
# Description : Generic script to process CLINVAR related VCF files: nstd45, nstd51, 
#                           nstd101 for both GRCh37, GRCh38 reference genome built
#
#
# Usage: python processCLINVAR.py -i <input-file> -o <output-file-string>
#
#
#
#############################################################################################################

import re,sys,os,getopt
from vcf_parser import VCFParser

script_path = os.path.dirname(os.path.abspath(__file__))

sys.path.append(script_path+"/CONFIG/")
sys.path.append(script_path+"/DB/")

from CONFIG import *
from DB import *


if __name__=="__main__":


    #manifest_file,xml_file,db_type,ref_genome,lift_ov_flag,
    #proj_date = processArgs(sys.argv[1:])

    objC = CONFIG()
    cmd_dict = objC.parseDBCommandArgs()
    try:
        manifest_file = cmd_dict['manifest'][0];
    except:
        manifest_file = None
    xml_file = cmd_dict['analysis']
    ref_genome = cmd_dict['ref'];db_name = cmd_dict['db']
    try:
        lift_ov_flag = cmd_dict['lov'][0]
    except:
        lift_ov_flag = None
    try:
        proj_date = cmd_dict['project'][0]
    except:
        proj_date = None


    objC = CONFIG()
    objD = dbSV()

    configDict = objC.getConfigDict(xml_file)

    # Assign boolean True/False to family overlap fraction option
    a1 = configDict['overlapMerge']['famOverlapFrac']
    configDict['overlapMerge']['famOverlapFrac'] = objC.str2bool(a1)

    objD.displayArguments(configDict,db_name,ref_genome,manifest_file,
                                                proj_date,lift_ov_flag)
    
    

