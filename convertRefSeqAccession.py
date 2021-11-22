#!/usr/bin/python

import re,sys,os,getopt,gzip
from collections import OrderedDict

script_path = os.path.dirname(os.path.abspath(__file__))

sys.path.append(script_path+"/CONFIG/")

from CONFIG import *

if __name__=="__main__":

    objC = CONFIG()
    cmd_dict = objC.parseRefseqArgs()

    ncDict = OrderedDict()

    inp_file = cmd_dict['inp']
    out_file = cmd_dict['out']
    report_file = cmd_dict['report']

    fh = open(report_file,'r')
    wh = gzip.open(out_file,'wb')

    for lines in fh:
        lines = lines.strip()
        if not re.search('^#',lines):
            strs = re.split('\t',lines)
            strs = [x.strip() for x in strs]
            nc_id = strs[6]
            chrNum = strs[9]
            ncDict[nc_id] = chrNum

    fh.close()

    fh = gzip.open(inp_file,'rb')

    for lines in fh:
        lines = lines.strip()
        if re.search('^#',lines):
            print >>wh,lines
        else:
            strs = re.split('\t',lines)
            strs = [x.strip() for x in strs]
            nc_id = strs[0]
            chrNum = ncDict[nc_id]
            strs[0] = chrNum
            nc_id = strs[0]
            print >>wh, '\t'.join(strs)
    
    wh.close()
        


