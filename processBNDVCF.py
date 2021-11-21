#!/usr/bin/python

import re,sys,os,gzip
import getopt

def processArgs(inp_args):

    opts = {}
    inp_file = [];out_file=[]

    print inp_args
    try:
        opts,args = getopt.getopt(inp_args,'hi:o:',['help'])
    except getopt.GetoptError as err:
        usage(err)

    for opt,arg in opts:
        print opt,arg
        if opt =="-h" or opt =="-help":
            usage()
        if opt in ("-i","--input"):
            inp_file = arg
        if opt in ("-o","--output"):
            out_file = arg

    print "\n"
    print inp_file
    if inp_file and out_file:
        print "Entered input file is: ",inp_file
        print "Entered output file is: ",out_file
        print "\n"
        return inp_file,out_file
    else:
        usage("Some of the arguments are missing. Please check the correct format of usage below\n")

    

def usage(err=[]):
    if err:
        print "\n"
        print err
        print "python processBNDVCF.py -i <input-raw-vcf> -o <out-file-name>"
        print "\n"
        sys.exit()
    else:
        print "\n"
        print "python processBNDVCF.py -i <input-raw-vcf> -o <out-file-name> -f"
        print "\n"
        sys.exit()

def getEndCoord(info_strs):

    end_pos = []
    for e in info_strs:
        if re.search("^END",e):
            end_pos = re.split("\=",e)[1]
            break

    return end_pos


if __name__=="__main__":

    inp_file,out_file = processArgs(sys.argv[1:])

    fh = gzip.open(inp_file)
    wh = gzip.open(out_file,"wb")

    for lines in fh:
        lines = lines.strip()
        if re.search("^#",lines):
            print >>wh,lines
        else:
            if re.search("MantaBND",lines):
                strs = re.split("\t",lines)
                strs = [x.strip() for x in strs]

                chr_num = strs[0]
                pos = int(strs[1])
                sv_id = strs[2]
                ref = strs[3]
                alt_id = strs[4]
                alt_id = "<INS:"+alt_id+">"
                info_strs = re.split("\;",strs[7])
                sv_len = "END="+str(pos+50)
                
                info_strs.insert(0,sv_len)
                #print strs[4]
                #print info_strs
                qual = strs[5]
                filter_tag = strs[6]
                gt_tag = strs[8]
                gt_val = strs[9]

                pos_50 = pos-50
                out_line = [chr_num,str(pos_50),sv_id,ref,alt_id,qual,filter_tag,
                                               ";".join(info_strs),gt_tag,gt_val
                           ]
                print >>wh,"\t".join(out_line)

            elif re.search("MantaINS",lines):
                strs = re.split("\t",lines)
                strs = [x.strip() for x in strs]

                chr_num = strs[0]
                pos = int(strs[1])
                sv_id = strs[2]
                ref = strs[3]
                alt_id = strs[4]
                if not re.search("INS",alt_id):
                    alt_id = "<INS:"+alt_id+">"
                info_strs = re.split("\;",strs[7])
                end_pos = getEndCoord(info_strs)
                
                if int(end_pos) == int(pos):
                    sv_len = "END="+str(pos+50)
                else:
                    sv_len = "END="+str(int(end_pos)+50)
                
                for k,e in enumerate(info_strs):
                    if re.search("^END",e):
                        info_strs.pop(k)
                        info_strs.insert(k,sv_len)

                info_strs.insert(0,sv_len)
                #print strs[4]
                #print info_strs
                qual = strs[5]
                filter_tag = strs[6]
                gt_tag = strs[8]
                gt_val = strs[9]

                pos_50 = pos-50
                out_line = [chr_num,str(pos_50),sv_id,ref,alt_id,qual,filter_tag,
                                            ";".join(info_strs),gt_tag,gt_val
                           ]
                print >>wh,"\t".join(out_line)
            else:
                print >>wh,lines

    wh.close()
