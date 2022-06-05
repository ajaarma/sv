#!/usr/bin/python

import re,sys,os,gzip
from collections import OrderedDict
global script_path
script_path = os.path.dirname(os.path.abspath( __file__ ))

sys.path.append(script_path+"/CONFIG/")
sys.path.append(script_path+"/SLURM/")

from CONFIG import *
from SLURM import *

def getPARFlag(chr_num,sv_st,sv_en,par_file):
    ''' Check if the input coordinates fall in PAR region '''

    parFlagList = []

    fh = open(par_file)

    for lines in fh:
        lines = lines.strip()
        strs = re.split('\t',lines);strs=[x.strip() for x in strs]
        chr_par = strs[0]
        par_st = strs[1]
        par_en = strs[2]

        if re.search('X',chr_num) and re.search('X',chr_par):
            if sv_st <= par_en and sv_en >= par_st:
                parFlagList.append(True)
        
        if re.search('Y',chr_num) and re.search('Y',chr_par):
            if sv_st <= par_en and sv_en >= par_st:
                parFlagList.append(True)
    
    fh.close()

    parFlagValue = any(parFlagList)



def formatGT(sample_gt,qual,gender,chr_num,sv_id,par_file):
    
    out_gt_strs = []
   
    sv_id_strs = re.split('\+',sv_id)
    try:
        sv_st = int(sv_id_strs[1])
    except:
        print sv_id
        print 'Fatal Error'
        sys.exit()
    try:
        sv_en = int(sv_id_strs[2])
    except:
        sv_en = int(sv_st)

    gt_strs = re.split('-',sample_gt)
    
    for ele in gt_strs:
        if re.search('\/',ele): #All Manta related calls & genotypes
            if gender=='M' and re.search('X|Y',chr_num):
                mn_gtrs = re.split('\:',ele)
                mn_gt = mn_gtrs[0]
                par_flag = getPARFlag(chr_num,sv_st,sv_en,par_file)
                
                if par_flag==False: #Check for non-PAR region overlap
                    if re.search('1/1',mn_gt):
                        mn_gt_str=':'.join(['1',mn_gtrs[1],mn_gtrs[2]])
                    if re.search('0/1',mn_gt):
                        mn_gt_str=':'.join(['1',mn_gtrs[1],mn_gtrs[2]])
                    if re.search('0/0',mn_gt):
                        mn_gt_str=':'.join(['0',mn_gtrs[1],mn_gtrs[2]])
                    
                    out_gt_strs.append(mn_gt_str)
                else:
                    out_gt_strs.append(ele)
            else:
                out_gt_strs.append(ele)

        else: # All Canvas related calls and genotypes
            cv_gtrs = re.split('\:',ele)
            #print cv_gtrs
            cv_gt = float(cv_gtrs[2])
            cv_gt_str=[]
            if not re.search('X|Y',chr_num): #autosomal regions
                if cv_gt==0:
                    cv_gt_str = '1/1'+':'+str(cv_gt)+':'+qual
                elif cv_gt==1:
                    cv_gt_str = '0/1'+':'+str(cv_gt)+':'+qual
                elif cv_gt==2:
                    cv_gt_str = '0/0'+':'+str(cv_gt)+':'+qual
                elif cv_gt >=3:
                    cv_gt_str = '0/1'+':'+str(cv_gt)+':'+qual
            

            # XX: Genotype formatting autosomal style; Canvas calls
            if gender=='F' and re.search('X|Y',chr_num):
                if cv_gt==0:
                    cv_gt_str = '1/1'+':'+str(cv_gt)+':'+qual
                elif cv_gt==1:
                    cv_gt_str = '0/1'+':'+str(cv_gt)+':'+qual
                elif cv_gt==2:
                    cv_gt_str = '0/0'+':'+str(cv_gt)+':'+qual
                elif cv_gt >=3:
                    cv_gt_str = '0/1'+':'+str(cv_gt)+':'+qual

            #XY: Genotype formatting hemizygous; Canvas calls
            if gender=='M' and re.search('X|Y',chr_num):

                # If the Canvas call is in PAR region, GT is autosomal
                par_flag = getPARFlag(chr_num,sv_st,sv_en,par_file)
                
                if par_flag==True:
                    if cv_gt==0:
                        cv_gt_str = '1/1'+':'+str(cv_gt)+':'+qual
                    elif cv_gt==1:
                        cv_gt_str = '0/1'+':'+str(cv_gt)+':'+qual
                    elif cv_gt==2:
                        cv_gt_str = '0/0'+':'+str(cv_gt)+':'+qual
                    elif cv_gt >=3:
                        cv_gt_str = '0/1'+':'+str(cv_gt)+':'+qual
                
                # If the Canvas call is in non-PAR, GT is hemizygous encoded        
                else:
                    if cv_gt==0:
                        cv_gt_str = '1'+':'+str(cv_gt)+':'+qual
                    if cv_gt==1:
                        cv_gt_str = '0'+':'+str(cv_gt)+':'+qual
                    if cv_gt >=2:
                        cv_gt_str = '1'+':'+str(cv_gt)+':'+qual
            
            out_gt_strs.append(cv_gt_str)

    out_gt_strs = [x for x in out_gt_strs if x]
    return '-'.join(out_gt_strs)
         

if __name__=="__main__":
    
    # Get Command line arguments
    objC = CONFIG()
    cmd_dict = objC.parseFormatGTCommandArgs()

    manifest = cmd_dict['manifest'];inp_file=cmd_dict['inp']
    out_file = cmd_dict['out']; fam_id = cmd_dict['fam']
    par_file = cmd_dict['par']; ref_genome = cmd_dict['ref']

    # Get Manifest dictionary
    objS = SLURM()
    famDict = objS.getManifestDict(manifest)
    print famDict[fam_id]
    gender = famDict[fam_id]['gender'][0]

    if famDict[fam_id]['gender'] ==4:
        sib_gender = famDict[fam_id]['gender'][3]
    else:
        sib_gender = []

    # Start processing the input merged file
    fh = gzip.open(inp_file)
    wh = gzip.open(out_file,'w')

    for lines in fh:
        lines = lines.strip()
        #print lines
        if not re.search('^SV_ID',lines):
            strs = re.split('\t',lines)
            strs = [x.strip() for x in strs]

            out_strs = []
            sv_id = strs[0]
            chr_num = strs[1]
            qual = strs[9]
            #k = 26 # If NGC liftover is present
            pb_gt_index = 25
            if re.search('37',ref_genome):
                proband_gt = strs[pb_gt_index]
                mother_gt  = strs[pb_gt_index+1]
                father_gt  = strs[pb_gt_index+2]
                sib_gt     = strs[pb_gt_index+3]
            elif re.search('38',ref_genome):
                proband_gt = strs[pb_gt_index]
                mother_gt  = strs[pb_gt_index+1]
                father_gt  = strs[pb_gt_index+2]
                sib_gt     = strs[pb_gt_index+3]

            out_strs = [sv_id,qual,proband_gt,mother_gt,father_gt,sib_gt]
           
            #print '\t'.join([sv_id,proband_gt,qual,gender])
            proband_gt_fmt = formatGT(proband_gt,qual,gender,chr_num,sv_id,par_file)
            
            if re.search('NA',mother_gt):
                mother_gt_fmt = mother_gt
            else:
                mother_gt_fmt = formatGT(mother_gt,qual,'F',chr_num,sv_id,par_file)
            
            if re.search('NA',father_gt):
                father_gt_fmt = father_gt
            else:
                father_gt_fmt = formatGT(father_gt,qual,'M',chr_num,sv_id,par_file)
           
            if re.search('NA',sib_gt):
                sib_gt_fmt = sib_gt
            else:
                #Sib gender to be extracted from manifest file. TO-DO
                if sib_gender:
                    sib_gt_fmt = formatGT(sib_gt,qual,sib_gender,chr_num,sv_id,par_file)
                else:
                    sib_gt_fmt = sib_gt
           
            strs_fmt = strs
            
            if re.search('grch37',ref_genome):
                strs_fmt[pb_gt_index] = proband_gt_fmt
                strs_fmt[pb_gt_index+1] = mother_gt_fmt
                strs_fmt[pb_gt_index+2] = father_gt_fmt
                strs_fmt[pb_gt_index+3] = sib_gt_fmt
            elif re.search('grch38',ref_genome):
                strs_fmt[pb_gt_index] = proband_gt_fmt
                strs_fmt[pb_gt_index+1] = mother_gt_fmt
                strs_fmt[pb_gt_index+2] = father_gt_fmt
                strs_fmt[pb_gt_index+3] = sib_gt_fmt
      
            print >>wh,'\t'.join(strs_fmt)
        else:
            print >>wh,lines

    fh.close()
    wh.close()
        
             
