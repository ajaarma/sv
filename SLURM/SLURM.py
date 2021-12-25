#!/usr/bin/python

########################################################################################################
#
# Description: Generic class for creating the SLURM script for annotating Structural Variants in HPC
#
#
#
#
#
########################################################################################################

import re,sys,os,subprocess,datetime,gzip
import numpy as np
import time
from collections import OrderedDict

global script_path
script_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

class SLURM:

    global date_str

    now = datetime.datetime.now()
    date_str = str(now.year)+str(now.month)+str(now.day)

    def __init__(self, elements=[]):
        self.__elements={}
        for e in elements:
            self.__elements[e]=1

    def display(self):
        print "Inside SLURM class. Creating SLURM script for launching jobs in the cluster"

    def getSlurmWriteHandle(self,bin_dir,slurm_name,chr_num=[]):

        #wh = []

        if len(chr_num)!=0:
            #print "inside here "
            #print slurm_name
            slurm_file = bin_dir+"/"+slurm_name+"_chr"+chr_num+".sh"
            wh = open(slurm_file,"w")
        else:
            slurm_file = os.path.abspath(bin_dir+"/"+slurm_name+".sh")
            wh = open(slurm_file,"w")

        return slurm_file, wh

    def writeSlurmTop(self,config_dict,wh):

        print >>wh,"#!/bin/bash"
        print >>wh,"#!"

        gen_dict = config_dict["general"]
        print >>wh,"## author: "+gen_dict["author"]
        print >>wh,"## date: "+gen_dict["date"]
        print >>wh,"##\n\n"

        return wh

    def writeSlurmInit(self,config_dict,sb_log,exp_type,wh):
        
        tmp_flag = 0
        slurm_dict = config_dict["slurm"]
        script_out = sb_log+"/"+exp_type+".o.txt"
        script_err = sb_log+"/"+exp_type+".e.txt"
        
        for slr in slurm_dict:
            dash = slurm_dict["dash"]
            doubDash = slurm_dict["doubDash"]
            if slr=="params1":

                par1_dict = slurm_dict[slr]
                for par1 in par1_dict:
                    if par1=="J":
                        print >>wh,"#SBATCH "+dash+par1+" "+par1_dict[par1]+"."+exp_type
                    else:
                        print >>wh,"#SBATCH "+dash+par1+" "+par1_dict[par1]

            if slr=="params2":
                par2_dict = slurm_dict[slr]
                for par2 in par2_dict:
                    try:
                        if par2=="output":
                            print >>wh,"#SBATCH "+doubDash+par2+"="+script_out#sb_log+"/"+"gvcfGT_"+str(chr_num)+".o.txt"
                        elif par2=="error":
                            print >>wh,"#SBATCH "+doubDash+par2+"="+script_err#sb_log+"/"+"gvcfGT_"+str(chr_num)+".e.txt"
                        else:
                            print >>wh,"#SBATCH "+doubDash+par2+par2_dict[par2]

                    except:
                        print >>wh,"#SBATCH "+doubDash+par2

                print >>wh,"\n"
         
        return script_out, script_err, wh

    def writeSlurmSpecific(self,wh):
        print >>wh,"#! Number of nodes and tasks per node allocated by SLURM (do not change):"
        print >>wh,"numnodes=$SLURM_JOB_NUM_NODES"
        print >>wh,"numtasks=$SLURM_NTASKS"
        print >>wh,"mpi_tasks_per_node=$(echo \"$SLURM_TASKS_PER_NODE\" | sed -e  's/^\([0-9][0-9]*\).*$/\\1/')"
        print >>wh,"\n"
        
        return wh

    def writeSlurmModule(self,config_dict,wh):

        mod_list = config_dict["module"]["value"]
        print >>wh,". /etc/profile.d/modules.sh"
        for modKey in mod_list:
            print >>wh,"module "+modKey
        
        print >>wh,"\n\n"
        print >>wh,"JOBID=$SLURM_JOB_ID"
        print >>wh,"echo -e \"JobID: $JOBID\n======\""
        print >>wh,"echo \"Time: `date`\""
        print >>wh,"echo \"Running on master node: `hostname`\""
        print >>wh,"echo \"Current directory: `pwd`\""
        print >>wh,"echo -e \"numtasks=$numtasks, numnodes=$numnodes, mpi_tasks_per_node=$mpi_tasks_per_node (OMP_NUM_THREADS=$OMP_NUM_THREADS)\""
        print >>wh,"\n\n"

        return wh

    def getManifestDict(self,manifest):

        famDict = OrderedDict()

        fh = open(manifest)

        for lines in fh:
            lines = lines.strip()
            if re.search("^fam",lines):
                pass
            else:
                strs = re.split("\t",lines);strs = [x.strip() for x in strs]
                
                fam_id = strs[0]
                ngc_id = strs[1]
                father_id = strs[2]
                mother_id = strs[3]
                gender_dec = strs[4]
                affected = strs[5]
                sample_id = strs[6]
                fam_size = strs[7]
                load_hpc_date = strs[8]
                bam_path = strs[9]
                sample_path = strs[10]
                sample_dir = os.path.dirname(sample_path)
                #sample_sv_name = sample_dir+"/"+sample_id+".SV.vcf.gz"
                #sample_sv_name = re.split('.genome.vcf.gz',sample_path)[0]+'.SV.vcf.gz'
                sample_sv_name = strs[10] #re.split('.genome.vcf.gz',sample_path)[0]+'.SV.vcf.gz'
                
                if famDict.has_key(fam_id):
                   tmp_id = famDict[fam_id]['sample']
                   tmp_id.append(sample_id)

                   tmp_ngc = famDict[fam_id]['ngc_id']
                   tmp_ngc.append(ngc_id)
                
                   tmp_f = famDict[fam_id]['father_id']
                   tmp_f.append(father_id)
                    
                   tmp_m = famDict[fam_id]['mother_id']
                   tmp_m.append(mother_id)

                   tmp_g = famDict[fam_id]['gender']
                   tmp_g.append(gender_dec)

                   tmp_a = famDict[fam_id]['affected']
                   tmp_a.append(affected)

                   tmp_fz = famDict[fam_id]['fam_size']
                   tmp_fz.append(fam_size)

                   tmp_lhd = famDict[fam_id]['hpc_date']
                   tmp_lhd.append(load_hpc_date)

                   tmp_bp = famDict[fam_id]['bam']
                   tmp_bp.append(bam_path)

                   tmp_sv = famDict[fam_id]['spath']
                   tmp_sv.append(sample_sv_name)

                   famDict[fam_id] = {'ngc_id':tmp_ngc,'sample':tmp_id,
                                      'father_id':tmp_f,'mother_id':tmp_m, 
                                      'gender':tmp_g,'affected':tmp_a,
                                      'fam_size':tmp_fz,'hpc_date':tmp_lhd,
                                      'bam':tmp_bp,'spath':tmp_sv}
                else:

                    famDict[fam_id] = {'ngc_id':[ngc_id],'sample':[sample_id],
                                       'father_id':[father_id],'mother_id':[mother_id],
                                       'gender':[gender_dec],'affected':[affected],
                                       'fam_size':[fam_size],'hpc_date':[load_hpc_date],
                                       'bam':[bam_path],'spath':[sample_sv_name]}

        fh.close()

        return famDict

    def getEditVCF(self,inp_svFile,sv_path,ngc_id):

        inp_svFileName = re.split("\.SV.vcf.gz",os.path.basename(inp_svFile))[0]
        out_dir = sv_path
        out_svEditFile = out_dir+"/"+inp_svFileName+"."+ngc_id+".SV.edit.vcf.gz"

        fh = gzip.open(inp_svFile)
        wh = gzip.open(out_svEditFile,"wb")

        for lines in fh:
            lines = lines.strip()
            if re.search("^#",lines):
                print >>wh,lines
            elif re.search("BND",lines):
                strs = re.split("\t",lines);strs = [x.strip() for x in strs]

                chr_num = strs[0].strip()
                pos = int(strs[1].strip())
                sv_id = strs[2].strip()
                ref = strs[3].strip()

                alt_id = strs[4].strip()
                alt_id = "<INS:"+alt_id+">"

                info_strs = re.split("\;",strs[7].strip())
                sv_len = "END="+str(pos+1)
                info_strs.insert(0,sv_len)
                #print strs[4]
                #print info_strs
                qual = strs[5]
                filter_tag = strs[6]
                gt_tag = strs[8].strip()
                gt_val = strs[9].strip()
                out_line = [chr_num,str(pos),sv_id,ref,alt_id,qual,filter_tag,
                            ";".join(info_strs),gt_tag,gt_val]
                print >>wh,"\t".join(out_line)
            else:
                print >>wh,lines

        wh.close()
        fh.close()

        return out_svEditFile


    def writeSlurmVcfanno(self,manifest,configDict,famDict,out_dir,fam_id,ngc_id,
                          svcf,ilm_id,tmp_stat_file,script,genome_ref,wh):

        global script_path 
        script_path = script

        db_list = configDict["annoDB"]["db"]["value"]
       
        repos_path = configDict['general']['repos']
        vcfannoDict = configDict['vcfanno']
        vcfannoBin = vcfannoDict['binaries']
        if vcfannoDict['lua']:
            lua = script_path+"/"+vcfannoDict['lua']

        if re.search("37",genome_ref):
            toml = vcfannoDict['toml37']
        elif re.search("38",genome_ref):
            toml = vcfannoDict['toml38']

        params = vcfannoDict['params']['value']
        params = " ".join(["-"+str(x) for x in params])
        cores = str(vcfannoDict['p'])

        famSVFile = svcf
        famSVName = ilm_id
        
        inp_svFile = famSVFile
        #out_edit_svFile_gz = out_dir+"/"+famSVName+"."+ngc_id+".SV.edit.vcf.gz"
        out_edit_svFile_gz = out_dir+"/"+ngc_id+".SV.edit.vcf.gz"
        out_svFile = out_dir+"/"+ngc_id+".SV.edit.annot.vcf"

        #out_edit_svFile_gz = self.getEditVCF(inp_svFile,out_dir,ngc_id)
        #inp_svFile = out_edit_svFile_gz
        inp_svFile = famSVFile
       
        
        #Extracting information related to Editing Raw VCFfiles for processing Ins & BND 
        # points
        bndBin = script_path+"/"+configDict['processBND']['binaries']
        inpBndFile = inp_svFile
        outBndFile = out_edit_svFile_gz
        edit_cmd = ["python",bndBin,"-i",inpBndFile,"-o",outBndFile]
        
        print >>wh,"################# Begining of Editing Raw Input VCF File for BreakEnd Points ####################\n"
        print >>wh,"echo \"Launching processBND step\"\n"

        print >>wh," ".join(edit_cmd)
        wh = self.checkErrorFile(outBndFile,wh,"error",tmp_stat_file)
        print >>wh,"echo \"End of processBND step\"\n"
        print >>wh,"################# End of Editing Raw Input VCF File for BreakEnd Points ####################\n"
        
        print >>wh,"################ Begin of Normalizing multi allelic sites ##################################\n"
        
        print >>wh,'echo \"Begin of Normalizing Input VCF file\"'
        inp_svFile_gz = outBndFile
        #inp_svFile_gz = svcf
        inp_svFile_gz_name = os.path.basename(inp_svFile_gz)
        out_edit_norm_svFile_gz = out_dir+'/'+\
                                  re.split('.vcf.gz',inp_svFile_gz_name)[0]+\
                                  '.norm.vcf.gz' 
        norm_cmd = ['bcftools norm -m - ',inp_svFile_gz,' | bgzip -c > ',
                        out_edit_norm_svFile_gz
                   ]
        norm_index_cmd = 'tabix -p vcf '+out_edit_norm_svFile_gz
        print >>wh,' '.join(norm_cmd)
        print >>wh,norm_index_cmd
        wh = self.checkErrorFile(out_edit_norm_svFile_gz,wh,'error',tmp_stat_file)
        print >>wh,'echo \"End of normalizing input vcf file\"'
        print >>wh,"################ End of Normalizing multi allelic sites   ##################################\n"

        #tmp_cmd = [vcfannoBin,params,"-p",cores,"-lua",lua,toml,inp_svFile,">",out_svFile]
        inp_svFile_gz = out_edit_norm_svFile_gz#out_edit_svFile_gz
        out_edit_norm_annot_svFile_gz = out_dir+"/"+ngc_id+".SV.edit.norm.annot.vcf.gz"
        
        tmp_cmd = [vcfannoBin,params,"-p",cores,toml,inp_svFile_gz,"| bgzip -c >",
                                                    out_edit_norm_annot_svFile_gz]


        print >>wh,"################# Begining of Annotation using vcfanno ####################\n"
        print >>wh,"echo \"Launching vcfanno step\"\n"
        #print tmp_cmd
        print >>wh," ".join(tmp_cmd)

        wh = self.checkErrorFile(out_edit_norm_annot_svFile_gz,wh,"error",
                                                            tmp_stat_file)
        
        print >>wh,"tabix -p vcf "+out_edit_norm_annot_svFile_gz
        #print >>wh,"rm "+out_edit_norm_annot_svFile
        print >>wh,"\n"
        print >>wh,"echo \"End of vcfanno step\"\n"

        print >>wh,"################# End of Annotation using vcfanno ####################\n"

        ''' To do: Add multi-allelic site Normalization steps during internal 
            family overlaps\n 
        '''

        if 'fam' in db_list:
            print >>wh,'################# Begin Family Internal Overlap   ####################\n'

            famBin = script_path+'/'+configDict['vcfanno']['famOverlap']['binaries']
            winLen = configDict['vcfanno']['famOverlap']['w']

            in_svFamWinFile_gz = out_dir+"/"+ngc_id+".SV.edit.fam.vcf.gz"
            out_svFamWinFile = out_dir+"/"+ngc_id+".SV.edit.annot.fam.vcf"
            out_svFamWinFile_gz = out_svFamWinFile+".gz"
            
            in_svFamDBFile_toml = out_dir+'/famAnnoDB/'+fam_id+".fam.toml"
        
            fam_cmd = ["python",famBin,'-m',manifest,'-f',fam_id,'-o',out_dir,
                       '-w',str(winLen)] 

            print >>wh, " ".join(fam_cmd)
            wh = self.checkErrorFile(in_svFamWinFile_gz,wh,'error',tmp_stat_file)
            wh = self.checkErrorFile(in_svFamDBFile_toml,wh,'error',tmp_stat_file)

            print >>wh,'echo \"Starting VCFanno step on internal window overlap\"'
            tmp_cmd = [vcfannoBin,"-permissive-overlap","-p",cores,
                       in_svFamDBFile_toml,in_svFamWinFile_gz,">",out_svFamWinFile]
            print >>wh," ".join(tmp_cmd)

            wh = self.checkErrorFile(out_svFamWinFile,wh,"error",tmp_stat_file)
            print >>wh,"bgzip -c "+out_svFamWinFile+" > "+out_svFamWinFile_gz
            print >>wh,"tabix -p vcf "+out_svFamWinFile_gz
            print >>wh,"rm "+out_svFamWinFile
            print >>wh,"\n"
            print >>wh,"echo \"End of vcfanno step\"\n"

            print >>wh,"################# End of Family Internal Overlap annotation"

        return out_edit_norm_annot_svFile_gz,wh

        
    def writeSlurmExtractFields(self,configDict,manifest,ngc_id,fam_id,
                                inp_svFile,tmp_stat_file,genome_ref,wh):

        ''' Subroutine to extract vcfanno related fields from the annotated 
            VCF file
        '''
        db_list = configDict["annoDB"]["db"]["value"]

        for i,db in enumerate(db_list):
            if db == "ngc" and re.search("38",genome_ref):
                db_list[i] = "ngc38"
                db_list.append("ngclov")
            if db == "ngc" and re.search("37",genome_ref):
                db_list[i] = "ngc37"

        db_list = [str(x.strip()) for x in db_list]
        db_list.sort()
        db_list.insert(0,"rest")

        out_dir = os.path.dirname(inp_svFile)
        out_annotDB = out_dir+"/annotDB"
        out_overlap = out_dir+"/overlap"

        os.system("mkdir -p "+out_annotDB)
        os.system("mkdir -p "+out_overlap)

        coord_file = []

        for db in db_list:

            print >>wh,"echo \"Extracting Relevant fields for downstream analysis: "+\
                        db," \" \n"

            baseFile = os.path.basename(inp_svFile)
            inp_svFile_orig = inp_svFile

            if db != "fam":
                out_svFile_tsv = out_annotDB+"/"+re.split(".vcf.gz",baseFile)[0]+\
                                                            "."+db.lower()+".bed"
                out_svFile_tsv_gz = out_svFile_tsv+".gz"

            elif db == "fam":
                inp_svFile_orig = inp_svFile
                inp_svFile = re.split(".vcf.gz",inp_svFile)[0]+".fam.vcf.gz"
                out_svFile_tsv = out_annotDB+"/"+re.split(".vcf.gz",baseFile)[0]+\
                                                                       ".fam.bed" 
                out_svFile_tsv_gz = out_svFile_tsv+".gz"

            if db == "rest":
                if re.search("38",genome_ref):
                    cmd = "bcftools query -f \"%CHROM+%POS+%INFO/END+%INFO/SVTYPE+" \
                          "%ID+[%GT]\\t%REF\\t%ALT\\t%QUAL\\t%FILTER\\t%INFO/CIPOS"\
                          "\\t%INFO/CIEND\\t%INFO/SVLEN\\t%INFO/CNVLEN\\t[%GT:%FT]"\
                          "\\n\" "+inp_svFile+" | bgzip -c > "+out_svFile_tsv_gz
                
                elif re.search("37",genome_ref):
                    cmd = "bcftools query -f \"%CHROM+%POS+%INFO/END+%INFO/SVTYPE+"\
                          "%ID+[%GT]\\t%REF\\t%ALT\\t%QUAL\\t%FILTER\\t%INFO/CIPOS"\
                          "\\t%INFO/CIEND\\t%INFO/SVLEN\\t[%GT:%FT]"\
                          "\\n\" "+inp_svFile+" | bgzip -c > "+out_svFile_tsv_gz
                coord_file = out_svFile_tsv+".gz"

            elif db == "fam":
                cmd = "bcftools query -f \"%CHROM+%POS+%INFO/END+%INFO/SVTYPE+"\
                      "%ID+[%GT]\\t%INFO/CIPOS\\t%INFO/CIEND\\t%INFO/"+db+\
                       "_s\\t%INFO/"+db+"_e\\t%INFO/"+db+"_ov_id\\n\" "+inp_svFile+\
                       " | bgzip -c > "+out_svFile_tsv_gz
                inp_svFile = inp_svFile_orig

            elif db in ["dbvar","gnomad","ngc38","ngclov","decipher","ngc37",
                                 'ens_trans','ens_exon','ref_gene','ref_exon',
                                                                       'user']:

                cmd = "bcftools query -f \"%CHROM+%POS+%INFO/END+%INFO/SVTYPE+"\
                      "%ID+[%GT]\\t%INFO/CIPOS\\t%INFO/CIEND\\t%INFO/"+db+\
                       "_s\\t%INFO/"+db+"_e\\t%INFO/"+db+"_ov_id\\t%INFO/left_"+\
                       db+"_s\\t%INFO/left_"+db+"_e\\t%INFO/left_"+db+\
                       "_ov_id\\t%INFO/right_"+db+"_s\\t%INFO/right_"+db+\
                       "_e\\t%INFO/right_"+db+"_ov_id\\n\" " +inp_svFile+\
                       "| bgzip -c > "+out_svFile_tsv_gz
            
            elif db =="ensembl":
                
                ensDict = configDict['ensembl']['ensVarAnnot']

                coordBin = script_path+'/'+ensDict['ensCoordBin']
                grBin    = script_path+'/'+ensDict['ensGrangeBin']
                inp_coord = coord_file
                out_coord = re.split(".rest.bed.gz",inp_coord)[0]+".ensembl.bed"

                inp_gr_file = out_coord
                out_gr_file =  out_dir+"/overlap/"+ngc_id+".ensembl.overlap.same.bed"
                
                gene_file = ensDict['g']

                if re.search("37",genome_ref):
                    exon_file = ensDict['x37']
                    trans_file = ensDict['t37']

                elif re.search("38",genome_ref):
                    exon_file = ensDict['x38']
                    trans_file = ensDict['t38']
                    
                manifest = manifest
                pheno_file = ensDict['n']
                hpodb_file = ensDict['d']

                coord_cmd = " ".join(["python",coordBin,inp_coord, out_coord])
                cmd       = " ".join(["Rscript",grBin,"-i",inp_gr_file,"-o",out_gr_file,
                                      "-c",ngc_id,"-g",gene_file,"-x",exon_file,
                                      "-t",trans_file,"-n",pheno_file,"-d",hpodb_file,
                                      "-p",manifest])
                
                print >>wh,coord_cmd,"\n"
                wh = self.checkErrorFile(out_svFile_tsv,wh,"error",tmp_stat_file)
                out_svFile_tsv = out_gr_file
                #print >>wh,"bgzip -c "+out_svFile_tsv+" > "+out_svFile_tsv+".gz"
                
            #print >>wh,cmd
            #print >>wh,"bgzip -c "+out_svFile_tsv+" > "+out_svFile_tsv+".gz"
            #print >>wh,"tabix -p bed "+out_svFile_tsv+".gz"
            #print >>wh,"rm "+ out_svFile_tsv
            if db != "ensembl": 
                print >>wh,cmd
                wh = self.checkErrorFile(out_svFile_tsv+".gz",wh,"error",tmp_stat_file)
            else:
                print >>wh,cmd
                print >>wh,"bgzip -c "+out_svFile_tsv+" > "+out_svFile_tsv+".gz"
                wh = self.checkErrorFile(out_svFile_tsv+".gz",wh,"error",tmp_stat_file)
                print >>wh,"rm "+ out_svFile_tsv
            print >>wh,"\n"
        
        db_list = []
        
        print >>wh,'echo "Removing the annotation file\n"'
        print >>wh,'rm -rf '+inp_svFile

        print >>wh,"############### End of reformatting and indexing VCF ####################\n"
        
        return out_annotDB, wh

    def writeSlurmOverlapMerge(self,configDict,manifest,ngc_id,fam_id,xml_file,
                                      out_annotdb,tmp_stat_file,genome_ref,wh):

        db_list = configDict["annoDB"]["db"]["value"]
        db_list = [str(x.strip()) for x in db_list]

        for i,db in enumerate(db_list):
            if db == "ngc" and re.search("38",genome_ref):
                db_list[i] = "ngc38"
                db_list.append("ngclov")
            if db == "ngc" and re.search("37",genome_ref):
                db_list[i] = "ngc37"

        db_list.sort()
        db_list.insert(0,"rest")

        overlapBin = script_path+'/'+configDict['overlapMerge']['ovBin']
        mergeBin = script_path+'/'+configDict['overlapMerge']['mergeBin']
        ovFrac_list = re.split("\,",configDict['overlapMerge']['ovFrac'])

        out_overlap = os.path.dirname(out_annotdb)+"/overlap"
        out_annotdb = os.path.dirname(out_annotdb)+"/annotDB"

        #print annot_files
        rest_file = out_annotdb+"/"+ngc_id+".SV.edit.norm.annot.rest.bed.gz"

        print >>wh, "############# Begin of Overlap and Formatting Step ###############\n"

        out_ovp_list = []


        for db_type in db_list:
            ovFrac_list = re.split("\,",configDict['overlapMerge']['ovFrac'])
            
            if not db_type in ["ensembl","rest"]:
                
                print >>wh,"echo \"Start of Overlap and formatting step for "+\
                                                        "db: "+db_type+"\"\n"
                
                if not db_type in ['gnomad','ngc38','ngclov','ngc37']:
                    
                    ovFrac = '0.0'
                    print >>wh,"echo \"Overlap fraction used: "+str(ovFrac)+"\"\n" 
                    ele_file = out_annotdb+"/"+ngc_id+".SV.edit.norm.annot."+\
                                                db_type.lower()+".bed.gz"
                    overlap_cmd = ["python",overlapBin,"-i",ele_file,"-f",ngc_id,
                                            "-t","overlap","-a",db_type.lower(),
                                                      "-m",manifest,"-r",ovFrac,
                                                                  '-x',xml_file
                                  ]
                    print >>wh,' '.join(overlap_cmd)
                    out_svFile_tsv_gz = out_overlap+"/"+".".join([ngc_id,
                                                         db_type.lower(),
                                                            str(ovFrac)+\
                                                 ".overlap.same.bed.gz"])
                    wh = self.checkErrorFile(out_svFile_tsv_gz,wh,"error",tmp_stat_file)
                    #print >>wh,"echo \"End of Overlap and formatting step for db: "+\
                                                                          #db_type+"\"\n"
                    
                else:
                    cmd_list = []
                    out_tsv_gz_list = []
                    
                    # Added this step for computing NGC family internal overlap 
                    # with range of overlap fractions [0.0,1.0]
                    #
                    # For gnomad and ngclov-37 default is 0.70 overlap fraction

                    if db_type in ['ngc38','ngc37']:
                        ovFrac_list = ovFrac_list

                    elif db_type in ['gnomad','ngclov']:
                        ovFrac_list = ['0.7']

                    for ovFrac in ovFrac_list:
               
                        ele_file = out_annotdb+"/"+ngc_id+".SV.edit.norm.annot."+\
                                                            db_type.lower()+".bed.gz"
                        overlap_cmd = ["python",overlapBin,"-i",ele_file,"-f",ngc_id,
                                                "-t","overlap","-a",db_type.lower(),
                                                          "-m",manifest,"-r",ovFrac,
                                                                      '-x',xml_file
                                      ]
                        #print >>wh," ".join(overlap_cmd)
                        out_svFile_tsv_gz = out_overlap+"/"+".".join([ngc_id,
                                                             db_type.lower(),
                                                                str(ovFrac)+\
                                                     ".overlap.same.bed.gz"])
                        
                        cmd_list += ['echo \"'+' '.join(overlap_cmd)+'\"']
                        out_tsv_gz_list.append(out_svFile_tsv_gz)
                    
                    print >>wh,"echo \"Overlap fraction used: "+','.join(ovFrac_list)+"\"\n" 
                    print >>wh,'('+' ; '.join(cmd_list)+') | parallel -j15 --no-notice '
                    
                    for out_svFile_tsv_gz in out_tsv_gz_list:
                        wh = self.checkErrorFile(out_svFile_tsv_gz,wh,"error",tmp_stat_file)

                print >>wh,"echo \"End of Overlap and formatting step for db: "+\
                                                                                db_type+"\"\n"
                
        print >>wh, "############## End of Overlap and Formatting Step ###############\n"
            
        print >>wh, "\n############ Begin of Merge and Formatting Step ###############\n"
        print >>wh, "echo \"Start of Merge and Formatting step\"\n" 
           
        for ovFrac in ovFrac_list:
            merge_cmd = ["python",mergeBin,'-i',rest_file,'-m',manifest,
                            '-a',out_overlap,'-f',ngc_id,'-r',genome_ref, 
                                                             '-o',ovFrac
                        ] 
            print >>wh," ".join(merge_cmd)
        
            out_svFile_tsv = os.path.dirname(out_overlap)+"/"+ngc_id+\
                                        ".merged.all."+str(ovFrac)+".overlap.bed.gz"

        
            out_ovp_list.append(out_svFile_tsv)
            wh = self.checkErrorFile(out_svFile_tsv,wh,"error",tmp_stat_file)
        
        print >>wh, "echo \"End of Merge and Formatting step\"\n" 

        print >>wh, "\n############## End of Merge and Formatting Step ###############\n"

        return out_ovp_list, wh

    def writeSlurmFormatGT(self,configDict,manifest,fam_id,out_ovp_mrg_list,
                                                            tmp_stat_file,wh):
        ''' Subroutine to launch formating genotype script '''

        print >>wh, 'echo \" Begin of formatting of genotype of merged step \"'
        fmt_gt_script = script_path+'/'+configDict['formatGT']['scripts']
        refGenome = configDict['general']['genomeBuild']
        par_file = configDict['par'][refGenome]
        out_ovp_mrg_fmt_list = []


        print >>wh, "\n############## Begin of Merge and Formatting Step ###############\n"

        for inpFile in out_ovp_mrg_list:
            outFile = re.split('.bed.gz',inpFile)[0]+'.fmt.bed.gz'
            cmd = ['python',fmt_gt_script,'-m',manifest,'-i',inpFile,'-o',outFile,
                                                        '-f',fam_id,'-p',par_file,
                                                                    '-r',refGenome
                  ]
            print >>wh,' '.join(cmd)
            out_ovp_mrg_fmt_list.append(outFile)
            wh = self.checkErrorFile(outFile,wh,"error",tmp_stat_file)
        
        print >>wh, 'echo \" End of formatting of genotype of merged step \"'
        print >>wh, "\n############## End of Merge and Formatting Step ###############\n"

        return out_ovp_mrg_fmt_list,wh

    def writeSlurmFamFilter(self,configDict,manifest,ngc_id,fam_id,out_ovp_mrg_list,
                                                fam_ngc_dir,tmp_stat_file,wh):

        ''' Subroutine Function to launch Rscript family filtering. '''

        print >>wh,"########### Starting Family filtering for family: ",\
                    fam_id," ############\n"
        print >>wh, "\necho \"Begin of Family-Filtering step\"\n" 
       
        script_path   = os.path.dirname(os.path.dirname(os.path.abspath( __file__ )))
        famBin        = script_path+'/'+configDict['famFilter']['famBin']
        ovFrac_list   = re.split('\,',configDict['overlapMerge']['ovFrac'])
        imprint_genes = configDict['imprint']
        ref_genome    = configDict['general']['genomeBuild']

        filter_dir    = fam_ngc_dir+'/filter/'
        os.system('mkdir -p '+filter_dir)
        
        for inpFile in out_ovp_mrg_list:
            for ovFrac in ovFrac_list:
                if re.search(str(ovFrac)+'.overlap',inpFile):
                    fam_cmd = ['Rscript ',famBin,'-v',inpFile,'-o',filter_dir,'-f',fam_id,
                                                    '-n',ngc_id,'-p',manifest,'-c',ovFrac,
                                                    '-i',imprint_genes,'-r',ref_genome
                              ]

                    print >>wh,' '.join(fam_cmd)
                    
        
        print >>wh, "echo \"End of Family-Filtering step\"\n" 
        
        print >>wh,"########### End of Family filtering for family: ",fam_id,\
                    " ############\n"

        return wh


    def writeSlurmAnnotSV(self,configDict,famDict,out_dir,fam_id,ngc_id,
                                                    tmp_stat_file,k,wh):
        
        #famDict = OrderedDict()

        #Extract the relevant binaries and parameter for AnnotSV module
        svDict = configDict["annotSV"]
        svBin = svDict["binaries"]
        gbuilt = svDict["genomeBuild"]
        bedT = svDict["bedtools"]
        svCol = svDict["svtBEDCol"]
        typeAnno = svDict["typeOfAnnotation"]
        svInpInfo = str(svDict["SVinputInfo"])
        reciprocalValue = str(svDict["reciprocal"])
        overlapPercentage = str(svDict["overlap"])
       
        extGnDict =  configDict["extractGnomad"]
        extGnBin = extGnDict["binaries"]

        procGnDict =  configDict["processGnomad"]
        procGnBin = procGnDict["binaries"]

        sv_path = os.path.abspath(out_dir+"/sv")

        #fh = open(fam_file)
        cmd = []
        fam_ngc = {}

        chr_list = [str(x) for x in range(1,23)]+['X','Y','MT']
        # Creating AnnotSV command foreach of the family memebrs of given family IDs
        fam_sv_dir = os.path.abspath(sv_path+"/"+fam_id)
        ngc_list = [ngc_id]

        for i in range(0,len(ngc_list)):
            fam_ngc_dir = os.path.abspath(sv_path+"/"+fam_id+"/"+ngc_list[i])
            #os.system("mkdir -p "+fam_ngc_dir)
            ngc_id = ngc_list[i] #famDict[fam_id]['ngc_id'][i]
            famSVFile = famDict[fam_id]['spath'][k]
            famSVName = famDict[fam_id]['sample'][k]
            
            for chr_num in chr_list:
                fam_ngc_chr_dir = os.path.abspath(fam_ngc_dir+"/chr"+str(chr_num))
                os.system("mkdir -p "+fam_ngc_chr_dir)
                print >>wh,"\n"
                print >>wh,"############## Begin Analysis of Family Member: "+\
                            ngc_list[i]+" #################\n"

                out_chr_famSVName = famSVName+".genome.SV.chr"+str(chr_num)
                print >>wh, "bcftools view -r chr"+str(chr_num)+" "+famSVFile+\
                            " -Oz -o "+fam_ngc_chr_dir+"/"+out_chr_famSVName+".vcf.gz"
                print >>wh, "bcftools index"+" "+fam_ngc_chr_dir+"/"+\
                            out_chr_famSVName+".vcf.gz"
                wh = self.checkErrorFile(fam_ngc_chr_dir+"/"+out_chr_famSVName+
                                         ".vcf.gz",wh,"error",tmp_stat_file)

                print >>wh,"a="+"`bcftools view "+fam_ngc_chr_dir+"/"+\
                            out_chr_famSVName+".vcf.gz | grep -v '^#'| wc -l`"
                print >>wh,"if [ $a != 0 ]; then"
                out_chr_annotSVName = ngc_id+".chr"+str(chr_num)+".AnnotSV.output"
                tmp_cmd = svBin+" -SVinputFile "+fam_ngc_chr_dir+"/"+\
                          out_chr_famSVName+".vcf.gz"+" -SVinputInfo "+svInpInfo+\
                        " -outputDir "+fam_ngc_chr_dir+" -outputFile "+\
                        out_chr_annotSVName+" -svtBEDCol "+str(svCol)+\
                        " -typeOfAnnotation "+typeAnno+" -genomeBuild "+\
                        gbuilt+" -bedtools "+bedT+" -reciprocal "+reciprocalValue+\
                        " -overlap "+overlapPercentage
                
                print >>wh,tmp_cmd
                wh = self.checkErrorFile(fam_ngc_chr_dir+"/"+out_chr_annotSVName+
                                         ".tsv",wh,"error",tmp_stat_file)
                
                print >>wh,"fi\n" 
               
            print >>wh,"echo \"End of analyzing NGC: "+str(ngc_list[i])+\
                        " for family: "+fam_id+"\""
            print >>wh,"############## End Analysis of Family Member: "+\
                        ngc_list[i]+" #################\n\n"


        return fam_ngc_dir,ngc_list, wh    


    def writeSlurmCombChrOut(self,configDict,fam_ngc_dir,ngc_id,tmp_stat_file,wh):
        
        combDict = configDict["combChromOut"]
        mergeBin = script_path+'/'+combDict["binaries"]

        
        out_dir = inp_dir = fam_ngc_dir

        out_merge_file = ngc_id+".merged.all.output"
        comb_cmd = "Rscript "+mergeBin+" -i "+inp_dir+" -o "+out_dir+" -s "+out_merge_file+".tsv"
        
        print >>wh,"\n############## Begin Merge AnnotSV output for Family id : "+ngc_id+" #####################\n"
        print >>wh,comb_cmd

        wh = self.checkErrorFile(fam_ngc_dir+"/"+out_merge_file+".tsv",wh,"error",tmp_stat_file)
        print >>wh,"\n############# End of Merging of AnnotSV output for Family id: "+ngc_id+" ###############\n"

        return wh


    def writeSlurmExtractGnomad(self,configDict,sv_ext_dict,wh):

        extGnDict =  script_path+'/'+configDict["extractGnomad"]
        extGnBin = script_path+'/'+extGnDict["binaries"]
        
        cmd = []

        sv_gn_dict = {}

        for ngc_id in sv_ext_dict:
            gnomadOutStr = sv_ext_dict[ngc_id]['outStr']+".gnomad"
            ngc_path = sv_ext_dict[ngc_id]["path"]
            extGnOutFile = gnomadOutStr+".tsv"

            sv_gn_dict[ngc_id] = {'path':ngc_path,'file':extGnOutFile,'outStr':gnomadOutStr}
            tmp_cmd = " echo \""+extGnBin+" -i "+sv_ext_dict[ngc_id]["file"]+" -o "+ngc_path+" -s "+extGnOutFile#+"\""
            cmd.append(tmp_cmd)

        print >>wh,"("+";".join(cmd)+")|parallel"

        return sv_gn_dict, wh

    def writeSlurmProcessGnomad(self,configDict,sv_gn_dict,wh):
        
        procGnDict =  configDict["processGnomad"]
        procGnBin = procGnDict["binaries"]
        
        cmd = []

        for ngc_id in sv_gn_dict:
            filterOutStr = sv_gn_dict[ngc_id]['outStr']+".filter"
            ngc_path = sv_gn_dict[ngc_id]['path']
            extGnInFile = sv_gn_dict[ngc_id]['file']
            tmp_cmd = " echo \""+procGnBin+" -i "+extGnInFile+" -o "+ngc_path+" -s "+filterOutStr+"\""
            cmd.append(tmp_cmd)

        print >>wh,"("+";".join(cmd)+")| parallel"

        return wh

    def checkErrorFile(self,out_file,wh,flag,status_file=[]):
        
        if_cmd = "\nif [ ! -e "+out_file+" ]; then \n"
        if flag=="error":
            echo_1 = "   echo -e \"ERROR: "+out_file+" doesnot exist\" \n"
            echo_2 = "   echo \"0\" | cat > \""+status_file+"\"\n"
            echo_3 = "   exit 1\n"
            echo_4 = "   else\n"
            echo_5 = "          echo \"1\" | cat > \""+status_file+"\""
        end_cmd = "\nfi"
        print >>wh,if_cmd+echo_1+echo_2+echo_3+echo_4+echo_5+end_cmd
        print >>wh,"\n"
        return wh

    def writeSlurmEndStatus(self,tmp_stat_file,wh):

        print >>wh,"echo \"2\" | cat > "+tmp_stat_file

        return wh

class DEMO:
    ''' Class containing subroutines & methods for generating annoated SVs'''

    def __init__(self, elements=[]):
        self.__elements={}
        for e in elements:
            self.__elements[e]=1


    def writeDemoVcfanno(self,manifest,configDict,famDict,out_dir,fam_id,ngc_id,
                          svcf,ilm_id,tmp_stat_file,script,genome_ref,wh):

        ''' Subroutine for launching vcfanno based annotation '''

        global script_path 
        script_path = script

        db_list = configDict["annoDB"]["db"]["value"]
       
        vcfannoDict = configDict['vcfanno']
        vcfannoBin = vcfannoDict['binaries']
        if vcfannoDict['lua']:
            lua = script_path+"/"+vcfannoDict['lua']

        if re.search("37",genome_ref):
            toml = vcfannoDict['toml37']
        elif re.search("38",genome_ref):
            toml = vcfannoDict['toml38']

        params = vcfannoDict['params']['value']
        params = " ".join(["-"+str(x) for x in params])
        cores = str(vcfannoDict['p'])

        famSVFile = svcf
        famSVName = ilm_id
        
        inp_svFile = famSVFile
        #out_edit_svFile_gz = out_dir+"/"+ngc_id+".SV.edit.vcf.gz"
        out_svFile = out_dir+"/"+ngc_id+".SV.edit.annot.vcf"

        ''' 
        #Extracting information related to Editing Raw VCFfiles for processing Ins & BND 
        # points
        bndBin = script_path+"/"+configDict['processBND']['binaries']
        inpBndFile = inp_svFile
        outBndFile = out_edit_svFile_gz
        edit_cmd = ["python",bndBin,"-i",inpBndFile,"-o",outBndFile]
        
        print >>wh,"################# Begining of Editing Raw Input VCF File for BreakEnd Points ####################\n"
        print >>wh,"echo \"Launching processBND step\"\n"

        print >>wh," ".join(edit_cmd)
        wh = self.checkErrorFile(outBndFile,wh,"error",tmp_stat_file)
        print >>wh,"echo \"End of processBND step\"\n"
        print >>wh,"################# End of Editing Raw Input VCF File for BreakEnd Points ####################\n"
        '''
        print >>wh,"################ Begin of Normalizing multi allelic sites ##################################\n"
        
        print >>wh,'echo \"Begin of Normalizing Input VCF file\"'
        inp_svFile_gz = inp_svFile
        inp_svFile_gz_name = os.path.basename(inp_svFile_gz)
        out_edit_norm_svFile_gz = out_dir+'/'+\
                                  re.split('.vcf.gz',inp_svFile_gz_name)[0]+\
                                  '.norm.vcf.gz' 
        norm_cmd = ['bcftools norm -m - ',inp_svFile_gz,' | bgzip -c > ',
                        out_edit_norm_svFile_gz
                   ]
        norm_index_cmd = 'tabix -p vcf '+out_edit_norm_svFile_gz
        print >>wh,' '.join(norm_cmd)
        print >>wh,norm_index_cmd
        
        wh = self.checkErrorFile(out_edit_norm_svFile_gz,wh,'error',tmp_stat_file)
        print >>wh,'echo \"End of normalizing input vcf file\"'
        print >>wh,"################ End of Normalizing multi allelic sites   ##################################\n"

        inp_svFile_gz = out_edit_norm_svFile_gz#out_edit_svFile_gz
        out_edit_norm_annot_svFile_gz = out_dir+"/"+ngc_id+".SV.edit.norm.annot.vcf.gz"
        
        tmp_cmd = [vcfannoBin,params,"-p",cores,toml,inp_svFile_gz,"| bgzip -c >",
                                                    out_edit_norm_annot_svFile_gz]


        print >>wh,"################# Begining of Annotation using vcfanno ####################\n"
        print >>wh,"echo \"Launching vcfanno step\"\n"
        print >>wh," ".join(tmp_cmd)

        wh = self.checkErrorFile(out_edit_norm_annot_svFile_gz,wh,"error",
                                                            tmp_stat_file)
        
        print >>wh,"tabix -p vcf "+out_edit_norm_annot_svFile_gz
        #print >>wh,"rm "+out_edit_norm_annot_svFile
        print >>wh,"\n"
        print >>wh,"echo \"End of vcfanno step\"\n"

        print >>wh,"################# End of Annotation using vcfanno ####################\n"


        return out_edit_norm_annot_svFile_gz,wh

    def writeDemoExtractFields(self,configDict,manifest,ngc_id,fam_id,
                                inp_svFile,tmp_stat_file,genome_ref,wh):

        ''' Subroutine to extract vcfanno related fields from the annotated 
            VCF file
        '''
        db_list = configDict["annoDB"]["db"]["value"]

        for i,db in enumerate(db_list):
            if db == "ngc" and re.search("38",genome_ref):
                db_list[i] = "ngc38"
                #db_list.append("ngclov")
            if db == "ngc" and re.search("37",genome_ref):
                db_list[i] = "ngc37"

        db_list = [str(x.strip()) for x in db_list]
        db_list.sort()
        db_list.insert(0,"rest")

        out_dir = os.path.dirname(inp_svFile)
        out_annotDB = out_dir+"/annotDB"
        out_overlap = out_dir+"/overlap"

        print >>wh, "mkdir -p "+out_annotDB
        print >>wh, "mkdir -p "+out_overlap

        coord_file = []

        for db in db_list:

            print >>wh,"echo \"Extracting Relevant fields for downstream analysis: "+\
                        db," \" \n"

            baseFile = os.path.basename(inp_svFile)
            inp_svFile_orig = inp_svFile

            if db != "fam":
                out_svFile_tsv = out_annotDB+"/"+re.split(".vcf.gz",baseFile)[0]+\
                                                            "."+db.lower()+".bed"
                out_svFile_tsv_gz = out_svFile_tsv+".gz"

            elif db == "fam":
                inp_svFile_orig = inp_svFile
                inp_svFile = re.split(".vcf.gz",inp_svFile)[0]+".fam.vcf.gz"
                out_svFile_tsv = out_annotDB+"/"+re.split(".vcf.gz",baseFile)[0]+\
                                                                       ".fam.bed" 
                out_svFile_tsv_gz = out_svFile_tsv+".gz"

            if db == "rest":
                if re.search("38",genome_ref):
                    cmd = "bcftools query -f \"%CHROM+%POS+%INFO/END+%INFO/SVTYPE+" \
                          "%ID+[%GT]\\t%REF\\t%ALT\\t%QUAL\\t%FILTER\\t%INFO/CIPOS"\
                          "\\t%INFO/CIEND\\t%INFO/SVLEN\\t[%GT:%AB]"\
                          "\\n\" "+inp_svFile+" | bgzip -c > "+out_svFile_tsv_gz
                
                elif re.search("37",genome_ref):
                    cmd = "bcftools query -f \"%CHROM+%POS+%INFO/END+%INFO/SVTYPE+"\
                          "%ID+[%GT]\\t%REF\\t%ALT\\t%QUAL\\t%FILTER\\t%INFO/CIPOS"\
                          "\\t%INFO/CIEND\\t%INFO/SVLEN\\t[%GT:%AB]"\
                          "\\n\" "+inp_svFile+" | bgzip -c > "+out_svFile_tsv_gz
                coord_file = out_svFile_tsv+".gz"

            elif db in ["dbvar","gnomad","ngc38","decipher","ngc37",
                                 'ens_trans','ens_exon','ref_gene','ref_exon',
                                                        'promoter','blacklist']:

                cmd = "bcftools query -f \"%CHROM+%POS+%INFO/END+%INFO/SVTYPE+"\
                      "%ID+[%GT]\\t%INFO/CIPOS\\t%INFO/CIEND\\t%INFO/"+db+\
                       "_s\\t%INFO/"+db+"_e\\t%INFO/"+db+"_ov_id\\t%INFO/left_"+\
                       db+"_s\\t%INFO/left_"+db+"_e\\t%INFO/left_"+db+\
                       "_ov_id\\t%INFO/right_"+db+"_s\\t%INFO/right_"+db+\
                       "_e\\t%INFO/right_"+db+"_ov_id\\n\" " +inp_svFile+\
                       "| bgzip -c > "+out_svFile_tsv_gz
            
            if db != "ensembl": 
                print >>wh,cmd
                wh = self.checkErrorFile(out_svFile_tsv+".gz",wh,"error",tmp_stat_file)
            else:
                print >>wh,cmd
                print >>wh,"bgzip -c "+out_svFile_tsv+" > "+out_svFile_tsv+".gz"
                wh = self.checkErrorFile(out_svFile_tsv+".gz",wh,"error",tmp_stat_file)
                print >>wh,"rm "+ out_svFile_tsv
            print >>wh,"\n"
        
        db_list = []
        
        print >>wh,'echo "Removing the annotation file\n"'
        print >>wh,'rm -rf '+inp_svFile

        print >>wh,"############### End of reformatting and indexing VCF ####################\n"
        
        return out_annotDB, wh

    def writeDemoOverlapMerge(self,configDict,manifest,ngc_id,fam_id,xml_file,
                                      out_annotdb,tmp_stat_file,genome_ref,wh):

        db_list = configDict["annoDB"]["db"]["value"]
        db_list = [str(x.strip()) for x in db_list]

        for i,db in enumerate(db_list):
            if db == "ngc" and re.search("38",genome_ref):
                db_list[i] = "ngc38"
                #db_list.append("ngclov")
            if db == "ngc" and re.search("37",genome_ref):
                db_list[i] = "ngc37"

        db_list.sort()
        db_list.insert(0,"rest")

        overlapBin = script_path+'/'+configDict['overlapMerge']['ovBin']
        mergeBin = script_path+'/'+configDict['overlapMerge']['mergeBin']
        ovFrac_list = re.split("\,",configDict['overlapMerge']['ovFrac'])

        out_overlap = os.path.dirname(out_annotdb)+"/overlap"
        out_annotdb = os.path.dirname(out_annotdb)+"/annotDB"

        #print annot_files
        rest_file = out_annotdb+"/"+ngc_id+".SV.edit.norm.annot.rest.bed.gz"

        print >>wh, "############# Begin of Overlap and Formatting Step ###############\n"

        out_ovp_list = []


        for db_type in db_list:
            ovFrac_list = re.split("\,",configDict['overlapMerge']['ovFrac'])
            
            if not db_type in ["ensembl","rest"]:
                
                print >>wh,"echo \"Start of Overlap and formatting step for "+\
                                                        "db: "+db_type+"\"\n"
                
                if not db_type in ['gnomad','ngc38','ngclov','ngc37']:
                    
                    ovFrac = '0.0'
                    print >>wh,"echo \"Overlap fraction used: "+str(ovFrac)+"\"\n" 
                    ele_file = out_annotdb+"/"+ngc_id+".SV.edit.norm.annot."+\
                                                db_type.lower()+".bed.gz"
                    overlap_cmd = ["python",overlapBin,"-i",ele_file,"-f",ngc_id,
                                            "-t","overlap","-a",db_type.lower(),
                                                      "-m",manifest,"-r",ovFrac,
                                                                  '-x',xml_file
                                  ]
                    print >>wh,' '.join(overlap_cmd)
                    out_svFile_tsv_gz = out_overlap+"/"+".".join([ngc_id,
                                                         db_type.lower(),
                                                            str(ovFrac)+\
                                                 ".overlap.same.bed.gz"])
                    wh = self.checkErrorFile(out_svFile_tsv_gz,wh,"error",tmp_stat_file)
                    #print >>wh,"echo \"End of Overlap and formatting step for db: "+\
                                                                          #db_type+"\"\n"
                    
                else:
                    cmd_list = []
                    out_tsv_gz_list = []
                    
                    # Added this step for computing NGC family internal overlap 
                    # with range of overlap fractions [0.0,1.0]
                    #
                    # For gnomad and ngclov-37 default is 0.70 overlap fraction

                    if db_type in ['ngc38','ngc37']:
                        ovFrac_list = ovFrac_list

                    elif db_type in ['gnomad','ngclov']:
                        ovFrac_list = ['0.7']

                    for ovFrac in ovFrac_list:
               
                        ele_file = out_annotdb+"/"+ngc_id+".SV.edit.norm.annot."+\
                                                            db_type.lower()+".bed.gz"
                        overlap_cmd = ["python",overlapBin,"-i",ele_file,"-f",ngc_id,
                                                "-t","overlap","-a",db_type.lower(),
                                                          "-m",manifest,"-r",ovFrac,
                                                                      '-x',xml_file
                                      ]
                        #print >>wh," ".join(overlap_cmd)
                        out_svFile_tsv_gz = out_overlap+"/"+".".join([ngc_id,
                                                             db_type.lower(),
                                                                str(ovFrac)+\
                                                     ".overlap.same.bed.gz"])
                        
                        cmd_list += ['echo \"'+' '.join(overlap_cmd)+'\"']
                        out_tsv_gz_list.append(out_svFile_tsv_gz)
                    
                    print >>wh,"echo \"Overlap fraction used: "+','.join(ovFrac_list)+"\"\n" 
                    print >>wh,'('+' ; '.join(cmd_list)+') | parallel -j15 --no-notice '
                    
                    for out_svFile_tsv_gz in out_tsv_gz_list:
                        wh = self.checkErrorFile(out_svFile_tsv_gz,wh,"error",tmp_stat_file)

                print >>wh,"echo \"End of Overlap and formatting step for db: "+\
                                                                                db_type+"\"\n"
                
        print >>wh, "############## End of Overlap and Formatting Step ###############\n"
            
        print >>wh, "\n############ Begin of Merge and Formatting Step ###############\n"
        print >>wh, "echo \"Start of Merge and Formatting step\"\n" 
           
        for ovFrac in ovFrac_list:
            merge_cmd = ["python",mergeBin,'-i',rest_file,'-m',manifest,
                            '-a',out_overlap,'-f',ngc_id,'-r',genome_ref, 
                                                             '-o',ovFrac
                        ] 
            print >>wh," ".join(merge_cmd)
        
            out_svFile_tsv = os.path.dirname(out_overlap)+"/"+ngc_id+\
                                        ".merged.all."+str(ovFrac)+".overlap.bed.gz"

        
            out_ovp_list.append(out_svFile_tsv)
            wh = self.checkErrorFile(out_svFile_tsv,wh,"error",tmp_stat_file)
        
        print >>wh, "echo \"End of Merge and Formatting step\"\n" 

        print >>wh, "\n############## End of Merge and Formatting Step ###############\n"

        return out_ovp_list, wh

    def writeDemoFormatGT(self,configDict,manifest,fam_id,out_ovp_mrg_list,
                                                            tmp_stat_file,wh):
        ''' Subroutine to launch formating genotype script '''

        print >>wh, 'echo \" Begin of formatting of genotype of merged step \"'
        fmt_gt_script = script_path+'/'+configDict['formatGT']['scripts']
        resource_path = configDict['general']['resourceDir']
        refGenome = configDict['general']['genomeBuild']
        par_file = '/'.join([resource_path,
                             configDict['par'][refGenome]
                            ])
        out_ovp_mrg_fmt_list = []


        print >>wh, "\n############## Begin of Merge and Formatting Step ###############\n"

        for inpFile in out_ovp_mrg_list:
            outFile = re.split('.bed.gz',inpFile)[0]+'.fmt.bed.gz'
            cmd = ['python',fmt_gt_script,'-m',manifest,'-i',inpFile,'-o',outFile,
                                                        '-f',fam_id,'-p',par_file,
                                                                    '-r',refGenome
                  ]
            print >>wh,' '.join(cmd)
            out_ovp_mrg_fmt_list.append(outFile)
            wh = self.checkErrorFile(outFile,wh,"error",tmp_stat_file)
        
        print >>wh, 'echo \" End of formatting of genotype of merged step \"'
        print >>wh, "\n############## End of Merge and Formatting Step ###############\n"

        return out_ovp_mrg_fmt_list,wh

    def writeDemoFamFilter(self,configDict,manifest,ngc_id,fam_id,out_ovp_mrg_list,
                                                fam_ngc_dir,tmp_stat_file,wh):

        ''' Subroutine Function to launch Rscript family filtering. '''

        print >>wh,"########### Starting Family filtering for family: ",\
                    fam_id," ############\n"
        print >>wh, "\necho \"Begin of Family-Filtering step\"\n" 
        
        resource_path = configDict['general']['resourceDir']
        script_path   = os.path.dirname(os.path.dirname(os.path.abspath( __file__ )))
        famBin        = script_path+'/'+configDict['famFilter']['famBin']
        ovFrac_list   = re.split('\,',configDict['overlapMerge']['ovFrac'])
        imprint_genes = '/'.join([resource_path,
                                  configDict['imprint']
                                 ])
        ref_genome    = configDict['general']['genomeBuild']

        filter_dir    = fam_ngc_dir+'/filter/'
        print >>wh, 'mkdir -p '+filter_dir
        
        for inpFile in out_ovp_mrg_list:
            for ovFrac in ovFrac_list:
                if re.search(str(ovFrac)+'.overlap',inpFile):
                    fam_cmd = ['Rscript ',famBin,'-v',inpFile,'-o',filter_dir,'-f',fam_id,
                                                    '-n',ngc_id,'-p',manifest,'-c',ovFrac,
                                                    '-i',imprint_genes,'-r',ref_genome
                              ]

                    print >>wh,' '.join(fam_cmd)
                    
        
        print >>wh, "echo \"End of Family-Filtering step\"\n" 
        
        print >>wh,"########### End of Family filtering for family: ",fam_id,\
                    " ############\n"

        return wh

    def checkErrorFile(self,out_file,wh,flag,status_file=[]):
        
        if_cmd = "\nif [ ! -e "+out_file+" ]; then \n"
        if flag=="error":
            echo_1 = "   echo -e \"ERROR: "+out_file+" doesnot exist\" \n"
            echo_2 = "   echo \"0\" | cat > \""+status_file+"\"\n"
            echo_3 = "   exit 1\n"
            echo_4 = "   else\n"
            echo_5 = "          echo \"1\" | cat > \""+status_file+"\""
        end_cmd = "\nfi"
        print >>wh,if_cmd+echo_1+echo_2+echo_3+echo_4+echo_5+end_cmd
        print >>wh,"\n"
        return wh

    def writeSlurmEndStatus(self,tmp_stat_file,wh):

        print >>wh,"echo \"2\" | cat > "+tmp_stat_file

        return wh


