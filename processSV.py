#!/usr/bin/python

############################################################################################
# Description: Program to process structural variants (SVs) of NGC Trios called
#              by Illumina v38 pipeline using Manta,Canvas. Generates shell scripts that can
#              be run on SLURM/LSF clusters.
#
# Author: Ajay A. Kumar 
#
#
###########################################################################################


import re,sys,getopt,os,datetime,subprocess
from collections import OrderedDict

global script_path 
script_path = os.path.dirname(os.path.abspath( __file__ ))

sys.path.append(script_path+"/CONFIG/")
sys.path.append(script_path+"/SLURM/")
sys.path.append(script_path+"/VCFANNO/")
from CONFIG import *
from SLURM import *
from VCFANNO import *


if __name__=="__main__":


    chr_list = [str(x) for x in range(1,23)]+['X','Y','MT']

    # Get the command line arguments
    objC = CONFIG()

    # Parse Input command line arguments
    cmd_dict = objC.parseSVCommandArgs()
    manifest = cmd_dict['manifest']; xml_file = cmd_dict['xml']; 
    proj_date = cmd_dict['proj'];  anal_type = cmd_dict['analType']; 
    work_dir = cmd_dict['workDir'];fam_file  = cmd_dict['famFile']
    genome_ref = cmd_dict['ref']
    launch_flag = cmd_dict['launch']

    # Create the temporary directories
    cmd_args = sys.argv
    tmp_dict = objC.processInit(work_dir,sys.argv,proj_date)
    tmp_dir = tmp_dict['tmpDir']; tmp_bin = tmp_dict['tmpBin']; 
    tmp_data = tmp_dict['tmpData']
    tmp_status = tmp_dict['tmpStat']; tmp_log = tmp_dict['tmpLog']; 
    sb_log = tmp_dict['sbLog']

    # Creating object from CONFIG class and processing Config or Analysis file.
    objC = CONFIG()
    configDict = objC.getConfigDict(xml_file)

    #Overidding XML default reference genome based on user input.
    if genome_ref == configDict['general']['genomeBuild']:
        pass
    else:
        configDict['general']['geneomeBuild'] = genome_ref

    ''' 
    # Printing the details of the configuration file for AnnotSV/vcfanno
    print " -- Parameters involved for Structural Variant Annotations using: "+anal_type+"\n"
    for keys in configDict[anal_type]:
        print "\t",keys,"\t",configDict[anal_type][keys]
    ''' 
    print "\n"

    slurm_str = [anal_type]
    
    # Create SLURM object to access the methods in SLURM Class
    objS = SLURM()

    # Create SLURM job per family member
    fh = open(fam_file)

    # Create master dictionary that store all the information present in manifest 
    # file with keys being family id
    famDict = objS.getManifestDict(manifest)

    # Create Toml file
    objV = VCFANNO()
    configDict = objV.createTomlFile(configDict,genome_ref,tmp_dir)

    # Processing for each family ID
    for fam_id in fh:
    #with open(fam_file) as fh:
        
        fam_id  = fam_id.strip()
        fam_id  = re.split("\t|\s",fam_id)[0].strip()
        fam_dir = os.path.abspath(tmp_data+"/"+fam_id)
        for ngc_id, af_status, svcf, ilm_id in zip(famDict[fam_id]['ngc_id'],
                                                   famDict[fam_id]['affected'],
                                                   famDict[fam_id]['spath'],
                                                   famDict[fam_id]['sample']):
            if int(af_status) == 2:
                #print ngc_id
                #print ngc_id,af_status,svcf,ilm_id

                fam_ngc_dir = os.path.abspath(fam_dir+"/"+ngc_id)
                #os.system("mkdir -p "+fam_ngc_dir)

                tmp_stat_file = tmp_status+"/"+ngc_id+".Job_status_"+slurm_str[0]+".txt"
                slurm_file, swh = objS.getSlurmWriteHandle(tmp_bin,ngc_id+"."+slurm_str[0])

                #SLURM specific parameters

                swh = objS.writeSlurmTop(configDict,swh)
               
                script_out, script_err, swh = objS.writeSlurmInit(configDict,
                                                            sb_log,ngc_id,swh)
                #swh = objS.writeSlurmSpecific(swh)
                #swh = objS.writeSlurmModule(configDict,swh)

                print >>swh, 'mkdir -p '+fam_ngc_dir

                if anal_type=="annotsv":
                    fam_ngc_dir,ngc_list,swh = objS.writeSlurmAnnotSV(configDict,
                                                              famDict,tmp_data,
                                                              fam_id,ngc_id,
                                                              tmp_stat_file,j,swh
                                                                     )
                    
                    swh = objS.writeSlurmCombChrOut(configDict,fam_ngc_dir,ngc_id,
                                                                tmp_stat_file,swh
                                                   )
                elif anal_type=="vcfanno_ngc":
                                        
                    #Step-1: Write VCFAnno related command in Slurm
                    out_svFile,swh = objS.writeSlurmVcfanno(manifest,configDict,
                                                           famDict, fam_ngc_dir,
                                                           fam_id, ngc_id, svcf,
                                                          ilm_id, tmp_stat_file,
                                                     script_path,genome_ref,swh
                                                           )

                    #Step-2: Extract relevant fields related to annotation sources
                    out_overlap,swh = objS.writeSlurmExtractFields(configDict,
                                                       manifest,ngc_id,fam_id,
                                                     out_svFile,tmp_stat_file,
                                                               genome_ref,swh
                                                                  )
                    '''
                    out_overlap = os.path.abspath(fam_ngc_dir)+'/annoDB'
                    ''' 
                    #Step-3: Overlap percentage merging
                    out_ovp_mrg_list,swh = objS.writeSlurmOverlapMerge(configDict,
                                                           manifest,ngc_id,fam_id,
                                                             xml_file,out_overlap,
                                                         tmp_stat_file,genome_ref,
                                                                    anal_type,swh
                                                                      )

                    #Optional: sub-step 3; Format GT if they are missing specially
                    # for Canvas related calls
                    if configDict['formatGT']['flag']=='True':
                        out_ovp_mrg_list,swh = objS.writeSlurmFormatGT(configDict,
                                                    manifest,fam_id,out_ovp_mrg_list,
                                                                 tmp_stat_file,swh
                                                                      
                                                                      )
                    #Step-4: Family filtering 
                    swh = objS.writeSlurmFamFilter(configDict,manifest,ngc_id,
                                                      fam_id,out_ovp_mrg_list,
                                                      fam_ngc_dir,tmp_stat_file,swh
                                                  )
                
                elif anal_type=='vcfanno_demo':
                    objD = DEMO()
                    print 'here'
                    #Step-1: Write VCFAnno related command in Slurm
                    out_svFile,swh = objD.writeDemoVcfanno(manifest,configDict,
                                                           famDict, fam_ngc_dir,
                                                           fam_id, ngc_id, svcf,
                                                          ilm_id, tmp_stat_file,
                                                     script_path,genome_ref,swh
                                                           )

                    #Step-2: Extract relevant fields related to annotation sources
                    out_overlap,swh = objD.writeDemoExtractFields(configDict,
                                                       manifest,ngc_id,fam_id,
                                                     out_svFile,tmp_stat_file,
                                                               genome_ref,swh
                                                                  )
                    '''
                    out_overlap = os.path.abspath(fam_ngc_dir)+'/annoDB'
                    ''' 
                    #Step-3: Overlap percentage merging
                    out_ovp_mrg_list,swh = objD.writeDemoOverlapMerge(configDict,
                                                           manifest,ngc_id,fam_id,
                                                             xml_file,out_overlap,
                                                         tmp_stat_file,genome_ref,
                                                                    anal_type,swh
                                                                      )
                    ''' 
                     
                    #Optional: sub-step 3; Format GT if they are missing specially
                    # for Canvas related calls
                    if configDict['formatGT']['flag']=='True':
                        out_ovp_mrg_list,swh = objD.writeDemoFormatGT(configDict,
                                                    manifest,fam_id,out_ovp_mrg_list,
                                                                 tmp_stat_file,swh
                                                                      
                                                                      )
                    #Step-4: Family filtering 
                    swh = objD.writeDemoFamFilter(configDict,manifest,ngc_id,
                                                      fam_id,out_ovp_mrg_list,
                                                      fam_ngc_dir,tmp_stat_file,swh
                                                  )
                    '''
                    

                swh.close()
                if launch_flag:
                    #print 0
                    os.system("sbatch "+ slurm_file) 
    
    fh.close()
