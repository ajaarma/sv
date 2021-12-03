#########################################################
#                                                       #
# Class to process Annotation sources objects           #
#                                                       #
#                                                       #
#                                                       #
#########################################################

import re,os,gzip,sys,time
from collections import OrderedDict
from vcf_parser import VCFParser


class dbSV:

    def __init__(self,elements=[]):
        self.__elements={}
        for e in elements:
            self.__elements[e]=1

    def getFormatNGCVariant(self,variant,ngc_id,ilm_id,famOvpFracFlag=False,
                                                           queryInsOffset=0,
                                                           queryBndOffset=0):

        ''' Subroutine to process variants line in the VCF file. 
        Extracts: Chromosme Num, Start-pos, End-pos, CIPOS,CIEND,Sample ID
                  SV-Type, Genotype
        '''

        #Extract Illumina-Id 
        for keys in variant:
            if re.search(ilm_id,keys):
                ilm_id = keys
            elif re.search("\_Proband",keys):
                ilm_id = keys
        
        #Variant Information dictionary
        varInfo = variant['info_dict']
        
        #Chromosome number
        chr_num = str(variant['CHROM']);
        
        #Start Position
        start_pos = int(variant['POS'])
        if start_pos != 0:
            start_pos = start_pos - 1
        else:
            start_pos = start_pos

        #Extract variant END position
        try:
            end_pos = int(varInfo['END'][0])
        except:
            #end_pos = start_pos+1
            end_pos = start_pos

        # Caller type: Manta,Canvas or NC: No Caller
        try:
            # removes rs ids attached to SV callers
            sv_caller = re.split("\;",variant['ID'])[0]
        except:
            sv_caller = "NC"

        # Canvas related calls dont have Genotypes
        # Format style: GT:CN:FT; For Manta calls CN=.
        try:
            var_fmt = variant['FORMAT']
            if re.search("GT",var_fmt):
                gt_strs = re.split("\:",variant[ilm_id])
                gt = gt_strs[0].strip()

                if re.search("CN",var_fmt):
                    #genoType = gt_strs[0]+":"+gt_strs[1]
                    genoType = gt_strs[0]+":"+gt_strs[3]+":"+gt_strs[1]
                else:
                    genoType = gt_strs[0]+":"+"."+":"+gt_strs[1]
            else:
                gt_strs = re.split("\:",variant[ilm_id])
                genoType = str(gt_strs[0])+":"+str(gt_strs[1])+":"+str(gt_strs[2])
        except:
            var_format = "NF"

        # Assign uniform SV-type for Canvas calls: GAIN=> DUP; LOSS=>DEL; 
        try:
            sv_type = varInfo['SVTYPE'][0]
            if re.search("GAIN",sv_caller):
                sv_type = "DUP"
            elif re.search("LOSS",sv_caller):
                sv_type = "DEL"
        except:
           sv_type = "IMP"

        # For illumina-vcf files; Some INS sv-types have same start & end position.
        # Made such coordinates 1-based
        if (sv_type == "INS" or sv_type == "BND") and (start_pos == end_pos):
            end_pos = start_pos+1

        # For  NGC samples add user declared offset to INS & BND start and end 
        # points
        if (sv_type == "INS"):
            start_pos = start_pos - int(queryInsOffset)
            end_pos = end_pos + int(queryInsOffset)

        if (sv_type == "BND"):
            start_pos = start_pos - int(queryBndOffset)
            end_pos = end_pos + int(queryBndOffset)

        
        # Parsing CIPOS and CIEND; Confidence interval (Both Manta and Canvas calls)
        ci_pos_list = []
        ci_end_list = []

        try:
            ci_pos = varInfo['CIPOS']
            for e in ci_pos:
                if not re.search("\.",e):
                    e = int(e)
                    ci_pos_list.append(e)
        except:
            ci_pos_list=[0]

        try:
            ci_end = varInfo['CIEND']
            for e in ci_end:
                if not re.search("\.",e):
                    e = int(e)
                    ci_end_list.append(e)
        except:
            ci_end_list = [0]

        # Adding CIPOS and CIEND to original start and end position.
        # If start position is -ve (generally for MT coordinates) then assign as 0
        
        if famOvpFracFlag == True:
            start_pos = start_pos
            end_pos   = end_pos
        else:
            start_pos = start_pos+int(min(ci_pos_list))
            end_pos = end_pos+int(max(ci_end_list))
        
        #For MT coordinates if start pos is < 0
        if start_pos < 0:
            print start_pos
            start_pos = 0

        # Combine all the information as: Chrom, Start, end, Information 
        # field and filter out IMPRECISE SVs generally from the Canvas:REF IDs
        
        if sv_type != "IMP":
            #print start_pos, end_pos, ci_pos_list, ci_end_list
            out_line = [chr_num,str(start_pos),str(end_pos),ngc_id+"|"+sv_caller+\
                                                         "|"+sv_type+"|"+genoType
                       ]
        else:
            out_line = []
        
        return out_line
        

    def displayArguments(self,configDict,db_type,ref_genome,manifest_file,
                                proj_date,lift_over_flag):
        
        ''' Subroutine function to process individual annotation sources'''

        # Extracting general binaries
        db_type = db_type.lower()
        ref_genome = ref_genome.lower() #str(configDict['general']['genomeBuild'].lower())
        resource_path = os.path.abspath(configDict['general']['resourceDir'])
        bgzip = 'bgzip' #os.path.abspath(configDict['general']['samtools'])+'/bgzip'
        bcftools = 'bcftools' #os.path.abspath(configDict['general']['samtools'])+'/bcftools'
        tabix = 'tabix' #os.path.abspath(configDict['general']['samtools'])+'/tabix'
        extName = configDict['general']['extName']

        print "\n\nProcessing the database: ",db_type,"\n"
        if db_type=="dbvar":
            
            dbDict = configDict[db_type]
            clingen_file = dbDict[ref_genome]['clingen']
            user_file = dbDict[ref_genome]['user']
            clinvar_file = dbDict[ref_genome]['clinvar']
           
            cv_out_file = dbDict[ref_genome]['annoFile']
            
            print " --Entered CLINGEN file :",clingen_file,"\n"
            print " --Entered User curated file :",user_file,"\n"
            print " --Entered CLINVAR file :",clinvar_file,"\n"
            '''
            out_dir = os.path.split(os.path.dirname(clingen_filea)[0]
            out_file = out_dir+"/"+db_type+"_"+ref_genome+extName+"CG.UR.CV.bed"
            out_sorted_file = out_dir+"/"+db_type+"_"+ref_genome+extName+"CG.UR.CV.sorted.bed"
            out_sorted_gz = out_dir+"/"+db_type+"_"+ref_genome+extName+"CG.UR.CV.sorted.merged.bed.gz"
            '''
            cv_dir = os.path.dirname(cv_out_file)
            tmp_out_file = '/'.join([cv_dir,'tmp.cv.bed'])

            wh = open(tmp_out_file,"w")
            wh = self.processDBVar(clingen_file,wh,"CG")
            wh = self.processDBVar(user_file,wh,"UR")
            wh = self.processDBVar(clinvar_file,wh,"CV")
            wh.close()

            print " -- Sorting by chromosome and coordinates"
            self.getCoordSorted(tmp_out_file,cv_out_file,db_type)
            
            print "-- Finished processing: ",db_type+"\n"
            #os.system("sort -k1,1n -k2n "+out_file+" > "+out_sorted_file)

        elif db_type =="ngc":
        
            # family overlap flag (boolean). If TRUE then don't include CIPOS and CIEND for 
            # adjusting start and end coordinates of input SV.
            
            print manifest_file,db_type,ref_genome
            print "Processing database: ",db_type
            dbDict = configDict[db_type]
            print dbDict
            #gnomad_raw_file = '/'.join([resource_path,dbDict[ref_genome]['rawFile']])
            #gnomad_out_file = '/'.join([resource_path,dbDict[ref_genome]['annoFile']])
            
            # Naming output directory and output SV coordinates for NGC database.
            if ref_genome=="grch37":
                
                ngc_out_file = '/'.join([resource_path,dbDict[ref_genome]['annoFile']])
                ngc_basename = os.path.basename(ngc_out_file)
                print ngc_out_file
                print proj_date
                ngc_dirname = '/'.join([os.path.dirname(ngc_out_file),proj_date])
                #out_file = 'NGC.samples'+extName+'grch37.bed' #os.path.basename(configDict['ngc'][ref_genome])
                #out_dir = os.path.dirname(configDict['ngc'][ref_genome])+"/"+proj_date
                print ngc_dirname
                os.system('mkdir -p '+ngc_dirname)
                ngc_out_file = '/'.join([ngc_dirname,ngc_basename])

            elif ref_genome=="grch38":
                ngc_out_file = '/'.join([resource_path,dbDict[ref_genome]['annoFile']])
                ngc_basename = os.path.basename(ngc_out_file)
                ngc_dirname = '/'.join([os.path.dirname(ngc_out_file),proj_date])
                #out_file = 'NGC.samples'+extName+'grch37.bed' #os.path.basename(configDict['ngc'][ref_genome])
                #out_dir = os.path.dirname(configDict['ngc'][ref_genome])+"/"+proj_date
                print ngc_dirname
                os.system('mkdir -p '+ngc_dirname)
                ngc_out_file = '/'.join([ngc_dirname,ngc_basename])

            print " -- The SV coordinates for NGC samples will be written to: "+ngc_out_file
            #sys.exit()
            tmp_out_file = '/'.join([os.path.dirname(ngc_out_file),'tmp.bed.gz'])
            
            # Processing the all the VCF files in the manifest to create NGC database
            wh = open(tmp_out_file,"w")
            wh = self.processNGC(configDict,manifest_file,ngc_dirname,wh)
            wh.close()

            #sys.exit()
            # Sorting and merging the SV coordinates
            print "\n"
            print " -- Sorting by chromosome and coordinates"
            self.getCoordSorted(tmp_out_file,ngc_out_file,db_type)

            #os.system('rm '+out_file+';'+'rm '+outFile_sm_all)

            if int(lift_over_flag)==1:
                print " -- Processing Liftover from GRCh37 to GRCh38"
                self.getLiftOver_37to38(configDict,ngc_out_file,db_type,proj_date)
            print " -- Finished processing "+db_type,"\n"
        
        elif db_type =="gnomad":
            
            dbDict = configDict[db_type]
            gnomad_raw_file = '/'.join([resource_path,dbDict[ref_genome]['rawFile']])
            gnomad_out_file = '/'.join([resource_path,dbDict[ref_genome]['annoFile']])
             
            print " -- Entered GNOMAD file is: ",gnomad_out_file,"\n"
            #out_dir = os.path.dirname(gnomad_file)
            #out_file = out_dir+"/"+db_type+"_"+ref_genome+extName+"AC.AN.AF.bed"
            #out_sorted_gz = out_dir+"/"+db_type+"_"+ref_genome+extName+"AC.AN.AF.sorted.bed.gz"
            tmp_out_file = '/'.join([os.path.dirname(gnomad_out_file),'tmp.bed.gz'])
            print tmp_out_file
            print " -- Processing Raw GNOMAD File for selected fields/columns:",gnomad_raw_file
            wh = open(tmp_out_file,"w")
            wh = self.processGNOMAD(gnomad_raw_file,wh)
            wh.close()
            
            print " -- Sorting by chromosome and coordinates"
            self.getCoordSorted(tmp_out_file,gnomad_out_file,db_type)
            print " -- Finished processing ",db_type,"\n"

        elif db_type =="ensembl":

            dbDict = configDict[db_type]
            ensg_file = '/'.join([resource_path,dbDict[ref_genome]])

            print " -- Entered Ensembl file is: ",ensg_file,"\n"
            out_dir = os.path.dirname(ensg_file)

            
            ensembl_type_list = configDict['ensembl']['ensType']['offset'].keys()
            print ensembl_type_list

            for ens_type in ensembl_type_list:
                
                if re.search('transcript',ens_type):
                    ens_out_file = '/'.join([resource_path,
                                             configDict['ens_trans'][ref_genome]['annoFile']
                                            ])
                elif re.search('exons',ens_type):
                    ens_out_file  = '/'.join([resource_path,
                                              configDict['ens_exon'][ref_genome]['annoFile']
                                             ])

                offset = configDict['ensembl']['ensType']['offset'][ens_type]
                tmp_out_file = '/'.join([out_dir,'tmp.'+ens_type+'.bed'])

                print " -- Processing the Raw Ensembl file: ",ensg_file
                wh = open(tmp_out_file,"w")
                wh = self.processEnsembl(ensg_file,ens_type,int(offset),wh)
                wh.close()

                print " -- Sorting by chromosome and coordinates"
                self.getCoordSorted(tmp_out_file,ens_out_file,db_type,ens_type)
                print " -- Compressing and Indexing sorted file"
                print " -- Finished processing: ",db_type,' for ',ens_type,"\n"

        elif db_type =="refseq":

            dbDict = configDict[db_type]
            refseq_file = '/'.join([resource_path,dbDict[ref_genome]])

            print " -- Entered Refseq raw file is: ",refseq_file,"\n"
            out_dir = os.path.dirname(refseq_file)
            
            refseq_type_list = configDict['refseq']['refType']['offset'].keys()

            for ref_type in refseq_type_list:
                
                if re.search('gene',ref_type):
                    ref_out_file = '/'.join([resource_path,
                                             configDict['ref_gene'][ref_genome]['annoFile']
                                            ])
                elif re.search('exons',ref_type):
                    ref_out_file  = '/'.join([resource_path,
                                              configDict['ref_exon'][ref_genome]['annoFile']
                                             ])

                tmp_out_file = '/'.join([out_dir,'tmp.'+ref_type+'.bed'])
                offset = configDict['refseq']['refType']['offset'][ref_type]

                print " -- Processing the Raw RefSeq file: ",refseq_file,\
                                                ' and ref-type: ',ref_type
                wh = open(tmp_out_file,"w")
                wh = self.processRefSeq(refseq_file,ref_type,int(offset),wh)
                wh.close()

                #sys.exit()
                print "     -- Sorting by chromosome and coordinates"
                self.getCoordSorted(tmp_out_file,ref_out_file,db_type,ref_type)

                print "     -- Compressing and Indexing sorted file"
                print " -- Finished processing: ",db_type,' for ',ref_type,"\n"
 
        elif db_type =="decipher":

            dbDict = configDict[db_type]
            decp_raw_file = dbDict[ref_genome]['rawFile']
            decp_anno_file = dbDict[ref_genome]['annoFile']
            
            print " -- Entered Decipher file is: ",decp_anno_file,"\n"
            decp_out_dir = os.path.dirname(decp_anno_file)
            tmp_out_file = '/'.join([decp_out_dir,'tmp.decp.bed'])
            #out_file = out_dir+"/"+db_type+"."+ref_genome+extName+"bed"
            #out_sorted_gz = out_dir+"/"+db_type+"."+ref_genome+extName+"sorted.merged.bed.gz"

            print " -- Processing the Raw Decipher File: ",decp_raw_file
            wh = open(tmp_out_file,"w")
            wh = self.processDECIPHER(decp_raw_file,wh)
            wh.close()
            
            print " -- Sorting by chromosome and coordinates"
            self.getCoordSorted(tmp_out_file,decp_anno_file,db_type)
            
            print " -- Finished processing: ",db_type,"\n"

        elif db_type=='user':

            # Retrieve the configuration files for User sepcific annotation source
            dbDict = configDict[db_type]
            user_file = dbDict[ref_genome]['rawFile']
            user_out_file = dbDict[ref_genome]['annoFile']
            db_tag = dbDict['tag']

            print " -- Entered User provided annotation file is: ",user_file,'\n'
            out_dir = os.path.dirname(user_file)
            tmp_out_file = '/'.join([out_dir,'tmp.bed'])
            #out_file = out_dir+"/"+db_type+"_"+ref_genome+extName+db_tag+'.bed'
            #out_sorted_gz = out_dir+"/"+db_type+"_"+ref_genome+'.'+extName+'.'+db_tag+'.'+\
            #                                                        "sorted.merged.bed.gz"


            print " -- Processing the raw annotation file is: ",user_file,'\n'
            wh = open(tmp_out_file,'w')
            wh = self.processUSER(user_file,wh)
            wh.close()

            print " -- Sorting by chromosome and coordinates"
            self.getCoordSorted(tmp_out_file,user_out_file,db_type,db_tag)

            print " -- Finished processing: ",db_type,'\n'



    ##########################
    #                        #
    #       METHODS          #
    #                        #
    ##########################

    def processUSER(self,user_file,wh):
        
        ''' Subroutine to process the USER provided customized annotation source '''
        
        if re.search('.gz',user_file):
            fh = gzip.open(user_file)
        else:
            fh = open(user_file)

        for lines in fh:
            lines = lines.strip()
            if not re.search("^#",lines):
                strs = re.split("\t|\s",lines); strs = [str(x.strip()) for x in strs]
                
                if re.search('gene',strs[2]):
                    
                    chr_num = strs[0]
                    st_pos = int(strs[3])-1000
                    if st_pos <0:
                        st_pos = 0
                    en_pos = int(strs[3])

                    ens_gene_id = re.split("\W+",strs[strs.index(
                                                                "gene_id")+1]
                                                                )[1]
                    hgnc_gene_name = re.split("\W+",strs[strs.index(
                                                                "gene_name")+1]
                                                                )[1]
                    ens_bio_type = re.split("\W+",strs[strs.index(
                                                                "gene_biotype")+1]
                                                                )[1]

                    wh.write(str(chr_num)+"\t"+str(st_pos)+"\t"+str(en_pos)+"\t"
                             +ens_gene_id+'*'+hgnc_gene_name+'*'+ens_bio_type+"\n"
                            )

        return wh


    def processDECIPHER(self,inp_file,wh):

        fh = open(inp_file)

        for lines in fh:
            lines = lines.strip()

            if re.search("^#",lines):
                pass
            else:
                strs = re.split("\t",lines); strs = [str(x.strip()) for x in strs]
                
                chr_num = strs[0]
                start_pos = int(strs[1])
                if start_pos !=0:
                    start_pos = str(start_pos -1)
                else:
                    start_pos = str(start_pos)

                end_pos = strs[2]
                sv_type = strs[3]
                syndrome = strs[4]

                out_list = [chr_num,start_pos,end_pos,syndrome+"|"+sv_type]

                print >>wh,"\t".join(out_list)

        fh.close()
        
        return wh

    def getLiftOver_37to38(self,configDict,inp_file,db_type,proj_date):
        ''' Subroutine to convert GRCh37 coordinates to GRCh38 '''

        #bgzip = 'bgzip' #os.path.abspath(configDict['general']['samtools'])+'/bgzip'
        #bcftools = 'bcftools' #os.path.abspath(configDict['general']['samtools'])+'/bcftools'
        #tabix = 'tabix' #os.path.abspath(configDict['general']['samtools'])+'/tabix'
        resource_path = configDict['general']['resourceDir']

        lovBin = configDict[db_type]['lovbin']
        lovChain = configDict[db_type]['lovchain']
        extName = configDict['general']['extName']
       
        lov38out_file = '/'.join([resource_path,configDict['ngclov']['annoFile']])
        lov38out_dir = os.path.dirname(lov38out_file)

        lov38_tmp_out_file = '/'.join([lov38out_dir,'lov.tmp.bed.gz'])
        lov38unmap_file    = '/'.join([lov38out_dir,'lov.tmp.unmap'])
        #lov38unmap_file = os.path.abspath(lov38out_file+".unmap")

        cmd = [lovBin,"-multiple","-noSerial","-bedPlus=3","-tab",inp_file,
                                    lovChain,lov38_tmp_out_file,lov38unmap_file]
        
        print " -- Lift Over from GRCh37 to GRCh38"
        print ' '.join(cmd)
        os.system(" ".join(cmd))
    
        print " -- Sorting by Chromosome and Coordinates"
        self.getCoordSorted(lov38_tmp_out_file,lov38out_file,db_type)
        os.system("rm "+lov38unmap_file)


    def getCoordSorted(self,tmp_file,sort_file,db_type=[],ens_type=[]):

        tmp_merge_file = '/'.join([os.path.dirname(tmp_file),'tmp.merge.bed'])
        wh = open(tmp_merge_file,'w')
        self.mergeDupCoordBed(tmp_file,wh,db_type)
        wh.close()

        print 'sort -k1,1V -k2,2n -k3,3n '+tmp_merge_file+' |bgzip -c > '+sort_file
        os.system('sort -k1,1V -k2,2n -k3,3n '+tmp_merge_file+' |bgzip -c > '+sort_file)
        time.sleep(5.5)
        os.system('tabix -p bed '+sort_file)
        #os.system('rm '+tmp_merge_file)
        #os.system('rm '+tmp_file)
        #return outFile_s_all,outFile_sm_all

    def mergeDupCoordBed(self,inp_file,wh,db_type=[],ens_type=[]):

        bed_hash = OrderedDict()

        fh = open(inp_file)
              
        for lines in fh:
            lines = lines.strip()
            strs = re.split("\t",lines)
            strs = [str(x) for x in strs]
            key_id = strs[0]+"+"+strs[1]+"+"+strs[2]
            val_id = strs[3]
            if strs[0]=='chr17_ctg5':
                print lines
                sys.exit()

            if bed_hash.has_key(key_id):
                tmp = bed_hash[key_id]
                tmp.append(val_id)
                bed_hash[key_id] = list(set(tmp))
                tmp = []
            else:
                bed_hash[key_id] = [val_id]

        fh.close()

        for keys in bed_hash:
            coords = re.split("\+",keys)
            vals = bed_hash[keys]
        
            if db_type == 'user' and ens_type == 'promoter':
                
                ensg_gene = []
                hgnc_gene = []
                gene_bt = []

                for ele in vals:
                    strs = re.split("\*",ele)
                    strs = [x.strip() for x in strs]
                    ensg_gene.append(strs[0])
                    hgnc_gene.append(strs[1])
                    gene_bt.append(strs[2])

                ensg_gene_uniq = "/".join(list(set(ensg_gene)))
                hgnc_gene_uniq = "/".join(list(set(hgnc_gene)))
                gene_bt_uniq = "/".join(list(set(gene_bt)))
                vals_list = [ensg_gene_uniq,hgnc_gene_uniq,gene_bt_uniq]
                print >>wh,"\t".join(coords)+"\t"+"*".join(vals_list)
  
            elif db_type == "ensembl" and ens_type=='transcript':

                ensg_gene = []
                hgnc_gene = []
                gene_bt = []
                trans_id = []
                trans_bt = []
               
                
                for ele in vals:
                    strs = re.split("\*",ele)
                    strs = [x.strip() for x in strs]
                    ensg_gene.append(strs[0])
                    hgnc_gene.append(strs[1])
                    gene_bt.append(strs[2])
                    trans_id.append(strs[3])
                    trans_bt.append(strs[4])

                ensg_gene_uniq = "/".join(list(set(ensg_gene)))
                hgnc_gene_uniq = "/".join(list(set(hgnc_gene)))
                gene_bt_uniq = "/".join(list(set(gene_bt)))
                trans_id_uniq = "/".join(list(set(trans_id)))
                trans_bt_uniq = "/".join(list(set(trans_bt)))
                #vals_list = [exon_ids,tr_ids,ensg_ids,gene_ids,gene_bt]
                vals_list = [ensg_gene_uniq,hgnc_gene_uniq,gene_bt_uniq,
                                                 trans_id_uniq,trans_bt_uniq]
                print >>wh,"\t".join(coords)+"\t"+"*".join(vals_list)
            
            elif db_type =='ensembl' and ens_type=='exons':
                print >>wh,'\t'.join(coords)+'\t'+'/'.join(vals)

            elif db_type == 'refseq' and ens_type=='gene':
                
                ref_gene_id = []
                ref_gene_name = []
                ref_gene_bt = []
                ref_gene_syn = []

                for ele in vals:
                    strs = re.split('\*',ele)
                    strs = [x.strip() for x in strs]
                    ref_gene_id.append(strs[0])
                    ref_gene_name.append(strs[1])
                    ref_gene_syn.append(strs[2])
                    ref_gene_bt.append(strs[3])

                ref_gene_uniq = '/'.join(list(set(ref_gene_id)))
                ref_gene_name_uniq = '/'.join(list(set(ref_gene_name)))
                ref_gene_bt_uniq = '/'.join(list(set(ref_gene_bt)))
                ref_gene_syn_uniq = '/'.join(list(set(ref_gene_syn)))

                vals_list = [ref_gene_uniq,ref_gene_name_uniq,ref_gene_syn[0],
                                                             ref_gene_bt[0]]
                print >>wh,'\t'.join(coords)+'\t'+'*'.join(vals_list)

            elif db_type == 'refseq' and ens_type=='exons':
                
                ref_gene_id = []
                ref_gene_name = []
                ref_trans_id = []
                ref_exon_num = []

                for ele in vals:
                    strs = re.split('\*',ele)
                    strs = [x.strip() for x in strs]
                    ref_gene_id.append(strs[0])
                    ref_gene_name.append(strs[1])
                    ref_trans_id.append(strs[2])
                    ref_exon_num.append(strs[3])

                ref_gene_uniq = '/'.join(list(set(ref_gene_id)))
                ref_gene_name_uniq = '/'.join(list(set(ref_gene_name)))
                ref_trans_uniq = '/'.join(list(set(ref_trans_id)))
                ref_exon_uniq = '/'.join(list(set(ref_exon_num)))

                vals_list = [ref_gene_uniq,ref_gene_name_uniq,ref_trans_uniq,
                                                             ref_exon_uniq]
                print >>wh,'\t'.join(coords)+'\t'+'*'.join(vals_list) 

            elif db_type=="ngc":        
                print >>wh,"\t".join(coords)+"\t"+"$".join(vals)
            elif db_type=='user' and len(ens_type)==0:
                print >>wh,"\t".join(coords)+"\t"+"$".join(vals)
            else:
                print >>wh,"\t".join(coords)+"\t"+"$".join(vals)

        return wh

    def processGNOMAD(self,inp_file,wh):
      
        if re.search('gz$',inp_file):
            fh = gzip.open(inp_file)
        else:
            fh = open(inp_file)

        for lines in fh:
            lines = lines.strip()

            if not re.search("^#",lines):
                strs = re.split("\t",lines)
                strs = [str(x.strip()) for x in strs]
                chr_num = strs[0]
                st_pos = int(strs[1]) - 1
                if st_pos<0:
                    st_pos = 0
                
                en_pos = strs[2]
                gnomad_id = strs[3]
                sv_type = strs[4]
                ac = re.split("\,",strs[29]);ac = "/".join(ac)
                an = re.split("\,",strs[28]);an = "/".join(an)
                af = re.split("\,",strs[30]);
                af = [str(round(float(x),6)) for x in af]; af = "/".join(af)
                
                try:
                    if re.search("\,",strs[38].strip()):
                        print gnomad_id
                    popmax_af = str(round(float(strs[38]),6))
                except:
                    popmax_af = "NA"
                
                rest = [gnomad_id,sv_type,ac,an,af,popmax_af]
                                                        
                wh.write(str(chr_num)+"\t"+str(st_pos)+"\t"+str(en_pos)+"\t"+
                                                            "|".join(rest)+"\n")

        fh.close()

        return wh


    def processNGC(self,configDict,manifest_file,out_dir,wh):
        ''' Subroutine to process all the NGC samples present in Manifest file'''
      
        
        famOvpFracFlag = configDict['overlapMerge']['famOverlapFrac']
        queryInsOffset = configDict['query']['offset']['ins']
        queryBndOffset = configDict['query']['offset']['bnd']
        bgzip = 'bgzip' #os.path.abspath(configDict['general']['samtools'])+'/bgzip'
        bcftools = 'bcftools' #os.path.abspath(configDict['general']['samtools'])+'/bcftools'

        # Processing each of the NGC family ids present in the manifest file
        print manifest_file
        fh = open(manifest_file)

        count = 0
        for lines in fh:
            #Skip the header
            count = count+1

            if count >=2:
                lines = lines.strip()
                if re.search("^family_id",lines):
                    pass
                else:
                    strs = lines.split("\t")
                    ngc_id = strs[1].strip()
                    vcf_file_strs = re.split(".genome.vcf.gz",strs[10].strip())
                    vcf_file_orig_gz = vcf_file_strs[0]+".SV.vcf.gz"
                    vcf_file_norm_gz = out_dir+'/'+\
                                       os.path.basename(vcf_file_strs[0])+\
                                       '.SV.norm.vcf.gz'

                    # Normalizing SV-vcf file for splitting multi-allelic sites
                    cmd_norm = bcftools+' norm -m - '+vcf_file_orig_gz+\
                                    ' | '+bgzip+' -c > '+vcf_file_norm_gz

                    cmd_norm_index = bcftools+' index '+vcf_file_norm_gz
                   
                    os.system(cmd_norm)
                    os.system(cmd_norm_index)
                   
                    # Extracting Illumina Id
                    ilm_id = strs[6].strip()
               
                    # Processing NGC sample vcf file.
                    print " -- Processing File: "+str(ngc_id)+" "+\
                            os.path.basename(vcf_file_norm_gz)+" "+str(count)
                    my_parser = VCFParser(infile=vcf_file_norm_gz,
                                                split_variants=True,check_info=True)
               
                    # Start parsing VCF file
                    for variant in my_parser:
                        out_list = self.getFormatNGCVariant(variant,ngc_id,ilm_id,
                                                                    famOvpFracFlag,
                                                                    queryInsOffset,
                                                                    queryBndOffset
                                                           )
                        if out_list:
                            print >>wh,"\t".join(out_list)
                    
                    # Removing the normalized vcf file
                    os.system('rm '+vcf_file_norm_gz)
                    os.system('rm '+vcf_file_norm_gz+'.csi')
                    #sys.exit()

        fh.close()

        return wh

    def processEnsembl(self,inp_file,ens_type,offset,wh):

        
        fh = gzip.open(inp_file,'rb')

        ''' 
        if ens_type == 'transcript':
            print >>wh, '\t'.join(['Chromosome','Transcript_start','Transcript_end',
                             'Gene_id_Gene_name_Gene_type_Transcript_id_Transcript_biotype'
                                  ])
        elif ens_type == 'exons':
            print >>wh, '\t'.join(['Chromosome','Exon_region_start',
                                  'Exon_region_end','Exon_stable_ID'])

        '''
 
        for lines in fh:
            lines = lines.strip()

            if not re.search("^#",lines):
                strs = re.split("\t|\s",lines); strs = [x.strip() for x in strs if x]

                if ens_type== 'exons' and re.search("exon",strs[2]):
                    chr_num = strs[0]
                    st_pos = int(strs[3])-offset #8
                    en_pos = int(strs[4])+offset #8
   
                    #For zero based coordinate
                    if st_pos !=0:
                        st_pos = st_pos -1
                    else:
                        st_pos = st_pos

                    if st_pos <0:
                        st_pos = 0

                    exon_id = re.split("\W+",strs[strs.index("exon_id")+1])[1]

                    wh.write(str(chr_num)+"\t"+str(st_pos)+"\t"+str(en_pos)+"\t"+
                                                                    exon_id+"\n")

                elif ens_type == 'transcript' and re.search('transcript',strs[2]):

                    chr_num = strs[0]
                    st_pos = int(strs[3]) - offset
                    en_pos = int(strs[4]) + offset

                    #For zero based coordinate
                    if st_pos !=0:
                        st_pos = st_pos -1
                    else:
                        st_pos = st_pos

                    if st_pos <0:
                        st_pos = 0


                    ens_gene_id = re.split("\W",strs[strs.index(
                                                            "gene_id"
                                                           )+1
                                                    ]
                                          )[1]

                    hgnc_gene_name = re.split("\W",strs[strs.index(
                                                                   "gene_name"
                                                                  )+1
                                                       ]
                                             )[1]

                    gene_bio_type = re.split("\W+",strs[strs.index(
                                                                "gene_biotype"
                                                                  )+1
                                                       ]
                                            )[1]

                    transcript_id = re.split("\W",strs[strs.index(
                                                                "transcript_id"
                                                                 )+1
                                                      ]
                                            )[1] 

                    trans_bio_type = re.split("\W+",strs[strs.index(
                                                            "transcript_biotype"
                                                                 )+1
                                                        ]
                                             )[1]

                    wh.write(str(chr_num)+"\t"+str(st_pos)+"\t"+str(en_pos)+"\t"
                             +ens_gene_id+'*'+hgnc_gene_name+'*'+gene_bio_type+
                                      '*'+transcript_id+"*"+trans_bio_type+"\n"
                            )
       
        return wh

    def processRefSeq(self,inp_file,ref_type,offset,wh):

        fh = gzip.open(inp_file,'rb')

        ''' 
        if ens_type == 'transcript':
            print >>wh, '\t'.join(['Chromosome','Transcript_start','Transcript_end',
                             'Gene_id_Gene_name_Gene_type_Transcript_id_Transcript_biotype'
                                  ])
        elif ens_type == 'exons':
            print >>wh, '\t'.join(['Chromosome','Exon_region_start',
                                  'Exon_region_end','Exon_stable_ID'])

        '''
 
        for lines in fh:
            lines = lines.strip()

            if not re.search("^#",lines):
                strs = re.split("\t",lines); strs = [x.strip() for x in strs if x]

                if ref_type== 'exons' and re.search("exon",strs[2]):
                    chr_num = strs[0]
                    st_pos = int(strs[3])-offset #8
                    en_pos = int(strs[4])+offset #8
       
                    #For 0-based coordinate
                    if st_pos !=0:
                        st_pos = st_pos -1
                    else:
                        st_pos = st_pos

                    if st_pos <0:
                        st_pos = 0


                    gene_strs = re.split('\;',strs[8])
                    gene_strs = [x.strip() for x in gene_strs]
                  
                    ref_gene_strs = [x for x in gene_strs if re.search(
                                                                '^gene_id',x)]
                    ref_gene_id = [re.split('\s|\"',x)[2] 
                                                    for x in ref_gene_strs][0]
                   
                    ref_gene_strs = [x for x in gene_strs if re.search(
                                                                  '^gene ',x)]
                    ref_gene_name = [re.split('\s|\"',x)[2] 
                                                    for x in ref_gene_strs][0]

                    ref_gene_strs = [x for x in gene_strs if re.search(
                                                          '^transcript_id',x)]
                    ref_trans_id = [re.split('\s|\"',x)[2] 
                                                    for x in ref_gene_strs][0] 
                    
                    ref_gene_strs = [x for x in gene_strs if re.search(
                                                            '^exon_number',x)]
                    ref_exon_number = ['Exon_'+str(re.split('\s|\"',x)[2]) 
                                                     for x in ref_gene_strs][0] 
                    if not re.search('^H',chr_num):
                        wh.write(str(chr_num)+"\t"+str(st_pos)+"\t"+str(en_pos)+"\t"+
                                 ref_gene_id+'*'+ref_gene_name+'*'+ref_trans_id+
                                 '*'+ref_exon_number+"\n"
                                )

                elif ref_type == 'gene' and re.search('gene',strs[2]):

                    chr_num = strs[0]
                    st_pos = int(strs[3]) - offset
                    en_pos = int(strs[4]) + offset
                    
                    #For 0-based coordinate
                    if st_pos != 0:
                        st_pos = st_pos -1
                    else:
                        st_pos = st_pos

                    if st_pos <0:
                        st_pos = 0

                    gene_strs = re.split('\;',strs[8])
                    gene_strs = [x.strip() for x in gene_strs]
                  
                    ref_gene_strs = [x for x in gene_strs if re.search('^gene_id',x)]
                    ref_gene_id = [re.split('\s|\"',x)[2] for x in ref_gene_strs][0]
                   
                    ref_gene_strs = [x for x in gene_strs if re.search('^gene ',x)]
                    ref_gene_name = [re.split('\s|\"',x)[2] for x in ref_gene_strs][0]

                    ref_gene_strs = [x for x in gene_strs if re.search('^gene_biotype',x)]
                    ref_gene_biotype = [re.split('\s|\"',x)[2] for x in ref_gene_strs][0]

                    try:
                        #syn_gene_strs = re.split("\;|\s|\"",gene_strs)
                        
                        syn_strs = [re.split('\s|\"',gene_strs[i]) 
                                        for i,val in enumerate(gene_strs) 
                                                if re.search('gene_synonym',val)
                                    ]
                        if syn_strs:
                            syn_gene_name = [ x[2] for x in syn_strs ]
                        else:
                            syn_gene_name = ['NA']
                    except ValueError:
                        syn_gene_name = ['NA']
                    
                    if not re.search('^H',chr_num):
                        wh.write(str(chr_num)+"\t"+str(st_pos)+"\t"+str(en_pos)+"\t"
                             +ref_gene_id+'*'+ref_gene_name+'*'+'/'.join(syn_gene_name)+
                                      '*'+ref_gene_biotype+"\n"
                                )
                    
        return wh 


    def processDBVar(self,inp_file,wh,db_file_name):

        print " -- Processing the database: ",db_file_name
        my_parser = VCFParser(infile=inp_file,split_variants=True,check_info=True)

        for variant in my_parser:
            #count=count+1
            #if count>0:
            chr_num = variant['CHROM']
            start_pos = int(variant['POS'])
            if start_pos != 0:
                start_pos = start_pos - 1
            else:
                start_pos = start_pos
                    
            var_info = variant['info_dict']
            #print var_info
            sv_type = var_info['SVTYPE'][0]
            clinsig = var_info['CLNSIG'][0]
            dbvar_id = variant['ID']
            
            ci_pos_list = []
            ci_end_list = []
            
            try:
                ci_pos = var_info['CIPOS']
                for e in ci_pos:
                    if not re.search("\.",e):
                        e = int(e)
                        ci_pos_list.append(e)
            except:
                ci_pos_list=[0]

            try:
                ci_end = var_info['CIEND']
                for e in ci_end:
                    if not re.search("\.",e):
                        e = int(e)
                        ci_end_list.append(e)
            except:
                ci_end_list = [0]

            try:
                end_pos = int(var_info['END'][0])
            except:
                end_pos = start_pos

            start_pos = start_pos+int(min(ci_pos_list))
            end_pos = end_pos+int(max(ci_end_list))
            
            out_list = [chr_num,start_pos,end_pos,dbvar_id+"|"+sv_type+"|"+clinsig+"|"+db_file_name]
            out_list = map(str,out_list)
            print >>wh,"\t".join(out_list)
       
        return wh

    def processNGCOverlap(self,famDict,fam_id,out_dir,winLen):
        ''' Function to output: 
            (a) Parse the SV vcf file for all NGC ids of given family id.
            (b) Write toml file for use with vcfanno
            (c) Edit the affected NGC vcf(s) for the overlap '''
        
        #Initialize the directories and output file for all the action
        out_dir_fam = out_dir+"/famAnnoDB"
        os.system("mkdir -p "+out_dir_fam)
        outNGC_svFile = out_dir_fam+"/"+fam_id+".famAnnoDB.raw.bed"
        outNGC_svFile_gz = out_dir_fam+"/"+fam_id+".famAnnoDB.raw.sorted.merged.bed.gz"
        outNGC_toml = out_dir_fam+"/"+fam_id+".fam.toml"

        with open(outNGC_svFile,'w') as wh:

            ngc_ids = famDict[fam_id]['ngc_id']
            for ngc_id,ilm_id,svcf in zip(famDict[fam_id]['ngc_id'],
                                          famDict[fam_id]['sample'],
                                          famDict[fam_id]['spath']):
                print " -- Processing File: ", ngc_id," ",os.path.basename(svcf)

                my_parser = VCFParser(infile=svcf,split_variants=True,check_info=True)
                for variant in my_parser:
                    out_list = self.getFormatNGCVariant(variant,ngc_id,ilm_id)
                    if out_list:
                        print >>wh,"\t".join(out_list)

        print " -- Sorting by chromosome and coordinates"
        outNGC_sm_all = self.getCoordSorted(outNGC_svFile)

        print " -- Compressing and Indexing sorted file "
        os.system("module load samtools/1.3; bgzip -c "+outNGC_sm_all+ " > "+
                  outNGC_svFile_gz+"; tabix -p bed "+outNGC_svFile_gz)
        
        print " -- Finished creating the annotation Database for family id: "+fam_id+"\n"

        #Editing the coordinates of the affected NGC samples for computing internal overlap
        for ngc_id,svcf,af_status in zip(famDict[fam_id]['ngc_id'],famDict[fam_id]['spath'], 
                                         famDict[fam_id]['affected']):
            if int(af_status) == 2:
                inp_svFile = svcf
                inp_svFileName = re.split("\.SV.vcf.gz",os.path.basename(inp_svFile))[0]

                inp_svFile = out_dir+"/"+inp_svFileName+"."+ngc_id+".SV.edit.vcf.gz"
                #out_svEditFile = out_dir+"/"+inp_svFileName+"."+ngc_id+".SV.edit."+str(winLen)+".vcf.gz"
                #out_svEditFile = out_dir+"/"+ngc_id+".SV.edit."+str(winLen)+".vcf.gz"
                out_svEditFile = out_dir+"/"+ngc_id+".SV.edit."+"fam.vcf.gz"
                
                print " -- Editing the VCF file and adding "+str(winLen)+"bp for family id: "+fam_id+"\n"

                fh = gzip.open(inp_svFile)
                wh = gzip.open(out_svEditFile,"wb")

                for lines in fh:
                    lines = lines.strip()
                    if re.search("^#",lines):
                        print >>wh,lines
                    else:
                        strs = re.split("\t",lines);strs = [x.strip() for x in strs]
                        
                        chr_num = strs[0]
                        start_pos = int(strs[1]) - int(winLen)
                        if (start_pos) < 0:
                            start_pos = start_pos
                        else:
                            start_pos = start_pos

                        sv_id = strs[2]
                        ref = strs[3]

                        alt_id = strs[4]
                        if re.search("DEL",sv_id):
                            sv_type = "DEL"
                            if re.search("DEL",alt_id):
                                alt_id = alt_id
                            else:
                                alt_id = "<DEL:"+alt_id+">"

                        #if re.search("MantaBND",sv_id):
                        #    alt_id = "<INS:"+alt_id+">"
                        #elif re.search("MantaINS",sv_id):
                        #    alt_id = "INS:"+alt_id+">"
                        
                        qual = strs[5]
                        filter_tag = strs[6]
                         
                        info_strs = re.split("\;",strs[7].strip())
                       
                        for i,e in enumerate(info_strs):
                            if re.search("^END",e):
                                e_strs = re.split("\=",e)
                                end_pos = int(e_strs[1])+int(winLen)
                                info_strs[i] = "END="+str(end_pos)
                                
                            if re.search("SVLEN",e):
                                e_strs = re.split("\=",e)
                                sv_len = abs(int(e_strs[1]))+2*int(winLen)
                                sv_len = (-1)*sv_len
                                info_strs[i] = "SVLEN="+str(sv_len)

                        gt_tag = strs[8].strip()
                        gt_val = strs[9].strip()
                        out_line = [chr_num,str(start_pos),sv_id,ref,alt_id,qual,
                                    filter_tag,";".join(info_strs),gt_tag,gt_val]
                        print >>wh,"\t".join(out_line)
                        

                wh.close()
                fh.close()

                wh = open(outNGC_toml,"w")
                print >>wh,'[[annotation]]'
                print >>wh,'file=\"'+outNGC_svFile_gz+'\"'
                print >>wh,'names=[\"fam_s\",\"fam_e\",\"fam_ov_id\"]'
                print >>wh,'columns=[2,3,4]'
                print >>wh,'ops=[\"self\",\"self\",\"self\"]'

                wh.close()
                            

        return outNGC_svFile

