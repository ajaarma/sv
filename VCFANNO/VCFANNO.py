#!/usr/bin/python

import re,os,gzip,sys
from collections import OrderedDict


class VCFANNO:
   
    #global VCFANNO-class variable to hold number of NGC sample count
    ngc_id_list = []

    def __self__(self,elements=[]):
        self.__elements={}
        for e in elements:
            self.__elements[e] = 1
            
    def query_db_overlap(self,q_start,q_end,db_start,db_end,flag=[]):

        # Computes the length of overlap between query and db
        diff = min(db_end,q_end) - max(db_start,q_start)+1

        # Computes only the % of query overlap
        query_ov = float(diff)/float(q_end-q_start+1)
        
        # Computes only the % of query overlap
        #db_ov = float(diff)/float(db_end - db_start+1) # deafult
        db_ov = float(diff)/float(db_end - db_start)    # O-based 
        query_ov = round(query_ov,3)
        db_ov = round(db_ov,3)

        return query_ov,db_ov


    def getQueryDBOverlap(self,st_pos,en_pos,db_s_list,db_e_list,db_info_list,
                                                    ovFrac=[],debug_flag=[]):

        db_sp_list = []; db_ep_list=[];db_infop_list=[]
        qd_list = []; dq_list=[]

        for i in range(0,len(db_s_list)):
            try:
                dbs = int(db_s_list[i])
                dbe = int(db_e_list[i])
            except ValueError:
                break

            print "ALL: ",st_pos,en_pos,dbs,dbe #query_ov,db_ov #,db_info_list[i]
            query_ov,db_ov = self.query_db_overlap(st_pos,en_pos,dbs,dbe)

            #if debug_flag:
                #print "ALL: ",st_pos,en_pos,dbs,dbe,query_ov,db_ov #,db_info_list[i]
            if query_ov >=ovFrac and db_ov >=ovFrac:
                db_sp_list.append(str(dbs))
                db_ep_list.append(str(dbe))
                db_infop_list.append(str(db_info_list[i]))
                qd_list.append(str(query_ov))
                dq_list.append(str(db_ov))
                #print "FILT: ",st_pos,en_pos,dbs,dbe,query_ov,db_ov,db_info_list[i]
                #if debug_flag:
                #    print "ALL: ",st_pos,en_pos,dbs,dbe,query_ov,db_ov,db_info_list[i]

        #print qd_list,dq_list
        if len(qd_list) ==0 and len(dq_list) ==0:
            db_sp_list=["NA"];db_ep_list=["NA"];db_infop_list=["NA"];
            qd_list=["NA"];dq_list=["NA"]

        return db_sp_list,db_ep_list,db_infop_list,qd_list,dq_list

    def getQueryDBFamOverlap(self,st_pos,en_pos,ci_pos,ci_end,db_s_list,
                             db_e_list,db_info_list,ovFrac,debug_flag=[]):

        db_sp_list = []; db_ep_list=[];db_infop_list=[]
        qd_list = []; dq_list=[]

        for i in range(0,len(db_s_list)):
            try:
                dbs = int(db_s_list[i])
                dbe = int(db_e_list[i])
            except ValueError:
                break

            stM_low = st_pos+int(min(ci_pos))
            stM_up  = st_pos+int(max(ci_pos))

            enM_low = en_pos+int(min(ci_end))
            enM_up  = en_pos+int(max(ci_end))
            
            if debug_flag:
                print "ST-CI ",stM_low,stM_up,st_pos,ci_pos,dbs
                print "EN-CI ",enM_low,enM_up,en_pos,ci_end,dbe
                print "\n"

            if (stM_low <= dbs <= stM_up) and (enM_low <= dbe <= enM_up):
                #query_ov,db_ov = self.query_db_overlap(st_pos,en_pos,dbs,dbe)
                query_ov,db_ov = self.query_db_overlap(st_pos,en_pos,dbs,dbe)
                
                #if debug_flag:
                #    print "ALL: ",stM_low,stM_up,st_pos,en_pos,dbs,dbe,query_ov,db_ov #,db_info_list[i]
                #print "ALL: ",st_pos,en_pos,dbs,dbe,query_ov,db_ov #,db_info_list[i]
                #if query_ov >=ovFrac and db_ov >= ovFrac:
                
                if 0 == 0:    
                    #print "ALL: ",st_pos,en_pos,dbs,dbe,query_ov,db_ov #,db_info_list[i]
                    db_sp_list.append(str(dbs))
                    db_ep_list.append(str(dbe))
                    db_infop_list.append(str(db_info_list[i]))
                    qd_list.append(str(query_ov))
                    dq_list.append(str(db_ov))
                    #print "FILT: ",st_pos,en_pos,dbs,dbe,query_ov,db_ov,db_info_list[i]
                    #if debug_flag:
                    #    print "FAM-ALL: ",st_pos,en_pos,dbs,dbe,query_ov,db_ov,db_info_list[i]


        #print qd_list,dq_list
        if len(qd_list) ==0 and len(dq_list) ==0:
            db_sp_list=["NA"];db_ep_list=["NA"];db_infop_list=["NA"];
            qd_list=["NA"];dq_list=["NA"]

        return db_sp_list,db_ep_list,db_infop_list,qd_list,dq_list          
    
    def getGT(self,variant):

        for ele in variant.keys():
            if re.search("Proband",ele):
                gt = variant[ele]

        return gt

    def getDBList(self,variant,db_type):
        try:
            db_s_list = var_info[db_type+'_s']
            db_e_list = var_info[db_type+'_e']
            db_ov_list = var_info[db_type+'_ov_id']
        except:
            db_s_list = ["NA"]
            db_e_list = ["NA"]
            db_ov_list = ["NA"]

        return db_s_list,db_e_list,db_ov_list

    def printOverlap(self,st_pos,en_pos,db_st_list,db_en_list,db_ov_list,
                                        qDB=[],DBq=[],db_type=[],ov_flag=[]):
            

        if ov_flag==1:
            print db_type+" overlap: ",st_pos,en_pos,db_st_list,db_en_list,db_ov_list,qDB,DBq
        else:
            print db_type+" raw: ",st_pos,en_pos,db_st_list,db_en_list,db_ov_list

    def getAlldbOverlap(self,st_pos,en_pos,stM_pos,enM_pos,strs,ovFrac,db_type,
                                         sv_type,ngc_fam_flag=[],debug_flag=[]):

        #if db_type.lower() in ["gnomad","ngc38","ngclov","dbvar","ensembl",
        if db_type.lower() in ['gnomad','ngc38','ngclov','ngc37','dbvar',
                               'ens_trans','ens_exon','ref_gene','ref_exon',
                                                         'decipher','user']:
            i = 3
            li = i+3
            ri = li+3

            '''
            #if db_type in ["dbVar","ensembl","decipher"]:
            if db_type.lower() in ['dbvar','ens_trans','ens_exon','ref_gene',
                                                'ref_exon','decipher','user']:
                ovFrac = 0.0
            '''
        else:
            print "The entered database doesnot exist. Please check the command\n"
            sys.exit()
        

        db_s_list = re.split("\,",strs[i])
        db_e_list = re.split("\,",strs[i+1])
        db_ov_list = re.split("\,",strs[i+2])

        db_ls_list = re.split("\,",strs[li])
        db_le_list = re.split("\,",strs[li+1])
        db_lov_list = re.split("\,",strs[li+2])
        
        db_rs_list = re.split("\,",strs[ri])
        db_re_list = re.split("\,",strs[ri+1])
        db_rov_list = re.split("\,",strs[ri+2])

        qDB =[]; DBq = [];

        #printOverlap(st_pos,en_pos,db_s_list,db_e_list,db_ov_list,qDB,DBq,db_type,ov_flag=0)
        db_ovp_merge = []

        if sv_type!="BND" and sv_type!="INS":
            if debug_flag:
                print "Inside : ",debug_flag
                print "Inside func: ",st_pos,en_pos,stM_pos,enM_pos
                #print db_s_list
                #print db_ov_list
                #sys.exit()

            db_sp_list,db_ep_list,\
            db_ovp_list,qDB,DBq = self.getQueryDBOverlap(st_pos,en_pos,db_s_list,
                                                             db_e_list,db_ov_list,
                                                               ovFrac, debug_flag)

            db_lsp_list,db_lep_list,\
            db_lovp_list,lqDB,DBql = self.getQueryDBOverlap(stM_pos,enM_pos,db_ls_list,
                                                         db_le_list,db_lov_list,ovFrac,
                                                                            debug_flag)
            db_rsp_list,db_rep_list,\
            db_rovp_list,rqDB,DBqr = self.getQueryDBOverlap(stM_pos,enM_pos,db_rs_list,
                                                         db_re_list,db_rov_list,ovFrac,
                                                                            debug_flag)

            if ngc_fam_flag:
                db_ovp_merge = list(set(db_ovp_list))
            else:
                db_ovp_merge = list(set(db_ovp_list+db_lovp_list+db_rovp_list))

        else:
            if debug_flag:
                print "Inside : ",debug_flag
                print "Inside func: ",st_pos,en_pos,stM_pos,enM_pos
                print db_rs_list
                print db_re_list
                #print db_ov_list
                #sys.exit()

            db_sp_list = db_s_list
            db_ep_list = db_e_list
            db_ovp_list =  ['NA'] if '.' in db_ov_list else db_ov_list
            db_lovp_list = ['NA'] if '.' in db_lov_list else db_lov_list
            db_rovp_list = ['NA'] if '.' in db_rov_list else db_rov_list
            
            if ngc_fam_flag:
                db_ovp_merge = list(set(db_ovp_list))
            else:
                db_ovp_merge = list(set(db_ovp_list+db_lovp_list+db_rovp_list))
            
            qDB = ["NA"]
            DBq = ["NA"]

        #db_ovp_merge = self.rmDotNa(db_ovp_merge)
        #print 'db ovlap merge: ',db_ovp_merge 
        return db_ovp_merge
        #return db_s_list,db_e_list,db_ov_list,",".join(db_sp_list),",".join(db_ep_list),",".join(db_ovp_list),",".join(qDB),",".join(DBq)

    #def rmDotNa(self,dv_ovp_merge):

    def getFamilyOverlap(self,st_pos,en_pos,ci_pos,ci_end,strs,db_type,sv_type,
                                                            ovFrac,debug_flag=[]):

        ''' Subroutine to compute family overlap '''

        #coordiantes are 0-based (original file column number - 1)
        if db_type.lower() in ["ngc38","ngc37"]:
            i = 3
            li = i+3
            ri = li+3
        else:
            print "The entered database type:",db_type," doesnot exist. Please \
                                                            check the command\n"
            sys.exit()
        
        db_s_list = re.split("\,",strs[i])
        db_e_list = re.split("\,",strs[i+1])
        db_ov_list = re.split("\,",strs[i+2])

        db_ls_list = re.split("\,",strs[li])
        db_le_list = re.split("\,",strs[li+1])
        db_lov_list = re.split("\,",strs[li+2])
        
        db_rs_list = re.split("\,",strs[ri])
        db_re_list = re.split("\,",strs[ri+1])
        db_rov_list = re.split("\,",strs[ri+2])

        qDB =[]; DBq = [];

        #printOverlap(st_pos,en_pos,db_s_list,db_e_list,db_ov_list,qDB,DBq,db_type,ov_flag=0)
        db_ovp_merge = []

        # PROCESS SV-Types: DEL,DUP,INV,
        if sv_type!="BND" and sv_type!="INS":
        #if 0==0:
            if debug_flag:
                print "Inside FAM: ",debug_flag
                print "Inside FAM-func: ",st_pos,en_pos,ci_pos,ci_end
                print db_ls_list
                print db_le_list

            db_sp_list,db_ep_list,db_ovp_list,qDB,DBq = self.getQueryDBFamOverlap(
                  st_pos,en_pos,ci_pos,ci_end,db_s_list,db_e_list,db_ov_list,ovFrac,
                                                                        debug_flag)
            db_lsp_list,db_lep_list,db_lovp_list,lqDB,DBql = self.getQueryDBFamOverlap(
                  st_pos,en_pos,ci_pos,ci_end,db_ls_list,db_le_list,db_lov_list,ovFrac,
                                                                            debug_flag)
            db_rsp_list,db_rep_list,db_rovp_list,lqDB,DBql = self.getQueryDBFamOverlap(
                  st_pos,en_pos,ci_pos,ci_end,db_rs_list,db_re_list,db_rov_list,ovFrac,
                                                                            debug_flag)
            #db_ovp_merge = list(set(db_ovp_list+db_lovp_list+db_rovp_list))
            db_ovp_merge = list(set(db_ovp_list))

        # PROCESS: BND and INS
        else:
            if debug_flag:
                print "Inside Fam-func: ",db_ls_list
                print "Inside Fam-func: ",db_le_list
            db_sp_list = db_s_list
            db_ep_list = db_e_list
            db_ovp_list = db_ov_list
            db_lovp_list = db_lov_list
            db_rovp_list = db_rov_list
            db_ovp_merge = list(set(db_ovp_list+db_lovp_list+db_rovp_list))
            qDB = ["NA"]
            DBq = ["NA"]

        #db_ovp_merge = self.rmDotNa(db_ovp_merge)
            
        return db_ovp_merge
     

    def getSVType(self,caller_id):

        caller_strs = re.split("\:",caller_id)
        caller_strs = [str(x) for x in caller_strs]
        #print caller_strs
        sv_type = []

        if caller_strs[1]=="GAIN":
            sv_type = "DUP"
        elif caller_strs[1]=="LOSS":
            sv_type = "DEL"
        elif caller_strs[1]=="REF":
            sv_type = "IMP"

        return sv_type
    

    def writeDBSpecificHeader(self,db_type,wh):

        sample_header = ["#CHROM","POS","END","CIPOS","CIEND","MPOS","MEND",
                                                        "SV","SVLEN","SVID","GT"]
        gnomad_header = ["GNOMAD_ID","S_GN_FRAC","GN_S_FRAC","GN_SAME","GN_ALL"]
        ngc38_header = ["NGC_ID","S_NGC_FRAC","NGC_S_FRAC","NGC_SAME","NGC_SAME_AF",
                                "NGC_ALL","NGC_ALL_AF","NGC_TRIO","NGC_TRIO_AF"]
        ngcLOV_header = ["NGCL_ID","S_NGCL_FRAC","NGCL_S_FRAC","NGCL_SAME",
                                         "NGCL_SAME_AF","NGCL_ALL","NGC_ALL_AF"]
        dbvar_header = ["DBVAR_ID","S_DBVAR_FRAC","DBVAR_S_FRAC"]
        ensg_header = ["ENSG_ID","S_ENSG_FRAC","ENSG_S_FRAC"]
        dc_header = ["DC_ID","S_DC_FRAC","DC_S_FRAC"]

        if db_type=="gnomad":
            header = sample_header+gnomad_header
        elif db_type == "ngc38":
            header = sample_header+ngc38_header
        elif db_type == "ngcLOV":
            header = sample_header+ngcLOV_header
        elif db_type == "dbvar":
            header = sample_header+dbvar_header
        elif db_type == "ensembl":
            header = sample_header+ensg_header
        elif db_type == "decipher":
            header = sample_header+dc_header
        elif db_type == "all":
            header = sample_header+gnomad_header+ngc38_header+ngcLOV_header+\
                                            dbvar_header+ensg_header+dc_header

        print >>wh,"\t".join(header)

        return header,wh

    def getDBstrsStartIndex(self,db_type):

        if db_type=="gnomad":
            index = 9
        
        return index

    def getDBMainLeftRightOverlapDiff(self,strs,db_type):

        ind = self.getDBstrsStartIndex(db_type)
        #print strs[11]
        if not re.search("^\.",strs[ind+2]): 
            db_main = set(re.split("\,",strs[ind+2]))
        else:
            db_main = set([])
        if not re.search("^\.",strs[ind+3]): 
            db_left = set(re.split("\,",strs[ind+3]))
        else:
            db_left = set([])
        if not re.search("^\.",strs[ind+4]): 
            db_right = set(re.split("\,",strs[ind+4]))
        else:
            db_right = set([])
        #print db_main
        #print strs
        #print strs[ind+2],"\n"
        #print db_main
        #print strs[ind+3],"\n"
        #print db_left
        #print strs[ind+4],"\n"
        #print db_right

        db_main_left = list(db_main - db_left)
        db_main_right = list(db_main - db_right)
        db_left_main = list(db_left - db_main)
        db_right_main = list(db_right - db_main)

        #print "\n"
        #print db_main_left
        #print db_main_right
        #print count,exp,line_index
        if len(gnomad_left_main) !=0 or len(gnomad_right_main) !=0:
            #print strs
            print count,sv_id,len(gnomad_main),len(gnomad_main_left),\
                            len(gnomad_main_right),len(gnomad_left_main),\
                                                    len(gnomad_right_main)

        return list(db_main),list(db_main_left),list(db_main_right),list(db_left_main),list(db_right_main)


    def getRestDict(self,rest_svFile):

        restDict = {}

        fh = open(rest_svFile)

        for lines in fh:
            lines = lines.strip()
            strs = re.split("\t",lines); strs = [x.strip() for x in strs]
            sv_id = strs[0]
            restDict[sv_id] = lines

        return restDict

    def getFieldInfo(self,lines,ensembl_flag):
        
        ''' Subroutine to extract Field information'''
        
        lines = lines.strip()
        strs = re.split("\t",lines); strs = [str(x.strip()) for x in strs]

        sv_id = strs[0]
        sv_id_strs = re.split("\+",sv_id)

        chrom_num = sv_id_strs[0]
        start_pos = int(sv_id_strs[1])
        sv_type = sv_id_strs[3]
        caller_id = sv_id_strs[4]
        
        if re.search("Canvas",caller_id):
            sv_type = self.getSVType(caller_id)
        
        try:
            end_pos = int(sv_id_strs[2])
        except:
            end_pos = start_pos+1

        try:
            if ensembl_flag!=0:
                ci_pos = re.split("\,",strs[1])
            else:
                ci_pos = re.split("\,",strs[5])
            ci_pos = [int(x) for x in ci_pos]
        except:
            ci_pos = [0]

        try:
            if ensembl_flag !=0:
                ci_end = re.split("\,",strs[2])
            else:
                ci_end = re.split("\,",strs[6])

            ci_end = [int(x) for x in ci_end]
        except:
            ci_end = [0]

        #print "CIPOS: ",ci_pos,ci_end

        stM_pos = int(start_pos)+int(min(ci_pos))
        enM_pos = int(end_pos)+int(max(ci_end))

        ci_pos = [str(x) for x in ci_pos]
        ci_end = [str(x) for x in ci_end]

        gnomad_raw_ovp = strs[5]

        #return [sv_id,chrom_num,start_pos,end_pos,stM_pos,enM_pos,ci_pos[0],ci_end[0],sv_type,strs]
        return [sv_id,chrom_num,start_pos,end_pos,stM_pos,enM_pos,ci_pos,ci_end,
                                                                    sv_type,strs]

    def getGnomadAcAf(self,gnomad_list):

        af_list = []
        ac_list = []
    
        gn_strs =  gnomad_list
        for e in gn_strs:
            e_strs = re.split("\|",e)
            if re.search("\/",e_strs[4]):
                af_strs = re.split("\/",e_strs[4]); af_strs = [float(x) for x in af_strs]
                ac_strs = re.split("\/",e_strs[2]); ac_strs = [float(x) for x in ac_strs]
                af_list.append(min(af_strs))
                ac_list.append(min(ac_strs))
            else:
                af_list.append(e_strs[4])
                ac_list.append(e_strs[2])

        af_min = [str(min(af_list))] if af_list else ['NA']
        ac_min = [str(min(ac_list))] if ac_list else ['NA']
        
        return ac_min,af_min

    def getUserAnnoOverlap(self,configDict,db_ovp_list,db_type):

        ''' Subroutine to process User defined annotation source. Current 
            implementation is suited for any overlap to promoter region
        '''
        out_list = []
        db_tag = configDict[db_type]['tag']
        db_tag_type = '_'.join([db_type,db_tag])

        prom_flag = []
        db_ovp_list = [x.strip() for x in db_ovp_list if x !='NA']
        

        if re.search('promoter',db_tag_type):
            
            ensg_gene_list=[];hgnc_gene_list=[];gene_bt_list=[]
            for ovp_id in db_ovp_list:
                strs = re.split('\*',ovp_id);
                strs = [x.strip() for x in strs]
                try:
                    ensg_gene_list.append(strs[0].upper())
                    hgnc_gene_list.append(strs[1].upper())
                    gene_bt_list.append(strs[2])
                except IndexError:
                    ensg_gene_list.append('NA')
                    hgnc_gene_list.append('NA')
                    gene_bt_list.append('NA')
           
            out_list = '*'.join(['/'.join(list(set(ensg_gene_list))),
                        '/'.join(list(set(hgnc_gene_list))),
                        '/'.join(list(set(gene_bt_list))),
                       ])
          
            if len(db_ovp_list)==0:
                out_list = 'NA'
                prom_flag = 'FALSE'
            if len(db_ovp_list)==1:
                if re.search('^NA$',db_ovp_list[0]):
                    prom_flag = 'FALSE'
                else:
                    prom_flag = 'TRUE'
            elif len(db_ovp_list)>1:
                prom_flag = 'TRUE'

            if re.search('protein_coding',out_list):
                prom_flag = 'TRUE'
            else:
                prom_flag = 'FALSE'

            out_list = [out_list,prom_flag]
            
            return out_list


    def getGeneOverlap(self,configDict,db_ovp_list,curated_gene_dict,hpo_gene_dict,
                                                            imp_gene_dict,db_type):
        
        ''' Method to process Gene overlap information for Ensembl Transcript;\
            Ensembl Exon. Refseq Gene and Exons '''


        out_list = []
        
        if db_type == 'ens_trans':
            ensg_gene_list = []; hgnc_gene_list = []; trans_id_list  = []; 
            gene_bt_list = []; trans_bt_list = []; prot_flag_list = []
            #curated_hgnc_list_flag = []; hpo_hgnc_list_flag = []

            for ovp_id in db_ovp_list:
                strs = re.split('\*',ovp_id);strs = [x.strip() for x in strs]
                try:
                    ensg_gene_list.append(strs[0].upper())
                    hgnc_gene_list.append(strs[1].upper())
                    gene_bt_list.append(strs[2])
                    trans_id_list.append(strs[3].upper())
                    trans_bt_list.append(strs[4])

                except IndexError:
                    ensg_gene_list.append('NA')
                    hgnc_gene_list.append('NA')
                    gene_bt_list.append('NA')
                    trans_id_list.append('NA')
                    trans_bt_list.append('NA')
                
            
            curated_hgnc_list_flag = []; hpo_hgnc_list_flag = []
            imp_hgnc_list_flag = []
           
            for gene_id in list(set(hgnc_gene_list)):
                if curated_gene_dict.has_key(gene_id.upper()):
                    curated_hgnc_list_flag.append('TRUE')
                else:
                    curated_hgnc_list_flag.append('FALSE')

                if hpo_gene_dict.has_key(gene_id.upper()):
                    hpo_hgnc_list_flag.append('TRUE')
                else:
                    hpo_hgnc_list_flag.append('FALSE')

                if imp_gene_dict.has_key(gene_id.upper()):
                    imp_hgnc_list_flag.append('TRUE')
                else:
                    imp_hgnc_list_flag.append('FALSE')

           
            out_list = ['/'.join(list(set(ensg_gene_list))),
                        '/'.join(list(set(hgnc_gene_list))),
                        '/'.join(list(set(gene_bt_list))),
                        '/'.join(list(set(trans_id_list))),
                        '/'.join(list(set(trans_bt_list))),
                        #'/'.join(list(set(curated_hgnc_list_flag))),
                        '/'.join(list(curated_hgnc_list_flag)),
                        #'/'.join(list(set(hpo_hgnc_list_flag)))
                        '/'.join(list(hpo_hgnc_list_flag)),
                        '/'.join(list(imp_hgnc_list_flag))
                       ]
            
            if re.search('protein_coding',out_list[2]) or \
                                re.search('protein_coding',out_list[3]):
                
                prot_flag = ['TRUE']
            else:
                prot_flag = ['FALSE']

            out_list = out_list+prot_flag
            
            return out_list

        elif db_type == 'ref_gene':
            hgnc_gene_list = []; syn_gene_list = []; 
            gene_bt_list = []; prot_flag_list = []

            for ovp_id in db_ovp_list:
                strs = re.split('\*',ovp_id);strs = [x.strip() for x in strs]
                try:
                    hgnc_gene_list.append(strs[0].upper())
                    hgnc_gene_list.append(strs[1].upper())
                    syn_gene_list.append(strs[2].upper())
                    gene_bt_list.append(strs[3])

                except IndexError:
                    hgnc_gene_list.append('NA')
                    syn_gene_list.append('NA')
                    gene_bt_list.append('NA')
                
            curated_hgnc_list_flag = []; hpo_hgnc_list_flag = []
            imp_hgnc_list_flag = []
           
            for raw_gene_id in list(set(hgnc_gene_list+syn_gene_list)):
                gene_strs = re.split('\/',raw_gene_id)
                for gene_id in gene_strs:
                    if curated_gene_dict.has_key(gene_id.upper()):
                        curated_hgnc_list_flag.append('TRUE')
                    else:
                        curated_hgnc_list_flag.append('FALSE')

                    if hpo_gene_dict.has_key(gene_id.upper()):
                        hpo_hgnc_list_flag.append('TRUE')
                    else:
                        hpo_hgnc_list_flag.append('FALSE')

                    if imp_gene_dict.has_key(gene_id.upper()):
                        imp_hgnc_list_flag.append('TRUE')
                    else:
                        imp_hgnc_list_flag.append('FALSE')
  
            out_list = ['/'.join(list(set(hgnc_gene_list))),
                        '/'.join(list(set(syn_gene_list))),
                        '/'.join(list(set(gene_bt_list))),
                        '/'.join(list(curated_hgnc_list_flag)),
                        '/'.join(list(hpo_hgnc_list_flag)),
                        '/'.join(list(imp_hgnc_list_flag))
                       ]
            
            if re.search('protein_coding',out_list[2]) or \
                                re.search('protein_coding',out_list[2]):
                
                prot_flag = ['TRUE']
            else:
                prot_flag = ['FALSE']

            out_list = out_list+prot_flag
            
            return out_list

        elif db_type == 'ens_exon':

            # counting number of exon overlaps
            ensg_exon_list = []
            out_list = []

            for ovp_id in db_ovp_list:
                strs = re.split('\*',ovp_id);strs = [x.strip() for x in strs]
                #print strs
                try:
                    ensg_exon_list.append(strs[0])
                except IndexError:
                    ensg_exon_list.append('NA')

            out_list = list(set(ensg_exon_list))
            if 'NA' in ensg_exon_list:
                exon_count = len(out_list)-1
            else:
                exon_count = len(out_list)
            
            out_cmd = ['|'.join(out_list),str(exon_count)]
            
            return out_cmd

        elif db_type == 'ref_exon':
            
            hgnc_gene_list = []; trans_id_list = []
            exon_id_list = []

            out_exon_list = list(set(db_ovp_list))

            if 'NA' in out_exon_list:
                exon_count = len(out_exon_list)-1
            else:
                exon_count = len(out_exon_list)
           
            gene_hash = OrderedDict()

            for ovp_id in db_ovp_list:
                strs = re.split('\*',ovp_id);strs = [x.strip() for x in strs]
                try:
                    gene_id = strs[0]
                    trans_id = strs[2]
                    exon_id = strs[3]
                except IndexError:
                    gene_id = 'NA'
                    trans_id = 'NA'
                    exon_id = 'NA'

                if gene_hash.has_key(gene_id):
                    tmp_exon = gene_hash[gene_id]['exons']
                    tmp_exon.append(exon_id)

                    tmp_trans = gene_hash[gene_id]['trans']
                    tmp_trans.append(trans_id)

                    gene_hash[gene_id] = {'trans':list(set(tmp_trans)),
                                          'exons':list(set(tmp_exon))
                                         }
                else:
                    gene_hash[gene_id] = {'trans':[trans_id],'exons':[exon_id]}
            
            exon_count_list = []
            exon_list = []

            for gene_id in gene_hash:
                
                if gene_id !='NA':
                    exon_count_list.append( len(gene_hash[gene_id]['exons']))
                    exon_list.append(gene_id+'*'+'/'.join(gene_hash[gene_id]['trans'])+
                                             '*'+'/'.join(gene_hash[gene_id]['exons'])
                                    )
                else:
                    exon_list.append('NA')
                    exon_count_list.append(0)

            exon_count = sum(exon_count_list)
            
            out_list = ['|'.join(exon_list)]+[str(exon_count)]

            return out_list
            
        

    def getSameSVTypeDB(self,db_ovp_list,db_type,query_svtype,fam_id_list):

        db_sv_list = []
        fam_sv_list = []
        #fam_base_id = re.split("_",ngc_id)[0]
        #fam_base_id = re.split("_",ngc_id)[0]

        if db_type in ["gnomad","ngc38","ngclov","dbvar","decipher","fam",
                                                                  "ngc37"]:

            db_strs = [x for x in db_ovp_list if x !="NA" and x !="."]
        
            #
            for db_sv_id in db_strs:
                db_sv_id_strs = re.split("\$",db_sv_id)
                
                for e in db_sv_id_strs:
                    if re.search(query_svtype,e):
                        db_sv_list.append(e)
                        for sample_id in fam_id_list:
                            if re.search(sample_id,e):
                            #if re.search(fam_base_id,e):
                                fam_sv_list.append(e)

                    elif re.search("MCNV",e):
                        db_sv_list.append(e)
        
                #return ",".join(list(set(db_sv_list)))
            if db_type in ["ngc38","ngc37"]:
                if db_sv_list:
                    if fam_sv_list:
                        return list(set(db_sv_list)),list(set(fam_sv_list))
                    else:
                        return list(set(db_sv_list)),["NA"]
                else:
                    db_sv_list = ["NA"]
                    fam_sv_list = ["NA"]
                    return db_sv_list,fam_sv_list

            elif db_type=="gnomad":
                if db_sv_list:
                    ac_min,af_min = self.getGnomadAcAf(db_sv_list)
                    return list(set(db_sv_list)),ac_min,af_min
                else:
                    db_sv_list = ["NA"]
                    af_min = ["NA"]
                    ac_min =["NA"]
                    return list(set(db_sv_list)),ac_min,af_min
            else:
                if db_sv_list:
                    return list(set(db_sv_list))
                else:
                    db_sv_list = ["NA"]
                    return list(set(db_sv_list))


    def getFileNameTask(self,out_dir,ngc_id,task,db_type,ovFrac):

        if task=="overlap":
            out_file = out_dir+"/"+ngc_id+"."+db_type+"."+str(ovFrac)+".overlap.same.bed.gz"
        elif task=="freq":
            out_file = out_dir+"/"+ngc_id+"."+db_type+".freq.bed.gz"

        return out_file

    def getSVCount(self,db_sv_list):

        sv_list = []
        if not "NA" in db_sv_list:
            for e in db_sv_list:
                strs = re.split("\|",e); strs = [x for x in strs if x]
                ngc_id = strs[0]
                sv_list.append(ngc_id)
            ngc_num = len(list(set(sv_list)))
        else:
            ngc_num = 0

        return ngc_num

    def getParentsGT(self,fam_merge_m_same,ngc_id,ngc_mother,ngc_father):

        tmp_hash = {}

        for i,e in enumerate(fam_merge_m_same):
            if e != "NA":
                ngc_strs = re.split("\|",e)
                #print e,ngc_strs
                tmp_id = ngc_strs[0]
                gt = ngc_strs[3]

                if tmp_hash.has_key(tmp_id):
                    tmp_gt = tmp_hash[tmp_id]
                    tmp_gt.append(gt)
                    tmp_hash[tmp_id] = tmp_gt
                else:
                    tmp_hash[tmp_id] = [gt]
                
        proband_gt = []
        mother_gt = []
        father_gt = []
        try:
            proband_gt = "-".join(tmp_hash[ngc_id])
        except:
            proband_gt = "NA"
        try:
            father_gt = "-".join(tmp_hash[ngc_father])
        except:
            father_gt = "NA"

        try:
            mother_gt = "-".join(tmp_hash[ngc_mother])
        except:
            mother_gt = "NA"

        return proband_gt,mother_gt,father_gt 

    
    def executeTask(self,configDict,field_list,task,debug_flag,db_type,lineCount,
                           famDict,ngc_id,ovFrac,curated_gene_dict,hpo_gene_dict,
                                                                imp_gene_dict,wh):
        
        ''' Method to execute task of finding overlaps w.r.t annotation 
            sources 
        '''
       
        fam_id_list = []
        fam_base_id = re.split("_",ngc_id)[0]
        ngc_fam_list = famDict[fam_base_id]['ngc_id']
        ngc_father = famDict[fam_base_id]['father_id'][0]
        ngc_mother = famDict[fam_base_id]['mother_id'][0]
        #print 'fam size: ',famDict[fam_base_id]['fam_size']
        if 'quad' in famDict[fam_base_id]['fam_size']:
            ngc_sib = list(set(ngc_fam_list) - set([ngc_id,ngc_father,ngc_mother]))[0]
            fam_id_list = [ngc_id,ngc_father,ngc_mother,ngc_sib]
        else:
            fam_id_list = [ngc_id,ngc_father,ngc_mother]
            ngc_sib = 'NA'


        sv_id     = field_list[0]
        chrom_num = field_list[1]
        st_pos    = field_list[2]
        en_pos    = field_list[3]
        stM_pos   = field_list[4]
        enM_pos   = field_list[5]
        ci_pos    = field_list[6]
        ci_end    = field_list[7]
        sv_type   = field_list[8]
        strs      = field_list[9]
       
        ngc_fam_flag = []

        ovFrac = float(ovFrac)

        if sv_type !="IMP":

            if task =="overlap":
                
                if db_type=="gnomad":
                    
                    # Get gnomad SVs with reciprocal overlap (rof) percentage
                    db_merge_m = self.getAlldbOverlap(st_pos,en_pos,stM_pos,
                                                               enM_pos,strs,
                                                             ovFrac,db_type,
                                                       sv_type,ngc_fam_flag,
                                                                 debug_flag
                                                     )
                    # Get SV ids with same SV types
                    db_merge_m_same,gnomad_ac,gnomad_af = self.getSameSVTypeDB(
                                                             db_merge_m,db_type,
                                                             sv_type,fam_id_list
                                                          )
                    
                    out_cmd = [str(lineCount),sv_id,",".join(db_merge_m_same),
                    str(len(db_merge_m_same)) if not "NA" in db_merge_m_same else str(0),
                    ",".join(gnomad_ac),",".join(gnomad_af)]
                    print >>wh,'\t'.join(out_cmd)


                elif db_type in ["ngc38","ngc37"]:
                    
                    ######### 70% Overlap with NGC-v38/v37 cohort ###############
                    # Compute SV overlaps with Internal NGC-38, NGC-37 cohort
                    # CIPOS and CIEND are added to start and end position of SVs
                    db_merge_m_70 = self.getAlldbOverlap(st_pos,en_pos,stM_pos,
                                                             enM_pos,strs,0.70,
                                                               db_type,sv_type,
                                                               True,debug_flag
                                                        )
 
                    if debug_flag:
                        print 'Debug 1'
                        print db_merge_m_70
                    # Extract list of overlappng SVs with same SV-types.
                    # CIPOS and CIEND are added to start and end position of SVs
                    # global ngc-full cohort; family internal overlap 
                    db_merge_m_70_same,fam_merge_m_same = self.getSameSVTypeDB(
                                                             db_merge_m_70,db_type,
                                                             sv_type,fam_id_list)
                   
                    if debug_flag:
                        print 'Debug 2'
                        print db_merge_m_70_same
                        print fam_merge_m_same
                        sys.exit()
                    ############################################################
                    
                    ############### overlap fraction [0.0,1.0] #################
                    # Compute SV overlaps with Internal NGC-38, NGC-37 cohort
                    # CIPOS and CIEND are added to start and end position of SVs
                    db_merge_m = self.getAlldbOverlap(st_pos,en_pos,stM_pos,
                                                               enM_pos,strs,
                                                             ovFrac,db_type,
                                                               sv_type,True,
                                                                 debug_flag
                                                     )

                    # Extract list of overlappng SVs with same SV-types.
                    # CIPOS and CIEND are added to start and end position of SVs
                    # global ngc-full cohort; family internal overlap 
                    db_merge_m_same,fam_merge_m_same  = self.getSameSVTypeDB(
                                                             db_merge_m,db_type,
                                                             sv_type,fam_id_list)
                   
                    if debug_flag:
                        print 'Debug 3'
                        print db_merge_m_same
                        print fam_merge_m_same
                        sys.exit()
                    ############################################################ 
                    ''' 
                    # Compute SV internal family overlap
                    # Internal family overlap is computed when the SV end points fall
                    # exactly within the CIPOS and CIEND
                    ci_merge_m = self.getFamilyOverlap(st_pos,en_pos,ci_pos,ci_end,
                                                            strs, db_type,sv_type,
                                                            ovFrac,debug_flag) 
                  
                    # Extract list of overlapping SVs with same SV types.
                    # Internal family overlap is computed when the SV end points fall
                    # exactly within the CIPOS and CIEND
                    db_ci_merge_m, fam_ci_merge_m_same = self.getSameSVTypeDB(
                                                             ci_merge_m,db_type,
                                                             sv_type,ngc_id)
                    '''  
                    # Total count of overlapping SV types of same SV types
                    db_var_num = self.getSVCount(db_merge_m_70_same) 
                    
                    # Allele frequency within cohort
                    af_ngc_cohort = float(db_var_num)/float(len(VCFANNO.ngc_id_list))

                    # Count of overlapping SV counts within family
                    fam_var_num = self.getSVCount(fam_merge_m_same)
                    #fam_var_num = self.getSVCount(fam_ci_merge_m_same)

                    # Extract Genotypes from the family overlap
                    tmp_hash = {}
                    for i,e in enumerate(fam_merge_m_same):
                    #for i,e in enumerate(fam_ci_merge_m_same):
                        if e != "NA":
                            ngc_strs = re.split("\|",e)
                            #print e,ngc_strs
                            tmp_id = ngc_strs[0]
                            gt = ngc_strs[3]

                            if tmp_hash.has_key(tmp_id):
                                tmp_gt = tmp_hash[tmp_id]
                                tmp_gt.append(gt)
                                tmp_hash[tmp_id] = tmp_gt
                            else:
                                tmp_hash[tmp_id] = [gt]
                    
                    proband_gt = []
                    mother_gt = []
                    father_gt = []
                    sib_gt = []

                    try:
                        proband_gt = "-".join(tmp_hash[ngc_id])
                    except:
                        sv_id_strs = re.split('\+',sv_id)
                        proband_gt = sv_id_strs[5]+':NA'
                    try:
                        father_gt = "-".join(tmp_hash[ngc_father])
                    except:
                        father_gt = "NA"
                    try:
                        mother_gt = "-".join(tmp_hash[ngc_mother])
                    except:
                        mother_gt = "NA"
                    try:
                        sib_gt = "-".join(tmp_hash[ngc_sib])
                    except:
                        sib_gt = "NA"

                    out_cmd = [str(lineCount),sv_id,str(db_var_num),
                                                 str(af_ngc_cohort),
                              #",".join(fam_ci_merge_m_same),str(fam_var_num),
                                         ",".join(fam_merge_m_same),
                                                   str(fam_var_num),
                                               proband_gt,mother_gt,
                                                    father_gt,sib_gt
                              ]

                    print >>wh,"\t".join(out_cmd)

                elif db_type=="ngclov":
                    
                    db_merge_m = self.getAlldbOverlap(st_pos,en_pos,stM_pos,
                                                        enM_pos,strs,ovFrac,
                                                            db_type,sv_type,
                                                               ngc_fam_flag,
                                                                 debug_flag
                                                     )
                    db_merge_m_same = self.getSameSVTypeDB(db_merge_m,db_type,
                                                           sv_type,fam_id_list
                                                          )
                    db_var_num = self.getSVCount(db_merge_m_same)

                    out_cmd = [str(lineCount),sv_id,",".join(db_merge_m_same),
                            str(db_var_num),str(float(db_var_num)/float(847))]
                    
                    print >>wh,"\t".join(out_cmd)

                elif db_type == "fam":

                    sv_id_strs = re.split("\+",sv_id)
                    sv_id_strs[1] = str(int(sv_id_strs[1])+1000)
                    sv_id_strs[2] = str(int(sv_id_strs[2])-1000)
                    sv_id_orig_coord = "+".join(sv_id_strs)
                    
                    db_merge_m = re.split("\,",strs[5])
                    db_merge_m_same = self.getSameSVTypeDB(db_merge_m,db_type,
                                                           sv_type,fam_id_list)
                    proband_gt,mother_gt,father_gt = self.getParentsGT(
                                db_merge_m_same, ngc_id, ngc_mother, ngc_father)
                    
                    out_cmd = [str(lineCount),str(sv_id_orig_coord),sv_id,
                               ",".join(db_merge_m_same),proband_gt,mother_gt,
                               father_gt]
                    print >>wh,"\t".join(out_cmd)

                elif db_type in ['dbvar','decipher']:

                    #Get all the overlapping db SV-ids with rof=0.0%
                    db_merge_m = self.getAlldbOverlap(st_pos,en_pos,stM_pos,
                                                          enM_pos,strs,ovFrac,
                                                              db_type,sv_type,
                                                                 ngc_fam_flag,
                                                                   debug_flag
                                                       )

                    db_merge_m_same = self.getSameSVTypeDB(db_merge_m,db_type,
                                                           sv_type,fam_id_list
                                                          )

                    # Get all the overlapping db SV-ids with rof=70% or user 
                    # provided
                    db_merge_m_70 = self.getAlldbOverlap(st_pos,en_pos,stM_pos,
                                                                  enM_pos,strs,
                                                                  0.70,db_type,
                                                          sv_type,ngc_fam_flag,
                                                                    debug_flag
                                                        )
                    db_merge_m_same_70 = self.getSameSVTypeDB(db_merge_m_70,
                                                                    db_type,
                                                                    sv_type,
                                                                 fam_id_list
                                                             )

                    out_cmd = [str(lineCount),sv_id,','.join(db_merge_m_same),
                                                    ','.join(db_merge_m_same_70)
                              ]
                    print >>wh,'\t'.join(out_cmd)

                elif db_type in ['ens_trans','ens_exon','ref_gene','ref_exon']:
                    
                    ovFrac = 0.0 
                    db_merge_m = self.getAlldbOverlap(st_pos,en_pos,stM_pos,
                                                               enM_pos,strs,
                                                             ovFrac,db_type, 
                                                       sv_type,ngc_fam_flag,
                                                                 debug_flag
                                                     )

                    out_str = self.getGeneOverlap(configDict,db_merge_m,
                                                      curated_gene_dict,
                                                          hpo_gene_dict,
                                                          imp_gene_dict,
                                                                db_type
                                                 )

                    out_cmd = [str(lineCount),sv_id,'\t'.join(out_str)]
                    print >>wh,'\t'.join(out_cmd)

                elif db_type in ['user']:
                    ovFrac=0.0
                    db_merge_m = self.getAlldbOverlap(st_pos,en_pos,stM_pos,
                                                        enM_pos,strs,ovFrac,
                                                            db_type,sv_type,
                                                               ngc_fam_flag,
                                                                 debug_flag
                                                     )

                    
                    out_str = self.getUserAnnoOverlap(configDict,db_merge_m,
                                                                    db_type
                                                     )
                    out_cmd = [str(lineCount),sv_id,'\t'.join(out_str)]
                    print >>wh,'\t'.join(out_cmd)

                elif db_type == "ensembl":
                    if len(db_merge_m) ==1 and (db_merge_m[0] == "." or 
                                                db_merge_m[0] == "NA"):
                        db_len = 0
                    else:
                        db_len = len(db_merge_m)
                    
                    prot_flag = True if re.search("protein_coding", 
                                                  ",".join(db_merge_m)
                                                 ) else False
                    
                    out_cmd = [str(lineCount),sv_id,",".join(db_merge_m),
                                                    str(db_len),prot_flag]
                    print >>wh,"\t".join(out_cmd)
                else:
                    out_cmd = [str(lineCount),sv_id,",".join(db_merge_m_same)]
                    print >>wh,"\t".join(out_cmd)

                '''
                if debug_flag:
                    print  lineCount,sv_id,"\t",",".join(db_merge_m_same),\
                    "\t",len(db_merge_m_same) if not "NA" in db_merge_m_same else 0
                    print "\n",lineCount,sv_id,"\t",ci_pos,"\t",ci_end,"\t",\
                                     sv_type,"\t",",".join(db_merge_m),"\t",\
                                                              db_merge_m_same
                '''

        return wh


    def processAnnotatedVCF(self,configDict,inp_file,ngc_id,db_type,famDict,
                                         line_index,debug_flag,ovFrac,task):

        for keys in famDict:
            VCFANNO.ngc_id_list = VCFANNO.ngc_id_list+famDict[keys]['ngc_id']

        VCFANNO.ngc_id_list = list(set(VCFANNO.ngc_id_list))

        
        out_dir = os.path.dirname(os.path.dirname(inp_file))+"/"+task
        os.system("mkdir -p "+out_dir)
        out_file = self.getFileNameTask(out_dir,ngc_id,task,db_type,ovFrac)

        lineCount = 0
        
        # Read the file names for Curated gene-list and HPO-genelist
        curated_gene_file = configDict['ensembl']['ensVarAnnot']['g'] 
        hpo_gene_file = configDict['ensembl']['ensVarAnnot']['d'] 
        imp_gene_file = configDict['imprint']


        # Process the files and store in the genes
        curated_gene_dict = {}
        hpo_gene_dict = {}
        imp_gene_dict = {}

        fh = open(curated_gene_file)
        for lines in fh:
            lines = lines.strip()
            curated_gene_dict[lines.upper()] = 1
        fh.close()

        fh = open(hpo_gene_file)
        for lines in fh:
            lines = lines.strip()
            if not re.search('\#',lines):
                strs = re.split('\t',lines)
                hpo_gene = strs[3].strip()
                hpo_gene_dict[hpo_gene.upper()] = 1
        fh.close()

        ##### Dict: Imprinted genes ##############
        fh = open(imp_gene_file)
        for lines in fh:
            lines = lines.strip()
            if not re.search('\#',lines):
                strs = re.split('\t',lines)
                imp_gene = strs[0].strip()
                imp_par = strs[1].strip()
                imp_gene_dict[imp_gene] = imp_par
        fh.close()

 
        fh = gzip.open(inp_file)
        wh = gzip.open(out_file,"wb")

        # Start processing the Input annotation file
        lineCount = 0
 
        for i, lines in enumerate(fh):
            lineCount = i+1

            if debug_flag:
                exp = "".join([str(lineCount),"==",str(line_index)])
                exp_eval = eval(exp)
            else:
                exp = "".join([str(lineCount),">=",str(line_index)])
                exp_eval = eval(exp)

            if exp_eval:
                field_list  = self.getFieldInfo(lines,ensembl_flag=1)
                
                if debug_flag:
                    wh = self.executeTask(configDict,field_list,task,debug_flag,
                                               db_type,lineCount,famDict,ngc_id, 
                                                       ovFrac,curated_gene_dict,
                                                 hpo_gene_dict,imp_gene_dict,wh
                                         )
                    break

                else:
                    wh = self.executeTask(configDict,field_list,task,debug_flag,
                                               db_type,lineCount,famDict,ngc_id,
                                                       ovFrac,curated_gene_dict,
                                                 hpo_gene_dict,imp_gene_dict,wh
                                         )
                

        wh.close()
    
    def getLatestNGCDbDate(self,ngc38_dir):

        dir_list = os.listdir(ngc38_dir)
        if 'gmon.out' in dir_list:
            os.system('rm -rf '+ngc38_dir+'/gmon.out')

        os.chdir(ngc38_dir)
        proj_date = max(os.listdir(ngc38_dir),key=os.path.getmtime)
        return proj_date

    def createTomlFile(self,configDict,ref_genome,tmp_dir,ngc_db_date=[]):
        
        ''' Method to create Toml file with respective annotation sources.
            Input: configDict, genome built version GRCh37/GRCh38, tmp directory
            Output: Creates toml for specfic project date '''
            
        resource_path = os.path.abspath(configDict['general']['resourceDir'])
        annoDBList    = configDict['annoDB']['db']['value']

        if re.search('37',ref_genome):
            #dbPathDict = configDict['annoDB']['dbPath']['GRCh37']
            toml_out = tmp_dir+'/'+configDict['vcfanno']['toml37']
            wh = open(toml_out,'w')

        elif re.search('38',ref_genome):
            annoDBList.insert(2,'ngclov')
            toml_out = tmp_dir+'/'+configDict['vcfanno']['toml38']
            wh = open(toml_out,'w')

        #for keys,values in dbPathDict.items():
        for db_type in annoDBList:
            
            print >>wh,'#'+db_type.upper()
            print >>wh,"[[annotation]]"
            
            # For the NGC cohort choose the latest DB creation date else specified
            if db_type =='ngc':# or keys=='ngc37':
                
                values = '/'.join([resource_path,
                                   configDict[db_type][ref_genome]['annoFile']
                                  ])
                ngc_dir = os.path.dirname(values)
                ngc_basename = os.path.basename(values)

                if ngc_db_date:
                    proj_date = ngc_db_date
                else:
                    proj_date = self.getLatestNGCDbDate(ngc_dir)
               
                values = '/'.join([ngc_dir,proj_date,ngc_basename])
                if re.search('37',ref_genome):
                    db_type = 'ngc37'
                elif re.search('38',ref_genome):
                    db_type = 'ngc38'

                print >>wh,"file=\""+values+"\""

            elif db_type =='ngclov':
                values = '/'.join([resource_path,
                                   configDict[db_type]['annoFile']
                                  ])
                print >>wh,"file=\""+values+"\""
            elif db_type =='gnomad':
                values = '/'.join([resource_path,
                                   configDict[db_type][ref_genome]['annoFile']
                                  ])
                print >>wh,"file=\""+values+"\""
            else:
                values = '/'.join([resource_path,
                                   configDict[db_type][ref_genome]['annoFile']
                                  ])
                print >>wh,"file=\""+values+"\""

            print >>wh,"names=[\""+db_type+"_s\",\""+db_type+"_e\",\""+db_type+"_ov_id\"]"
            print >>wh,"columns=[2,3,4]"
            print >>wh,'ops=[\"self\",\"self\",\"self\"]'
            print >>wh,"\n"

        wh.close()
        
        if re.search('37',ref_genome):
            configDict['vcfanno']['toml37'] = toml_out
        elif re.search('38',ref_genome):
            configDict['vcfanno']['toml38'] = toml_out

        return configDict

