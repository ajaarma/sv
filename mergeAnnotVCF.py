#!/usr/bin/python

import re,sys,os,gzip
from collections import OrderedDict
import json

global script_path
script_path = os.path.dirname(os.path.abspath(__file__))

sys.path.append(script_path+"/CONFIG/")
sys.path.append(script_path+"/SLURM/")

from CONFIG import *
from SLURM import *



if __name__=="__main__":

    #print sys.argv
    objC = CONFIG()
    cmd_dict = objC.parseMergeAnnotCommandArgs()
    
    manifest = cmd_dict['manifest']; rest_file = cmd_dict['rest'];
    anno_dir = cmd_dict['annoDir']; fam_id = cmd_dict['fam']
    genome_ref = cmd_dict['ref']; ovFrac = str(cmd_dict['overlap'])

    objS = SLURM()
    famDict = objS.getManifestDict(manifest)

    #inp_rest = sys.argv[1]
    #anno_dir = os.path.abspath(sys.argv[2])
    #fam_id = sys.argv[3]
    #genome_ref = str(sys.argv[4])
    #ovFrac = str(sys.argv[5])

    dir_list = os.listdir(anno_dir)
    out_dir = os.path.dirname(anno_dir)

    rest_hash = OrderedDict()
    gn_hash = OrderedDict()
    ngc_hash = OrderedDict()
    #ngc_fam_hash = OrderedDict()

    ngclov_hash = OrderedDict()
    dbvar_hash =OrderedDict()
 
    decipher_hash =OrderedDict()
    ensembl_hash = OrderedDict()

    ens_trans_hash = OrderedDict()
    ens_exon_hash = OrderedDict()
    ref_gene_hash = OrderedDict()
    ref_exon_hash = OrderedDict()
    user_prom_hash = OrderedDict()

    wh = gzip.open(out_dir+"/"+fam_id+".merged.all."+ovFrac+".overlap.bed.gz","wb")

    header_list = []
    '''
    rest_header = ["SV_ID","REF","ALT","QUAL_SCORE","FILTER","CIPOS","CIEND",
                                  "SVLEN","CNVLEN","SVLEN_ALL","GT","SVTYPE"
                  ]
    '''

    a1 = True
    
    for ele in dir_list:

        fh = []
        if re.search('ngc38',ele) or re.search('ngc37',ele):# or \
            #re.search('gnomad',ele) or re.search('ngclov',ele):
            ele_strs = re.split('\.',ele)
            if re.search(ele_strs[1]+'.'+str(ovFrac),ele):
                #ele = ele_strs[0]+'.ngc38.'+str(ovFrac)+'.overlap.same.bed.gz'
                fh = gzip.open(anno_dir+"/"+ele)
            else:
                fh = []
        else:
            fh = gzip.open(anno_dir+"/"+ele)
       
        
        if fh:
            print ele
            for i,lines in enumerate(fh):
                lines = lines.strip()
                strs = re.split("\t",lines);strs = [x for x in strs if x]
                sv_id = strs[1].strip()

                if re.search("gnomad",ele):
                    gn_hash[sv_id] = strs[2:]
                    
                if re.search("ngc38",ele):
                    ngc_hash[sv_id] = strs[2:]

                if re.search("ngc37",ele):
                    ngc_hash[sv_id] = strs[2:]

                if re.search("fam",ele):
                    ngc_fam_hash[sv_id] = strs[3:]

                if re.search("ngclov",ele):
                    ngclov_hash[sv_id] = strs[2:]

                if re.search("decipher",ele):
                    decipher_hash[sv_id] = strs[2:]
                    
                if re.search("dbvar",ele):
                    #strs = re.split("\t",lines)
                    #sv_id = re.split("\t|\s",strs[0].strip())[1]
                    #sv_id = strs[1]
                    dbvar_hash[sv_id] = strs[2:]
                
                if re.search('ens_trans',ele):
                    ens_trans_hash[sv_id] = strs[2:]
                if re.search('ens_exon',ele):
                    ens_exon_hash[sv_id] = strs[2:]
                if re.search('ref_gene',ele):
                    ref_gene_hash[sv_id] = strs[2:] 
                if re.search('ref_exon',ele):
                    ref_exon_hash[sv_id] = strs[2:] 
                if re.search('user',ele):
                    user_prom_hash[sv_id] = strs[2:]
           
            fh.close()
            
        if re.search("ensembl",ele):
            fh = gzip.open(anno_dir+"/"+ele)
            for i, lines in enumerate(fh):
                lines = lines.strip()
                strs = re.split("\t|\s",lines);strs = [x for x in strs if x]
                sv_id = strs[0].strip()
                ensembl_hash[sv_id] = strs[4:]
            fh.close()


    # Processing REST file
    fh = gzip.open(rest_file)
    for i,lines in enumerate(fh):
        lines = lines.strip()
        strs = re.split("\t|\s",lines);strs = [x for x in strs if x]
        sv_id = strs[0].strip()
        # Filter out all the REF tagged variants
        if not re.search("REF",sv_id):
            rest_hash[sv_id] = strs[1:]
    fh.close()

    # Initializing the header information
    rest_header = ["SV_ID","CHROM","START_POS","END_POS","SVTYPE_ORIG","CALLER_ID",
                         "REF","ALT","QUAL_SCORE","FILTER","CIPOS","CIEND","SVLEN",
                                                "CNVLEN","SVLEN_ALL","GT","SVTYPE"
                  ]
    gn_header   = ["GN_SAME_ID","GN_SAME_COUNT","GN_AC","GN_AF"]
    
    if re.search("38",genome_ref):
        ngc_header = ["NGC38_ALL_SAME_COUNT","NGC38_ALL_SAME_FREQ",
                        "NGC38_FAM_SAME_ID","NGC38_FAM_SAME_COUNT",
                                "PROBAND","MOTHER","FATHER","SIB"
                     ]
        ngclov_header = ["NGC37_ALL_SAME_ID","NGC37_ALL_SAME_COUNT",
                                            "NGC37_ALL_SAME_FREQ"
                        ]
    elif re.search("37",genome_ref):
         ngc_header = ["NGC37_ALL_SAME_COUNT","NGC37_ALL_SAME_FREQ",
                       "NGC37_FAM_SAME_ID","NGC37_FAM_SAME_COUNT",
                       "PROBAND","MOTHER","FATHER","SIB"
                      ]

    #fam_header = ["NGC38_FAM_1K_SAME_ID","PROBAND_FAM_1K","MOTHER_FAM_1K","FATHER_FAM_1K"]
    decipher_header = ["DEC_SAME_ID","DEC_SAME_ROF_ID"]
    dbvar_header    = ["DBVAR_SAME_ID","DBVAR_SAME_ROF_ID"]
    
    ens_trans_header = ['ENSMBLT_ID','ENS_HGNC_ID','ENS_GENE_BIOTYPE',
                        'ENS_TRANSCRIPT_ID','ENS_TRANSCRIPT_BIOTYPE',
                        'ENS_IN_GENELIST','ENS_IN_HPOGENE','ENS_IN_IMPRINTED',
                        'ENS_IS_PROT_CODING'
                       ]
    ens_exon_header = ['ENSMBLE_ID','ENS_NUM_OVERLAP_GEXONS']
   
    ref_gene_header = ['REF_GENE_HGNC_ID','REF_GENE_HGNC_SYN','REF_GENE_BIOTYPE',
                       'REF_IN_GENE_LIST','REF_IN_HPOGENE','REF_IN_IMPRINTED',
                       'REF_IS_PROT_CODING'
                      ]
    ref_exon_header = ['REF_GEXON_HGNC_ID','REF_NUM_OVERLAP_GEXONS']
    
    user_prom_header = ['USER_PROMOTER_ID','IN_PROMOTER']

    extra_columns   = ['IN_IMPRINTED_COLL','IN_GENELIST_COLL','IN_HPOGENE_COLL',
                        'IS_PROT_CODING','NUM_OVLAP_EXONS'
                      ]

    if re.search("38",genome_ref):
        print >>wh,"\t".join(rest_header+gn_header+ngclov_header+ngc_header+
                              decipher_header+dbvar_header+ens_trans_header+
                            ens_exon_header+ref_gene_header+ref_exon_header+
                                              user_prom_header+extra_columns
                            )
    elif re.search("37",genome_ref):
        print >>wh,"\t".join(rest_header+gn_header+ngc_header+decipher_header+
                                dbvar_header+ens_trans_header+ens_exon_header+
                                              ref_gene_header+ref_exon_header+
                                                user_prom_header+extra_columns
                            )

    # Extracting and concatenating information from each of the overlap files
    print 'begin rest'
    for i,keys in enumerate(rest_hash):
        key_strs = re.split("\+",keys)
       
        if re.search("GAIN",keys):
            sv_type = ["DUP"]
        elif re.search("LOSS",keys):
            sv_type = ["DEL"]
        else:
            sv_type = [key_strs[3].strip()]
        
        rest_strs = rest_hash[keys]
        st_pos = int(key_strs[1])
        en_pos = int(key_strs[2])

        if re.search('37',genome_ref):
            rest_strs.insert(7,'.')
        
        try:
            svlen = abs(int(rest_strs[6]))
        except:
            svlen = en_pos - st_pos 
        try:
            cnvlen = abs(int(rest_strs[7]))
            #cnvlen = svlen
        except:
            cnvlen = en_pos - st_pos
    
        svlen_all = ''
       
        svlen_all = str((svlen+cnvlen)/2)
        #rest_strs.insert(7,str(cnvlen))
        rest_strs.insert(8,svlen_all)

        rest_list = ["\t".join(rest_strs)]
        gn_list   = ["\t".join(gn_hash[keys])]
        ngc_list  = ["\t".join(ngc_hash[keys])]

        # Extract NGC-37 lift-over information
        if re.search("38",genome_ref):
            ngclov_list   = ["\t".join(ngclov_hash[keys])]

        # Extract - Decipher-DDD information
        dec_list  = ["\t".join(decipher_hash[keys])]
        
        # Extract Ensembl based information - Albas' code
        if ensembl_hash:
            en_list   = ["\t".join(ensembl_hash[keys])]
        
        # Extract Ensembl and Refseq - gene and exon information
        ens_trans_list = ['\t'.join(ens_trans_hash[keys])]
        ens_exon_list  = ['\t'.join(ens_exon_hash[keys])]
        ref_gene_list  = ['\t'.join(ref_gene_hash[keys])]
        ref_exon_list  = ['\t'.join(ref_exon_hash[keys])]
        user_prom_list = ['\t'.join(user_prom_hash[keys])]

        ens_in_gene = any([json.loads(x.lower()) for x in re.split('\/',
                                          ens_trans_hash[keys][5])
                          ]
                         )
        ens_in_hpo  = any([json.loads(x.lower()) for x in re.split('\/',
                                          ens_trans_hash[keys][6])
                         ]
                        )
        
        ens_in_imp  = any([json.loads(x.lower()) for x in re.split('\/',
                                          ens_trans_hash[keys][7])
                         ]
                        )


        ens_in_prot = any([json.loads(x.lower()) for x in re.split('\/',
                                          ens_trans_hash[keys][8])
                          ]
                         )

        ref_in_gene = any([json.loads(x.lower()) for x in re.split('\/',
                                          ref_gene_hash[keys][3])
                          ]
                         )

        ref_in_hpo  = any([json.loads(x.lower()) for x in re.split('\/',
                                          ref_gene_hash[keys][4])
                         ]
                        )
        
        ref_in_imp  = any([json.loads(x.lower()) for x in re.split('\/',
                                          ref_gene_hash[keys][5])
                         ]
                        )


        ref_in_prot = any([json.loads(x.lower()) for x in re.split('\/',
                                          ref_gene_hash[keys][6])
                          ]
                        )
        
        # Extract and Combine Ensembl & RefSeq for genelist_coll, hpo_genelist,
        # presence in protein coding list

        ens_ref_in_imp_coll      = ens_in_imp or ref_in_imp
        ens_ref_in_genelist_coll = ens_in_gene or ref_in_gene
        ens_ref_in_hpo_coll      = ens_in_hpo  or ref_in_hpo
        ens_ref_in_prot          = ens_in_prot or ref_in_prot
        ens_ref_num_ovlap_exons  = max(int(ens_exon_hash[keys][1]),
                                       int(ref_exon_hash[keys][1])
                                      )
        
        # Extract dbVar information
        dbvar_list   = ["\t".join(dbvar_hash[keys])]
        
        # Concatenate the information based on reference genome build
        if re.search("38",genome_ref):
            merge_list = rest_list+sv_type+gn_list+ngclov_list+\
                         ngc_list+dec_list+dbvar_list+ens_trans_list+\
                         ens_exon_list+ref_gene_list+ref_exon_list+\
                         user_prom_list+\
                         [str(ens_ref_in_imp_coll).upper(),
                          str(ens_ref_in_genelist_coll).upper(),
                          str(ens_ref_in_hpo_coll).upper(),
                          str(ens_ref_in_prot).upper(),
                          str(ens_ref_num_ovlap_exons)
                         ]

        elif re.search("37",genome_ref):
            merge_list = rest_list+sv_type+gn_list+ngc_list+\
                         dec_list+dbvar_list+ens_trans_list+\
                         ens_exon_list+ref_gene_list+ref_exon_list+\
                         user_prom_list+\
                         [str(ens_ref_in_imp_coll).upper(),
                          str(ens_ref_in_genelist_coll).upper(),
                          str(ens_ref_in_hpo_coll).upper(),
                          str(ens_ref_in_prot).upper(),
                          str(ens_ref_num_ovlap_exons)
                         ]
 
        key_strs = re.split('\+',keys)[0:5]
        #print >>wh,keys+"\t"+"\t".join(merge_list)
        print >>wh,keys+'\t'+'\t'.join(key_strs)+"\t"+"\t".join(merge_list)

    wh.close()

    
