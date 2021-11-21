#loc3_5 = "/home/aak64/R/x86_64-pc-linux-gnu-library/3.5/"

#require("data.table",lib=loc3_5)
#require("plyr",lib=loc3_5)
#require("future.apply",lib=loc3_5)
#require('R.utils',lib=loc3_5)

require("data.table")
require("plyr")
require("future.apply")
require('R.utils')


`%ni%` <- Negate(`%in%`)

initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, 
                                                        initial.options)])
script.basename <- dirname(script.name)

source(file.path(script.basename,"/familyFilterFunctions.R"))
#source("/home/aak64/rds/hpc-work/REPOS/ngc_sv37/familyFilterFunctions.R")
#suppressMessages(require(optparse,lib=loc3_5))
suppressMessages(require(optparse))

#############################################################################
#                                                                           #
#  Formulating command line arguments and options for processing the data   #
#                                                                           #
#############################################################################

option_list = list(
            make_option(c("-v", "--vars"), type="character", default="NULL", 
                        help="Input variants file with SV coordinates", 
                        metavar="character"),
            make_option(c("-o", "--outDir"), type="character",default="NULL",
                        help="Output File with annotated Genes",
                        metavar="character"),
            make_option(c("-f","--family"), type="character", default=NULL, 
                        help="family id",metavar="character"),
            make_option(c("-n","--affected"), type="character", default=NULL, 
                        help="affected proband",metavar="character"), 
            make_option(c("-p","--ped"), type="character", default=NULL, 
                        help="pedigree file",metavar="character"),
            #make_option(c("-s","--filterFunc"),type="character", default=NULL,
            #            help="filtering subroutines",metavar="character")
            make_option(c("-c","--ovFrac"),type="character",default=NULL,
                        help="overlap fraction list",metavar="character"),
            make_option(c("-i","--imprint"),type="character",default=NULL,
                        help="imprinted gene list",metavar="character"),
            make_option(c("-r","--reference"),type="character",default=NULL,
                        help="reference genome build (grch37/grch38)",metavar="character"),
            make_option(c("-a","--aff"),type="character",default=NULL,
                        help="[mother|father], will not do regular filtering",
                        metavar="character") 
                );

opt = parse_args(OptionParser(option_list=option_list))

inpFile = opt$vars
outDir = opt$outDir
fam  = opt$family
ngc_id = opt$affected
ped_file = opt$ped
#ovp_list = unlist(strsplit(opt$ovpFrac,","))
ovFrac = opt$ovFrac
imprintedGenesFile = opt$imprint
ref_genome = opt$reference

#system(paste('cp -r ',outDir,' ',outDir,'_V1',sep=''))
#fam = "NGC00296"; ovFrac_list = c("0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0")
#fam = "NGC00296"; ovFrac_list = c("0.7")
#fam = "NGC00315"; ovFrac_list = c("1.0")
#ngc_id = paste(fam,"_01",sep="")

#for (ovFrac in ovFrac_list) {

    cat("Processing for overlap fraction:",ovFrac,"\n")

    #inpFile = paste("/home/aak64/rds/rds-flr24-wgs10k-ngc/WGS10K/data/NGC/us/misc/TEST_SV/20200630/tmp_data/",fam,"/",ngc_id,"/",ngc_id,".merged.all.",ovFrac,".overlap.fmt.bed.gz",sep="")
    #outDir = paste("/home/aak64/rds/rds-flr24-wgs10k-ngc/WGS10K/data/NGC/us/misc/TEST_SV/20190916/tmp_data/",fam,"/",ngc_id,"/filter/",sep="")
    outFile = paste(outDir,ngc_id,".merged.all.",ovFrac,".overlap.freqImpFilter.tsv",sep="")
    outSVType = paste(outDir,ngc_id,".merged.all.",ovFrac,".overlap.SVTypes.tsv",sep="")
    outSVSumm = paste(outDir,ngc_id,".merged.all.",ovFrac,".overlap.SVsummary.tsv",sep="")

    #system(paste("mkdir -p ",outDir,sep=""))

    #outFile = "/home/aak64/rds/rds-flr24-wgs10k-ngc/WGS10K/data/NGC/us/misc/TEST_SV/20190715/tmp_data/NGC00324/NGC00324_01/NGC00324_01.merged.all.overlap.freqImpFilter.tsv"

    #ped_file = "/home/aak64/rds/rds-flr24-wgs10k-ngc/WGS10K/data/NGC/us/release/20190614/manifest/manifest.txt"

    ############################################################
    #                                                          # 
    #  Processing family structure                              # 
    #                                                          #
    ############################################################

    ped = fread(ped_file,header=T,stringsAsFactors=F)
    ped = data.frame(ped)

    ped_fam <- subset(ped, family_id == fam)
    samples <- ped_fam$ilmn_id
    family_structure <- length(samples)

    affected <- subset(ped_fam, affected == 2)$ngc_id
    non_affected <- subset(ped_fam, affected == 1)$ngc_id

    mother <- unique(subset(ped_fam, affected == 2)$mother_id)
    father <- unique(subset(ped_fam, affected == 2)$father_id)

    imprintedgenes = read.table(imprintedGenesFile,sep='\t',head=T)

    if(mother != "None" & father != "None"){
            parents <- c(mother, father)
    }else if(mother != "None"){
            parents <- c(mother)
    }else if(father != "None"){
            parents <- c(father)
    }

#############################################
#                                           #
# Variants filtering by quality and impact  #
#                                           #
#############################################

    vars = fread(inpFile,sep="\t",head=T,stringsAsFactor=F)
    #vars = fread(inpFile,sep="\t",stringsAsFactor=F)
    vars = as.data.frame(vars)
    #sv_len = abs(as.numeric(vars[,'SVLEN'])); sv_len[is.na(sv_len)] <- 0
    #vars$SVLEN_ALL <- as.numeric(sv_len)

    #save.image(paste(outDir,'.RData',sep=''))
    
    vars_qual = qualFilter(vars,ngc_id,ref_genome)
    vars_qual_imp = impFilter(vars_qual)
    vars_qual_imp_freq = recessive_filter(vars_qual_imp,ref_genome)
    vars_qual_imp1 = vars_qual_imp[apply(data.frame(vars_qual_imp[,'PROBAND'])
                                         , MARGIN=1
                                         , function(x) {
                                            checkIn(x,'^0/1:|^1/1:|^1/0:|^1:','-') 
                                           }
                                        )
                                 ,]
    #write.table(vars_qual_imp_freq,file=outFile,sep="\t",row.names=F,col.names=T,quote=F)
    vars_qual_imp$MOI <- NA

    affected = affected[affected %in% ngc_id]

    if(family_structure %in% c(1,2,3,4) & length(affected) == 1){

        affected <- subset(ped_fam, affected == 2)$ngc_id
        non_affected_ped <- subset(ped_fam, affected == 1)$ngc_id
        mother <- unique(subset(ped_fam, affected == 2)$mother_id)
        father <- unique(subset(ped_fam, affected == 2)$father_id)
        non_affected = c()

        for (ele in non_affected_ped) {
            if (ele==mother) {
                non_affected = c(non_affected,'MOTHER')
            }else if (ele == father){
                non_affected = c(non_affected,'FATHER')
            }else{
                non_affected = c(non_affected,'SIB')
            }
        }

        parents = c()
        if(mother != "None" & father != "None"){
            parents <- c(mother, father)
            parents <- c('MOTHER','FATHER')
            mother = 'MOTHER'
            father = 'FATHER'
        }else if(mother != "None"){
            parents <- c(mother)
            mother = 'MOTHER'
            parents = mother
        }else if(father != "None"){
            parents <- c(father)
            father = 'FATHER'
            parents = father
        }

        if(family_structure %in% c(1,2)){
                #vars_filt <- subset(vars_qual_imp1, IN_GENELIST_COLL |
                #                                    IN_HPOGENE_COLL
                #                    )
                vars_filt <- vars_qual_imp1
        }else{
                vars_filt <- vars_qual_imp1
        }
    
        
        vars_filt = as.data.frame(vars_filt)
        #save.image(paste(outDir,'.RData',sep=''))

        dnv = denovo(vars_filt,affected, non_affected,mother,father,
                                         imprintedgenes,'ensembl',ref_genome)
        ar_hom = AR_hom(vars_filt,affected, non_affected, parents,ref_genome)
        xl = XL(vars_filt,affected, non_affected,'ensembl',ref_genome)
        ar_c = AR_comphet(as.data.frame(vars_filt),affected, non_affected, 
                                    parents,mother,father,'ensembl',ref_genome)
        #d_c = denovo_comphet(vars_filt,non_affected)
        ih = inherited(subset(vars_filt,IN_HPOGENE_COLL),affected,non_affected,
                              parents,ref_genome
                      )
        
        message(" --Combining all the output filters")
        #outData = rbind(dnv,ar_hom,xl,ar_c,d_c,ih)
        outData = rbind(dnv,ar_hom,xl,ar_c,ih)
    
        # Removal of NGC37_ALL_SAME_ID, ENSEMBLT_ID,ENS_TRANSCRIPT_ID, ENSEMBLE_ID
        sp_col_names = c('NGC37_ALL_SAME_ID','ENSMBLT_ID','ENS_TRANSCRIPT_ID','ENSMBLE_ID')
        sp_col_ind = which(colnames(outData) %in% sp_col_names)
        outSummary = outData[,-c(sp_col_ind)]

        write.table(outSummary,sep="\t",file=outSVSumm,col.names=T,row.names=F,quote=F)
        write.table(outData,sep="\t",file=outFile,col.names=T,row.names=F,quote=F)

        if (file.exists(outSVType)){
            system(paste("rm ",outSVType,sep=""))
        }
        
        message("Denovo")
        #cat("Denovo\n",file=outSVType,append=T)
        a1 = as.data.frame(table(dnv[,'SVTYPE']))
        #if(dim(a1)[1]!=0){ 
        #    a2 = rep('Denovo',dim(a1)[1])
        #}else{
        #    a1 = as.data.frame(0)
        #    a2 = rep('Denovo',1)
        #}
        a2 = rep('Denovo',dim(a1)[1])
        a3 = cbind(a2,a1)
        #write.table(table(dnv[,'SVTYPE']),file=outSVType,append=T,row.names=F,col.names=F,quote=F)
        write.table(a3,file=outSVType,append=T,row.names=F,col.names=F,quote=F)
        #cat("\n",file=outSVType,append=T)
        
        message("AR-hom")
        #cat("AR-homozygous\n",file=outSVType,append=T)
        a1 = as.data.frame(table(ar_hom[,'SVTYPE']))
        a2 = rep('AR-Hom',dim(a1)[1])
        a3 = cbind(a2,a1)
        write.table(a3,file=outSVType,append=T,row.names=F,col.names=F,quote=F)
        #cat("\n",file=outSVType,append=T)

        #message("X-linked")
        #cat("X-linked\n",file=outSVType,append=T)
        a1 = as.data.frame(table(xl[,'SVTYPE']))
        a2 = rep('X-linked',dim(a1)[1])
        a3 = cbind(a2,a1)
        write.table(a3,file=outSVType,append=T,row.names=F,col.names=F,quote=F)
        #cat("\n",file=outSVType,append=T)

        message("AR-comphet")
        #cat("AR-comphet\n",file=outSVType,append=T)
        a1 = as.data.frame(table(ar_c[,'SVTYPE']))
        a2 = rep('AR-comphet',dim(a1)[1])
        a3 = cbind(a2,a1)
        write.table(a3,file=outSVType,append=T,row.names=F,col.names=F,quote=F)
        #cat("\n",file=outSVType,append=T)

        #message("Denovo-comphet")
        #cat("Denovo comphet\n",file=outSVType,append=T)
        #write.table(table(d_c[,'SVTYPE']),file=outSVType,append=T,row.names=F,col.names=F,quote=F)
        #cat("\n",file=outSVType,append=T)

        message("Inherited")
        #cat("Inherited\n",file=outSVType,append=T)
        a1 = as.data.frame(table(ih[,'SVTYPE']))
        a2 = rep('Inherited',dim(a1)[1])
        a3 = cbind(a2,a1)
        write.table(a3,file=outSVType,append=T,row.names=F,col.names=F,quote=F)
        #cat("\n",file=outSVType,append=T)
    }
#}
