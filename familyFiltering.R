require("data.table")
require("plyr")
require("future")
require("future.apply")
require('optparse')
require('R.methodsS3')
require('R.oo')
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
ovFrac = opt$ovFrac
imprintedGenesFile = opt$imprint
ref_genome = opt$reference


    cat("Processing for overlap fraction:",ovFrac,"\n")

    outFile = paste(outDir,ngc_id,".merged.all.",ovFrac,".overlap.freqImpFilter.tsv",sep="")
    outSVType = paste(outDir,ngc_id,".merged.all.",ovFrac,".overlap.SVTypes.tsv",sep="")
    outSVSumm = paste(outDir,ngc_id,".merged.all.",ovFrac,".overlap.SVsummary.tsv",sep="")

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

    vars = fread(inpFile,sep="\t",stringsAsFactor=F,header=TRUE,fill=T)
    vars = as.data.frame(vars)

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
                vars_filt <- vars_qual_imp1
        }else{
                vars_filt <- vars_qual_imp1
        }
    
        
        vars_filt = as.data.frame(vars_filt)

        dnv = denovo(vars_filt,affected, non_affected,mother,father,
                                         imprintedgenes,'ensembl',ref_genome)
        ar_hom = AR_hom(vars_filt,affected, non_affected, parents,ref_genome)
        xl = XL(vars_filt,affected, non_affected,'ensembl',ref_genome)
        ar_c = AR_comphet(as.data.frame(vars_filt),affected, non_affected, 
                                    parents,mother,father,'ensembl',ref_genome)
        ih = inherited(subset(vars_filt,IN_HPOGENE_COLL),affected,non_affected,
                              parents,ref_genome
                      )
        
        message(" --Combining all the output filters")
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
        a1 = as.data.frame(table(dnv[,'SVTYPE']))
        a2 = rep('Denovo',dim(a1)[1])
        a3 = cbind(a2,a1)
        write.table(a3,file=outSVType,append=T,row.names=F,col.names=F,quote=F)
        #cat("\n",file=outSVType,append=T)
        
        message("AR-hom")
        a1 = as.data.frame(table(ar_hom[,'SVTYPE']))
        a2 = rep('AR-Hom',dim(a1)[1])
        a3 = cbind(a2,a1)
        write.table(a3,file=outSVType,append=T,row.names=F,col.names=F,quote=F)
        #cat("\n",file=outSVType,append=T)

        message("X-linked")
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

        message("Inherited")
        #cat("Inherited\n",file=outSVType,append=T)
        a1 = as.data.frame(table(ih[,'SVTYPE']))
        a2 = rep('Inherited',dim(a1)[1])
        a3 = cbind(a2,a1)
        write.table(a3,file=outSVType,append=T,row.names=F,col.names=F,quote=F)
        #cat("\n",file=outSVType,append=T)
    }
#}
