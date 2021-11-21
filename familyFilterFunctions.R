#loc3_5 = "/home/aak64/R/x86_64-pc-linux-gnu-library/3.5/"

#require("data.table",lib=loc3_5)
#require("plyr",lib=loc3_5)
#require("future.apply",lib=loc3_5)

require("data.table")
require("plyr")
require("future.apply")

`%ni%` <- Negate(`%in%`)



#######################################

qualFilter = function(df,ngc_id,ref_genome) {

    message(' -- Processing quality filter')
    if (nrow(df) >=1) {
        #tmpA <- subset(df,grepl("Manta",NGC38_ALL_SAME_ID) & grepl("Canvas",NGC38_ALL_SAME_ID))
        if(ref_genome=='grch37'){
            tmpA <- subset(df,grepl("Manta",NGC37_FAM_SAME_ID) & grepl("Canvas",NGC37_FAM_SAME_ID))
            tmpA1 <- subset(tmpA,grepl(ngc_id,NGC37_FAM_SAME_ID))
        }
        else if(ref_genome=='grch38'){
            tmpA <- subset(df,grepl("Manta",NGC38_FAM_SAME_ID) & grepl("Canvas",NGC38_FAM_SAME_ID))
            tmpA1 <- subset(tmpA,grepl(ngc_id,NGC38_FAM_SAME_ID))
        }
        tmpB <- subset(df, SVTYPE %in% c('DEL','DUP','INV','INS','BND') & grepl("PASS|^MGE10kb$|^CLT10kb$",FILTER)) #Include MGE10KB, CLT10Kb
        #tmpC <- subset(df, SVTYPE %in% c('INV') & FILTER == 'PASS' )#& abs(as.numeric(SVLEN_ALL)) < 500000)
        #tmpD <- subset(df,SVTYPE %in% c('DEL','BND') & FILTER == 'PASS')

        #x = unique(rbind(tmpA,tmpB,tmpC,tmpD))
        x = unique(rbind(tmpA,tmpA1,tmpB))
    }
    else {
        x <- NULL
    }
    return(x)
}

impFilter <- function(df) {
  
    print(dim(df))
    message(' -- Processing impact filter')
    if (nrow(df) >=1) {
        tmpA <- subset(df, SVTYPE %in% c('DEL','DUP','INV','INS') & IS_PROT_CODING & ((as.numeric(NUM_OVLAP_EXONS) >0)| SVLEN_ALL > 1000 | IN_PROMOTER))
        tmpA1 <- subset(df, SVTYPE %in% c('DEL','DUP') & ((DEC_SAME_ID != 'NA' & SVLEN>=500000)|DEC_SAME_ROF_ID !='NA'))
        tmpA2 <- subset(df, SVTYPE %in% c('DEL','DUP') & ((grepl("pathogenic",DBVAR_SAME_ROF_ID) & SVLEN>=500000)))
        #tmpA2 <- subset(df, SVTYPE %in% c('DEL','DUP') & (grepl("pathogenic",DBVAR_SAME_ROF_ID))

        #tmpB = subset(df, SVTYPE %in% c('INV','INS') & IN_GENELIST_COLL & IS_PROT_CODING & ((as.numeric(NUM_OVLAP_EXONS) > 0) | SVLEN_ALL >1000)) #Get rid of Gene list and merge with tmpA
        #tmpA1 = subset(df,SVTYPE %in% c('DEL','DUP'))
        #tmpC = subset(df,SVTYPE %in% c('INS') & IN_GENELIST_COLL & IS_PROT_CODING & SVLEN_ALL >1000) #Get rid of this; make INS and BND same
        tmpC = subset(df,SVTYPE %in% c('BND') & IS_PROT_CODING ) #Get rid of gene list; Only BND
        #tmpD = subset(df,SVTYPE %in% c('DEL','DUP','INS','INV') & IN_PROMOTER & IS_PROT_CODING)
        #tmpD = subset(df,SVTYPE %in% c('BND') & IS_PROT_CODING)
        
        #x = unique(rbind(tmpA,tmpA1,tmpA2,tmpB,tmpC,tmpD))
        #x = unique(rbind(tmpA,tmpC,tmpD))
        x = unique(rbind(tmpA,tmpA1,tmpA2,tmpC))
    }else {
        x <- NULL
    }
    return(x)
}

red_pen_filter <- function(df,ref_genome) {
   
    message(' -- Processing penetrance filter')
    if(nrow(df) >=1) {
        
        if(ref_genome=='grch37'){
            x = subset(df, ((is.na(GN_AF)|GN_AF <=0.001)
                        & (is.na(NGC37_ALL_SAME_FREQ) | NGC37_ALL_SAME_FREQ <=0.10)
                        #& (is.na(NGC37_ALL_SAME_FREQ) | NGC37_ALL_SAME_FREQ <=0.01)
                       )
                   )
        }else if(ref_genome=='grch38'){
            x = subset(df, ((is.na(GN_AF)|GN_AF <=0.001)
                        & (is.na(NGC38_ALL_SAME_FREQ) | NGC38_ALL_SAME_FREQ <=0.10)
                        #& (is.na(NGC37_ALL_SAME_FREQ) | NGC37_ALL_SAME_FREQ <=0.01)
                       )
                   )
        }
    }else {
        x <- df
    }
    return(x)
}

dominant_inheritance_filter <- function(df) {
    
    if(nrow(df) >=1) {
        # Filter the variants based on OMIM reported genes
        tmp1 <- subset(df, is.na(df$OMIM_gene_moi)
                       | df$OMIM_gene_moi == "NA"
                       | df$OMIM_gene_moi == ""
                       )
        # Filter the variants based on OMIM reported genes and mode of inheritance
        tmp2 = df[apply(data.frame(df[,'OMIM_gene_moi'])
                        , MARGIN = 1
                        , function(x) {
                            FALSE %in% paste(
                                             unlist(strsplit(x," | ",fixed = TRUE))
                                             %in% c("Autosomal recessive",
                                                    "X-linked recessive")
                                            )
                                      }       
                        )
                  ,]
        x = unique(rbind(tmp1, tmp2))
    }else {
        x = df
    }
}


recessive_filter <- function(df,ref_genome) {

    message(' -- Processing recessive filter')
    if(nrow(df) >=1) {
        if(ref_genome=='grch37'){
            x = subset(df,((is.na(GN_AF)|GN_AF <=0.01) 
                   #& (is.na(NGC38_ALL_SAME_FREQ) | NGC38_ALL_SAME_FREQ <=0.20) 
                   & (is.na(NGC37_ALL_SAME_FREQ) | NGC37_ALL_SAME_FREQ <=0.10) 
                   #& (is.na(NGC37_ALL_SAME_FREQ) | NGC37_ALL_SAME_FREQ <=0.03)
                  )
                )
        }else if(ref_genome=='grch38'){
            x = subset(df,((is.na(GN_AF)|GN_AF <=0.01) 
                   #& (is.na(NGC38_ALL_SAME_FREQ) | NGC38_ALL_SAME_FREQ <=0.20) 
                   & (is.na(NGC38_ALL_SAME_FREQ) | NGC38_ALL_SAME_FREQ <=0.10) 
                   #& (is.na(NGC37_ALL_SAME_FREQ) | NGC37_ALL_SAME_FREQ <=0.03)
                  )
                )
        }
    }else {
        x <- df
    }
    return(x)
}

dominant_filter <- function(df,ref_genome){

    if(ref_genome == 'grch37'){
        x = subset(df,(
                   (is.na(GN_AF)| GN_AF <=0.001)
                    #x = subset(df,((is.na(GN_AF)| GN_AC <= 200)
                    & (is.na(NGC37_ALL_SAME_COUNT) | NGC37_ALL_SAME_COUNT <=5) 
                    #& (is.na(NGC37_ALL_SAME_COUNT) | NGC37_ALL_SAME_COUNT <=5)
                  )
              )
    }else if(ref_genome=='grch38'){
        x = subset(df,(
                   (is.na(GN_AF)| GN_AF <=0.001)
                    #x = subset(df,((is.na(GN_AF)| GN_AC <= 200)
                    & (is.na(NGC38_ALL_SAME_COUNT) | NGC38_ALL_SAME_COUNT <=5) 
                    #& (is.na(NGC37_ALL_SAME_COUNT) | NGC37_ALL_SAME_COUNT <=5)
                  )
              )
    }
    return(x)
}

checkNotIn <- function(x,g,s) {
    y = c()
    for(i in x) {
        i = as.character(i)
        y = append(y,all(!grepl(g,paste0(strsplit(i,s)[[1]]),perl=TRUE)))
    }
    return(y)
}

checkIn <- function(x,g,s) {
    y = c()
    for(i in x) {
        i = as.character(i)
        y = append(y,all(grepl(g,paste0(strsplit(i,s)[[1]]),perl=TRUE)))
    }
    return(y)
}

denovo <- function(df,affected,non_affected,mother,father,imprintedgenes,
                                                        gene_type,ref_genome) {
    
        message(" -- Applying Denovo filter")
        if (gene_type=='ensembl') {
            gene_col = which(colnames(df)=='ENS_HGNC_ID')
        }else if (gene_type=='refseq'){
            gene_col = which(colnames(df)=='REF_GENE_HGNC_ID')
        }

        df = as.data.frame(df)
        #Singleton; Affected is hetrozygous   
        if (length(non_affected)==0){
            #tmp2 = df[apply(data.frame(df[,affected])
            tmp2 = df[apply(data.frame(df[,'PROBAND'])
                           , MARGIN = 1
                           , function(x) all(checkIn(x,'^0/1:|^1/0:','-'))
                          )
                    ,]
        }else{
          
            #Affected is hetrozygous, homozygous or hemizygous 
            #tmp1 = df[apply(data.frame(df[,affected])
            tmp1 = df[apply(data.frame(df[,'PROBAND'])
                           , MARGIN = 1
                           , function(x) all(checkIn(x,'^0/1:|^1/0:|^1:|^1/1:|^[.]/1:','-'))
                          )
                    ,] 
            
            tmp1 = as.data.frame(tmp1)
            #tmp2 = tmp1[apply(data.frame(tmp1[,non_affected])
            tmp2 = tmp1[apply(data.frame(tmp1[,non_affected])
            #tmp2 = tmp1[apply(data.frame(tmp1[,which(colnames(tmp1) %in% non_affected==TRUE)])
                            , MARGIN = 1
                            , function(x) all(checkIn(x,"^0/0:|^0:|^[.]/[.]:|^NA","-"))
                            )
                      ,]
            }
            tmp3 = dominant_filter(tmp2,ref_genome)

            tmp4 = imprinted(subset(df,IN_IMPRINTED_COLL),affected,non_affected,
                                mother,father,imprintedgenes,gene_type,ref_genome
                            )
           
            
            x = unique(rbind(tmp3,tmp4))

            if (nrow(x) >=1) {x$MOI <- "denovo"}
            return(x)
}

imprinted <- function(df,affected,non_affected,mother,father,imprintedgenes,
                                                            gene_type,ref_genome){
    
        message(" -- Applying Imprinted filter")

        if (gene_type=='ensembl') {
            #gene_col = which(colnames(df)=='ENS_HGNC_ID')
            gene_col = 'ENS_HGNC_ID'
        }else if (gene_type=='refseq'){
            #gene_col = which(colnames(df)=='REF_GENE_HGNC_ID')
            gene_col = 'REF_GENE_HGNC_ID'
        }

        tmp = df[apply(data.frame(df[,'PROBAND'])
                     , MARGIN = 1
                     , function(x) all(checkIn(x,'^0/1:|^1/0:|^1/1:','-'))
                     )
                 ,]

        tmp2 <- tmp[apply(data.frame(tmp[,non_affected])
                        , MARGIN=1
                        , function(x) all(checkNotIn(x,'^1/1:','-'))
                        )
                  ,]

        if(nrow(tmp2) ==0){
            return(tmp2)
        }
       
        if(mother == "None" & father == "None"){
                message('--singleton')
                tmp3 = tmp2
        
        }else if(mother =="None" | father == "None"){
            if(mother !="None"){
                tmp3 = tmp2[apply(data.frame(tmp2[,c(gene_col,mother)])
                                  , MARGIN = 1
                                  , function(x) checkImprint(x,'mother',imprintedgenes))
                           ,]
            }else if(father !="None"){
                tmp3 = tmp2[apply(data.frame(tmp2[,c(gene_col,father)])
                                  , MARGIN = 1
                                  , function(x) checkImprint(x,'father',imprintedgenes))
                           ,]            
            }
        }else { #Trios
                print('Inside trios')
                tmp3 = tmp2[apply(data.frame(tmp2[,c(gene_col,parents)])
                                  , MARGIN = 1
                                  , function(x) checkImprint(x,'trio',imprintedgenes))
                           ,]
            }
        
        y = red_pen_filter(tmp3,ref_genome)
        return(y)
}

checkImprint <- function(t,tp,imprintedgenes) {
    #################################################################################
    # Function to check if for given SV overlapping gens is present in Imprinted genes
    #
    # Input: [sv-overlap-genes,'mother,father'],'singleton/duo/trio-mode',
    #                                           'imprinted-genelist'
    # Output: boolean true/false
    #################################################################################

    # Get list of imprinted genes
    #g <- imprintedgenes[imprintedgenes$GENE == t[1],]$EXPRESSION
    g <- imprintedgenes[imprintedgenes$GENE %in% 
                                        strsplit(t[1],'/')[[1]]
                        ,]$EXPRESSION
   
    if(tp =='trio'){ #Trio
        g_list = c()
        for (exp in g) {
            if(exp =='Both'){ #Check in both parents
                #return(any(checkIn(t[c(2,3)],'^0/1:|^1/0:','-')))
                g_list = c(g_list,any(checkIn(t[c(2,3)],'^0/1:|^1/0:','-')))
        
            }else if(exp == 'Maternal'){ #check only in Mother
                #return(any(checkIn(t[2],'^0/1:|^1/0:','-')))
                g_list = c(g_list,any(checkIn(t[2],'^0/1:|^1/0:','-')))

            }else if(exp == 'Paternal'){ #Check only in Father
                #return(any(checkIn(t[3],'^0/1:|^1/0:','-')))
                g_list = c(g_list,any(checkIn(t[3],'^0/1:|^1/0:','-')))
            }
        }
        if(length(g_list)!=0){
            return(any(g_list)) 
        }else{
            return(FALSE)
        }
    
    }else if(tp == 'mother'){
        g_list = c()
        for (exp in g) {
            if(exp == 'Maternal'){
                #return(any(checkIn(t[2],'^0/1:|^1/0:','-')))
                g_list = c(g_list,any(checkIn(t[2],'^0/1:|^1/0:','-')))
            }else {
                g_list = c(g_list,TRUE)
            }
        }
        if(length(g_list)!=0){
            return(any(g_list)) 
        }else{
            return(FALSE)
        }

    }else if (tp == 'father'){
        g_list = c()
        for (exp in g){
            if(exp == 'Paternal'){
                g_list = c(g_list,any(checkIn(t[2],'^0/1:|^1/0:','-')))
            }else {
                g_list = c(g_list,TRUE)
            }
        }
        if(length(g_list)!=0){
            return(any(g_list)) 
        }else{
            return(FALSE)
        }
    }
    
    return(FALSE)
}

AR_hom <- function(df,affected, non_affected,parents,ref_genome) {

        message(" --Applying AR Homozygous Filter")
        
        ## Variant is homozygous in affected
        
        #tmp = df[grepl("1/1:",df[,'GT']),]
        #tmp = df[grepl("1/1:",df[,'PROBAND']),]
        tmp = df[apply(data.frame(df[,'PROBAND'])
                       , MARGIN=1
                       , function(x) all(checkIn(x,'^1/1:','-'))
                       )
                 ,]

        if (length(non_affected)==0) {
            tmp3 = tmp

        }else{
            # Variant is non homozygous in non affected (Parents+ unaffected-Sib) if applicable
            tmp2A = tmp[apply(data.frame(tmp[,non_affected])
                                , MARGIN=1
                                , function(x) all(checkNotIn(x,'^1/1:|^1:','-'))
                              )
                        ,]
            # Parents are only het.
            if(length(parents)==2){
                tmp3 = tmp2A[apply(data.frame(tmp2A[,parents])
                                    , MARGIN = 1
                                    , function(x) any(checkIn(x,'^0/1:|^1/0:','-'))
                                   )
                             ,]
            }else{
                tmp3 = tmp2A
            }
        }

        #tmp2 = tmp[apply(data.frame(tmp[,c('MOTHER_FAM_1K','FATHER_FAM_1K')])
        #tmp2 = tmp[apply(data.frame(tmp[,c('MOTHER','FATHER')])
        #tmp2 = tmp[apply(data.frame(tmp[,non_affected])
        #                 ,MARGIN = 1
        #                 #,function(x) all(checkNotIn(x,c("1/1","1"),"-")))
        #                 ,function(x) all(checkNotIn(x,c("^1/1:|^1:"),"-")))
        #            ,]

        #x = recessive_filter(tmp2)
        #x = recessive_filter(tmp3,ref_genome)
        x = recessive_filter(tmp3,ref_genome)
        
        if (nrow(x) >=1) {x$MOI <- "AR_hom"}
        return(x)
}


XL <- function(df,affected, non_affected,gene_type,ref_genome){
        
        message(" --Applying X-Linked filter")
      
        # Gene type
        if (gene_type=='ensembl') {
                gene_col_name = which(colnames(df)=='ENS_HGNC_ID')
        }else if (gene_type=='refseq'){
                gene_col_name = which(colnames(df)=='REF_GENE_HGNC_ID')
        }

        ## Extract only X chrom
        df.X = df[grepl("^chrX|^X",df[,'SV_ID']),]

        ## Hemizygous variant in X-chrom affected; inherited from mother
        df.hemi.aff = df.X[apply(data.frame(df.X[,'PROBAND'])
                                    , MARGIN = 1
                                    #, function(x) all(checkIn(x,'^1:|^1/1:|0/1:','-'))
                                    , function(x) all(checkIn(x,'^1:','-'))
                                 )
                           ,]
        if(length(non_affected)==0) {
                df.filt = df.hemi.aff    

        }else {
                df.unaff = df.hemi.aff[apply(data.frame(df.hemi.aff[,non_affected])
                                            ,MARGIN = 1
                                            ,function(x) all(
                                                    checkNotIn(x,'^1/1:|^1:','-')
                                                    )
                                            )
                                       ,]
                df.filt = df.unaff
        }

        tmp4 = subset(df.filt,grepl('PCDH19',ENS_HGNC_ID) || grepl('PCDH19',REF_GENE_HGNC_ID))

        ##### To DO : implement hemi filter #####
        ## x = hemi_filter(df.filt)
        #########################################    
        
        #tmp = df[(grepl("chrX",df[,'SV_ID']) & grepl("^1:",df[,'GT'])),]
        #if (nrow(tmp) <1) {return(NULL)}


        #tmp2 = tmp[apply(data.frame(tmp[,c('MOTHER_FAM_1K','FATHER_FAM_1K')])
        #tmp2 = tmp[apply(data.frame(tmp[,c('MOTHER','FATHER')])
        ##tmp2 = tmp[apply(data.frame(tmp[,non_affected])
        ##                ,MARGIN = 1
        ##                ,function(x) all(checkNotIn(x,c("1/1:|1:"),"-")))
        ##            ,]
        
        x = recessive_filter(df.filt,ref_genome)

        if(nrow(x) >=1){x$MOI <- "XLR"}
        return(x)
}

AR_comphet <- function(df,affected,non_affected,parents,mother,father,
                                gene_type='ensembl',ref_genome){

    message(" --Applying AR comphet filter")

    ####################################################################
    # Count those SV-IDS which were called by both Manta and Canvas.   #
    # Keep only the Manta calls                                        #
    #################################################################### 

    if(ref_genome=='grch37'){
        tmpA <- subset(df,!(grepl("Manta",NGC37_FAM_SAME_ID) & grepl("Canvas",NGC37_FAM_SAME_ID)))
        tmpB <- subset(df,(grepl("Manta",NGC37_FAM_SAME_ID) & grepl("Canvas",NGC37_FAM_SAME_ID)))
        tmpC <- subset(tmpB,grepl('Manta',SV_ID))
        df = rbind(tmpA,tmpC)
    }else if(ref_genome=='grch38'){
        tmpA <- subset(df,!(grepl("Manta",NGC38_FAM_SAME_ID) & grepl("Canvas",NGC38_FAM_SAME_ID)))
        tmpB <- subset(df,(grepl("Manta",NGC38_FAM_SAME_ID) & grepl("Canvas",NGC38_FAM_SAME_ID)))
        tmpC <- subset(tmpB,grepl('Manta',SV_ID))
        df = rbind(tmpA,tmpC)
    }

    df.recs = recessive_filter(df,ref_genome)
    affected_hets = df.recs[apply(data.frame(df.recs[,'PROBAND'])
                                    , MARGIN = 1
                                    , function(x) all(checkIn(x,'^0/1:|^1/0:','-'))
                                  )
                            ,]

    if (gene_type=='ensembl') {
        gene_col = which(colnames(df)=='ENS_HGNC_ID')
    }else if (gene_type=='refseq'){
        gene_col = which(colnames(df)=='REF_GENE_HGNC_ID')
    }

    cat('geneCol: ',gene_col,'\n')
    if(length(parents) == 0) { #singleton
        
        message(' -- Working for Singletons')

        #Extract affected - het individual gene list and counts
        ah_genes_count = sort(
                              table(
                                    unlist(
                                       #apply(data.frame(affected_hets[,'SYMBOL'])
                                       apply(data.frame(affected_hets[,gene_col])
                                            , MARGIN=1
                                            , function(x) {
                                                    gene_list = strsplit(x,'/')[[1]];
                                                    return(gene_list)
                                                  }
                                            )
                                          )
                                  ),decreasing=T
                             )

        # Choose only those genes whose counts are > 1 i.e are comphet
        ah_comphet_genes = names(ah_genes_count[ah_genes_count>1]) 
       
        # Extract the variants that are affected-hets and comphets (AH-CH)
        #x = affected_hets[apply(data.frame(affected_hets[,'SYMBOL'])
        x = affected_hets[apply(data.frame(affected_hets[,gene_col])
                                        , MARGIN=1
                                        , function(g) {
                                            gene_list = strsplit(g,'/')[[1]];
                                            vals = any(gene_list %in% 
                                                        ah_comphet_genes
                                                      ) 
                                            return(vals)
                                            }
                                       )
                                 ,]
    }else{ #parents > 0
        message(' -- Working for Duo/Trios')

        # Extract list of non-homozygous variants present in non-affected
        nonaffected_nonhom = affected_hets[apply(
                                            data.frame(affected_hets[,non_affected])
                                              , MARGIN = 1
                                              , function(x) all(
                                                     checkNotIn(x,'^1/1:|^1:','-')
                                                     )
                                                )
                                           ,]

        # Extract list of denovo variants
        denovo = nonaffected_nonhom[apply(
                                      data.frame(nonaffected_nonhom[,non_affected])
                                        , MARGIN=1
                                        , function(x) all(
                                          checkIn(x,'^0/0:|^0:|^[.]/[.]:|^NA','-'))
                                        )
                                    ,]
     
        # Get List of denovo genes
        #denovo_genes = unique(unlist(apply(data.frame(denovo[,'SYMBOL'])
        denovo_genes = unique(unlist(apply(data.frame(denovo[,gene_col])
                                , MARGIN=1
                                , function(x) {
                                              gene_list = strsplit(x,'/')[[1]];
                                              return(gene_list)
                                           }
                                        )
                                    )
                             )

        if(length(parents)==1){ #duo

            # Get count of all nonaffected non homozygous genes
            nanh_genes_count = sort(
                                table(
                                 unlist(
                                   #apply(data.frame(nonaffected_nonhom[,'SYMBOL'])
                                   apply(data.frame(nonaffected_nonhom[,gene_col])
                                            , MARGIN=1
                                            , function(x) {
                                                    gene_list = strsplit(x,'/')[[1]];
                                                    return(gene_list)
                                                  }
                                            )
                                          )
                                      ),decreasing=T
                                   )
            # Choose only those genes whose counts are > 1 i.e are comphet
            nanh_comphet_genes = names(nanh_genes_count[nanh_genes_count>1])
    
            # Extract the variants that are Non-affected non-hom comphets (NANH-CH)
            #tmp2 = nonaffected_nonhom[apply(data.frame(nonaffected_nonhom[,'SYMBOL'])
            tmp2 = nonaffected_nonhom[apply(
                                        data.frame(nonaffected_nonhom[,gene_col])
                                           , MARGIN=1
                                           , function(g) {
                                               gene_list = strsplit(g,'/')[[1]];
                                               vals = any(gene_list %in% 
                                                           nanh_comphet_genes
                                                         ) 
                                                return(vals)
                                              }
                                           )
                                      ,]

            # Overlapping variants between NANH-CH and Denovo-genes
            #x = tmp2[apply(data.frame(tmp2[,'SYMBOL'])
            x = tmp2[apply(data.frame(tmp2[,gene_col])
                                        , MARGIN=1
                                        , function(g) {
                                                gene_list = strsplit(g,'/')[[1]];
                                                vals = any(gene_list %in% 
                                                            denovo_genes
                                                          ) 
                                                return(vals)
                                         }
                          )
                    ,] 

         
            }else { #trio

                message(' -- Processing Trio')
                # Extract variants that overlap with denovo genes from NANH
                nanh_denovo = nonaffected_nonhom[apply(
                                    data.frame(nonaffected_nonhom[,gene_col])
                                            ,MARGIN=1
                                            ,function(g) {
                                                gene_list = strsplit(g,'/')[[1]];
                                                vals = any(gene_list %in%
                                                              denovo_genes
                                                          )
                                                return(vals)
                                                }
                                            )
                                    ,]

                # Extract denovo genes comphet in trio; At least one variant is 
                #                                       het in exactly one parent.
                #                                       AND
                #                                       At leasr one variant is het
                #                                       in other parent
                gene_denovo_comphet = any(
                                          apply(nanh_denovo[,parents]
                                            , MARGIN = 1
                                            , function(x) {
                                              ( any(checkIn(x,'^0/1:|^1/0:','-'))
                                                &
                                                any(checkIn(x,
                                                  '^0/0:|^0:|^[.]/[.]:|^NA','-'
                                                           )
                                                   )
                                              )
                                                  }
                                              )
                                         )

                # Extract variants that overlaps with denovo comphet genes
                tmp3 = nonaffected_nonhom[apply(
                                           data.frame(nonaffected_nonhom[,gene_col])
                                                 , MARGIN=1
                                                 , function(g) {
                                                     gene_list=strsplit(g,'/')[[1]];
                                                     vals = any(gene_list %in%
                                                            gene_denovo_comphet
                                                            )
                                                    }
                                              )
                                         ,]
                tmp4 = AR_comphet_true(nonaffected_nonhom, affected, non_affected, 
                                        parents,mother,father,gene_col,gene_type)
                x <- unique(rbind(tmp3,tmp4))
            } #end trio
        } #end else

    if(nrow(x) >=1){x$MOI = 'AR_comphet'}
    return(x)
} #end function

AR_comphet_true <- function(nonaffected_nonhom, affected, non_affected,
                                parents,mother,father,gene_col,gene_type) {
##Genes with 2 variants, one from each parent
## Reports from autosomal only
## Require both parents
## Doesnot check if CH in unaffected
    message(' -- Inside AR Comphet true function')
    
    df = as.data.frame(nonaffected_nonhom)
    df$motherMask = apply(data.frame(df[,mother])
                            , MARGIN=1
                            , function(x) all(checkIn(x,'^0/1:|^1/0:','-'))
                                     #)
                         )
    df$fatherMask = apply(data.frame(df[,father])
                            , MARGIN=1
                            , function(x) all(checkIn(x,'^0/1:|^1/0:','-'))
                         )

    df$comphetFatherMask = !df$motherMask &  df$fatherMask
    df$comphetMotherMask =  df$motherMask & !df$fatherMask

    #symbolList <- sapply(unique(unlist(apply(data.frame(df[,"SYMBOL"]) 
    symbolList <- sapply(unique(unlist(apply(data.frame(df[,gene_col]) 
                                        ,MARGIN=1
                                        ,function(x) 
                                        strsplit(x, "/")[[1]]
                                        )
                                      )
                               ), function(ANNO_SYMBOL){
                                      #any(df[apply(data.frame(df[,c("SYMBOL")])
                                      any(df[apply(data.frame(df[,c(gene_col)])
                                              , MARGIN = 1
                                              , function(x) { 
                                                  ANNO_SYMBOL %in% strsplit(x[1], "/")[[1]]
                                                }
                                              )
                                            ,]$compHetMotherMask
                                        ) & 
                                      #any(df[apply(data.frame(df[,c("SYMBOL")])
                                      any(df[apply(data.frame(df[,c(gene_col)])
                                                , MARGIN = 1
                                                , function(x) {
                                                        ANNO_SYMBOL %in% 
                                                        strsplit(x[1], "/")[[1]]
                                                  }
                                                )
                                            ,]$compHetFatherMask
                                          ) 
                                  }
                         )


    #df$comphetMask <- apply(data.frame(df[,"SYMBOL"])
    df$comphetMask <- apply(data.frame(df[,gene_col])
                                , MARGIN = 1
                                , function(x) {
                                    any(
                                        strsplit(x, "/")[[1]] 
                                            %in% 
                                        attributes(
                                                   which(symbolList) == TRUE
                                                  )$names
                                       )
                                 }
                             )

    if(length(symbolList[symbolList == TRUE]) > 0){

            x <- df[df$comphetMask,
                    !names(df) %in% c( "fatherMask", "motherMask", 
                                      "compHetFatherMask","compHetMotherMask")
                    ]

    }else {
            y = NULL
    }
    return(y)

} #end function


inherited <- function(df,affected,non_affected,parents,ref_genome) {

    ## Variants where the non affected relatives are all ref
    ## The affected individual are all heterozygous
    ## Reports for autosomal, X-female (and incidentally for male XLR)
    ## Works for duo, trio, quad (both affected have to have variant)
    
    message(" --Applying inherited filter")
    
    # Affected individual are all hetrozygous
    tmp =    df[apply(data.frame(df[,'PROBAND'])
                    , MARGIN = 1
                    , function(x) all(checkIn(x,"^0/1:|^1:|^1/0:",'-'))
                    )
               ,]
    # Atleast one of the parents is het
    if(length(parents)==2) { #trio
        tmp2 = tmp[apply(data.frame(tmp[,parents])
                        , MARGIN = 1
                        , function(x) any(checkIn(x,"^0/1:|^1/0:|^1:",'-'))
                       )
                 ,]
    
    }else{ #singleton/duo
        tmp2 = tmp
    }


    tmp3 = red_pen_filter(tmp2,ref_genome)
    #*** To Do ***#
    ## Implement dominant inheritance filter which requires adding u of OMIM annotation ##
    x = tmp3

    if(nrow(x) >=1){x$MOI <- 'inherited'}
    return(x)

}

