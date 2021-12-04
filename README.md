# SV pipeline

Structural variant (SV) calling pipeline explicity developed to process individual or trio vcf files containing SVs aligned to either (grch37/grch38).

The pipeline requires annotation sources, customized datasets, available tools (vcfanno) and input set of SVs in vcf format. It generates analysis scripts that can be incorporated into any high performance clusters (HPC). This results in list of filtered variants per family that can be used for interpretation, reporting and downstream analysis/

For demonstration purpose below example is presented for GRCh37. However, the same can be replicated for GRCh38.

# Installation
	git clone https://github.com/ajaarma/sv.git
	
#
# Required installation packages
#

##### Install anaconda  #####
	Follow this link for installation: https://docs.anaconda.com/anaconda/install/linux/

##### Conda environment commands ##########

	$ conda create --name sv
	$ source activate sv
	$ conda install python=2.7.16
	$ conda install -c bioconda vcfanno
	$ pip install xmltodict
	$ pip install dicttoxml
	$ conda install numpy
	$ conda install -c bioconda ucsc-liftover
	$ pip install vcf_parser
	$ conda install samtools=1.3
	$ conda install vcftools=0.1.14
	$ conda install bcftools=1.9
	$ conda install gcc #(OSX)
	$ conda install gcc_linux-64 #(Linux)
	$ conda install parallel
	$ conda install -c r r-optparse
	$ conda install -c r r-dplyr
	$ conda install -c r r-plyr
	$ conda install -c r r-data.table
	$ conda install -c bioconda ucsc-liftover
	$ conda install -c bioconda vcfanno

#
# Annotation sources
#

##### Building up of Annotation sources ############
	* Annotation of SVs in this pipeline is done by vcfanno tool[Reference]. 
	This requires annotation sources to be processed in specific tab separated format. 
	The format is defined as:
		#Chrom	#Start-Pos	#End Position	#SV-IDs,Other information
		1	100	200	SV_ID_1
		1	105	500	SV_ID_2 

	* The pipeline requires following annotation sources categorised as:
		1. Population databases:
			1.1. gnomad SVs
			1.2. NGC samples (internal cohort)	
		2. Gene information:
			2.1. Ensembl: Promoter, Transcript, Exons
			2.2. Refseq: Gene, Exons 
		3. Pathogenecity
			3.1. Decipher
			3.2. dbVar: Clingen(nstd45),User(nstd51),Clinvar(nstd102)
		4. Others:
			4.1. Blacklist regions
	* 

					
