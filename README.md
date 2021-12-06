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
	
	* Processing annotation sources: All of the above annotation sources can be processed using following command:
		$ python processDB.py -a <user-config-xml-file>
				      -d <database-name>
				      -r <ref-genome-version>
				      -l <lift-over-flag>
				      -m <manifest-file-internal-cohort-NGC>
				      -p <project-name>
	* Processing gnomad SVs
		Example command to process gnomad SVs (grch37 version):
		- Download the gnomad SVs (bed format) from here: https://storage.googleapis.com/gcp-public-data--gnomad/papers/2019-sv/gnomad_v2.1_sv.sites.bed.gz 
		- Place in <resource-path>gnomad/grch37/ directory. 
		- Run the command:
		$ python processDB.py -a CONFIG/Analysis_user.grch37.xml -d gnomad 
				      -r grch37 -l False -m NA -p NA
		Choose liftover flag to True for lifting over to GRCh38.
		$ python processDB.py -a CONFIG/Analysis_user.grch37.xml -d gnomad 
				      -r grch37 -l False -m NA -p NA

	* Processing NGC samples (Internal Cohort)
		- The sample vcfs have not been provided with this repo. But can be replicated for any other  internal cohort.
		- The samples should be spcified as in example manifest file: 
			<resource-path>/examples/manifest.txt
		- Example command:
		$ python processDB.py -a CONFIG/Analysis_user.grch37.xml -d ngc
				      -r grch37 -l false -m <resource-path>/examples/manifest.txt
				      -p 20211203
	
	* Processing Ensembl
		- Download Ensembl gene information file from: http://ftp.ensembl.org/pub/grch37/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz
		- Example command:
		$ python processDB.py -a CONFIG/Analysis.v1.xml -r grch37 -d ensembl 
				      -l false -m NA -p NA
	* Processing RefSeq
		- Standardizing the raw Refseq gene file by mapping the GenBank-Accn to Chromosome number 
		- Download RefSeq gene information files 
		(1) Genomic : https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.gtf.gz
		(2) Assembly report: https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_assembly_report.txt
		- Put these files in following directory:
			<resource-path>/refseq/grch37/
		- Run this command for standardization:
		$ python convertRefSeqAccession.py 
			-i <resource-path>/refseq/grch37/GRCh37_latest_genomic.gtf.gz 
			-a <resource-path>/refseq/grch37/GRCh37_latest_assembly_report.txt
		- Example command to create the database:
		$ python processDB.py -a CONFIG/Analysis_user.grch37.xml -d refseq
				      -l false -m NA -p NA

	* Processing Decipher
		- Download the file: https://www.deciphergenomics.org/files/downloads/population_cnv_grch37.txt.gz
		- Put it in following directory: <resource-dir>/decipher/grch37/
		- Run the command:
		$ python processDB.py -a CONFIG/Analysis_user.grch37.xml -d decipher
				      -l false -m NA -p NA
	 
	* Processing dbVar
		- Download these files:
		 	(1) Clingen (nstd45):https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/vcf/nstd45.GRCh37.variant_call.vcf.gz
			(2) User curated (nstd51): https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/vcf/nstd51.GRCh37.variant_call.vcf.gz
			(3) Clinvar (nstd102): https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/vcf/nstd102.GRCh37.variant_call.vcf.gz
		- Put these files in: <resource-path>/dbvar/grch37/
		- Example run of the command
		$ python processDB.py -a CONFIG/Analysis_user.grch37.xml -d dbvar
				      -l false -m NA -p NA
	
	* Processing User specific
		- Promoter regions (Ensembl)
		- Download: http://ftp.ensembl.org/pub/grch37/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz
		- Put it here: <resource-path>/user/promoter/grch37/
		- Example command:
		$ python processDB.py -a CONFIG/Analysis_user.grch37.xml -d user_promoter
				      -l false -m NA -p NA
	
		
		- Blacklist regions
		- Download: https://github.com/Boyle-Lab/Blacklist/blob/master/lists/Blacklist_v1/hg38-blacklist.bed.gz
		- Put it here: <resource-path>/user/blacklist/grch37/
		- Example command:
		$ python processDB.py -a CONFIG/Analysis_user.grch37.xml -d user_blacklist
				      -l false -m NA -p NA
	 
	
	
	 

					
