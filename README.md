# SV pipeline

Structural variant (SV) calling pipeline explicity developed to process individual or trio vcf files containing SV classes (DEL/DUP/INS/INV/CNVs) aligned to either grch37/grch38 version of reference genome.

The pipeline requires annotation sources, customized datasets, available tools (vcfanno,bcftools) and input set of SVs in vcf format. It generates analysis scripts that can be incorporated into any high performance clusters (HPC). This results in list of filtered variants per family that can be used for interpretation, reporting and downstream analysis/

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
	$ pip install numpy
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
	$ conda install -c conda-forge r-future.apply
	$ conda install -c conda-forge r-r.utils 
	$ conda install -c bioconda ucsc-liftover
	$ conda install -c bioconda vcfanno

#
# Annotation sources
#

##### Building up of Annotation sources ############

	* With this release we are providing all the respective annotation sources downloaded and processed for both GRCh37 and GRCh38 version.
		- The example data can be found via this link: https://drive.google.com/drive/folders/1pCJj2KbQScpy1XSGTTDx70ewgwc5tirn?usp=sharing 
		- Please extract the demo.tar.gz file under your chosen directory. It contains:
			(a) resource directory [demo/resources/]
			(b) examples directory [demo/examples/]
		- Please define your own <resource-path> and append to the demo folder path. Edit the CONFIG/Analysis.xml file:
		<general>
			<resourceDir> your own path </resourceDir>
		</general>
		
	* Annotation of SVs in this pipeline is done by vcfanno tool [Brent et.al 2016]. 
	  This requires annotation sources to be processed in specific tab separated format. 
	  The format is defined as:
		#Chrom	#Start-Pos	#End Position	#SV-IDs,Other information
		1	100	200	SV_ID_1
		1	105	500	SV_ID_2 

	* The pipeline requires following annotation sources categorised as:
		1. Population databases:
			1.1. gnomAD SVs
			1.2. Internal cohort (NGC samples)	
		2. Gene information:
			2.1. Ensembl: Transcript, Exons
			2.2. Refseq: Gene, Exons 
		3. Pathogenecity
			3.1. Decipher
			3.2. dbVar: Clingen(nstd45),User(nstd51),Clinvar(nstd102)
		4. Others:
			4.1. Promoter regions (Ensembl)
			4.2. Blacklist regions

	
	* Processing annotation sources: All of the above annotation sources can be processed using following command:
		$ python processDB.py -a <user-config-xml-file>
				      -d <database-name>
				      -r <ref-genome-version>
				      -l <lift-over-flag>
				      -m <manifest-file-internal-cohort-NGC>
				      -p <project-name (optional)>

	All the annotation sources including raw and processed (grch37/grch38) are currently provided with this release. See in demo/resource/ folders. However for demonstration purpose following commands are used:

	* Processing gnomAD SVs
		Example command to process gnomad SVs (grch37 version):
		- Download the gnomad SVs (bed format) from here: https://storage.googleapis.com/gcp-public-data--gnomad/papers/2019-sv/gnomad_v2.1_sv.sites.bed.gz 
		- Place in <resource-path>gnomad/grch37/ directory. 
		- Run the command:
		
		$ python processDB.py -a CONFIG/Analysis.xml -d gnomad 
				      -r grch37 -l True 
		Choose liftover flag to True for lifting over to GRCh38 version
		$ python processDB.py -a CONFIG/Analysis.xml -d gnomad 
				      -r grch37 -l True

	* Processing Internal cohort (Example cases )
		- The sample vcfs have not been provided with this repo. But can be replicated for any other  internal cohort.
		- The samples should be spcified as in example manifest file: 
			demo/examples/manifest.txt
		
		- Example command:
		$ python processDB.py -a CONFIG/Analysis.xml -d example
				      -r grch38 -l false -m <resource-path>/demo/examples/manifest.txt
				      
	
	* Processing Ensembl
		- Download file: (grch37) http://ftp.ensembl.org/pub/grch37/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz
		
		- Download file: (grch38) http://ftp.ensembl.org/pub/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh38.98.gtf.gz
		- Put these files in respective grch37/grch38 directory: <resource-path>/demo/resources/ensembl/grch38/
		- Example run for grch38 version:
		$ python processDB.py -a CONFIG/Analysis.xml -r grch38 -d ensembl -l false 
		$ Offset
			- For transcript regions: ±1000 bp either ends.
			- For exonic regions: ±50 bp either ends.

	* Processing RefSeq
		- Standardizing the raw Refseq gene file by mapping the GenBank-Accn to Chromosome number 
		- Download RefSeq gene information files (grch37) 
		(1) Genomic : https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.gtf.gz
		(2) Assembly report: https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_assembly_report.txt

		- Download RefSeq gene information files (grch38) 
		(1) Genomic : https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gtf.gz
		(2) Assembly report: https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_assembly_report.txt

		- Put these files in respective grch37 or grch38 directory:
			<resource-path>/demo/resources/refseq/grch38/
		- Run this command for standardization:
		
		$ python convertRefSeqAccession.py 
			-i <resource-path>/refseq/grch38/GRCh38_latest_genomic.gtf.gz 
			-a <resource-path>/refseq/grch38/GRCh38_latest_assembly_report.txt
			-o <resource-path>/refseq/grch38/GRCh38_latest_genomic.map.gtf.gz
		
		- Example command to create the database:
		$ python processDB.py -a CONFIG/Analysis.xml -d refseq -r grch38 -l false

		$ Offset
			- For transcript regions: ±1000 bp either ends.
			- For exonic regions: ±50 bp either ends.

	* Processing Decipher
		- Compilation source: https://www.deciphergenomics.org/disorders/syndromes/list
		- Provided default with this repo directory: <resource-dir>/demo/resources/decipher/grch37/decipher_cnv_syndromes.bed
		and <resource-path>/demo/resources/decipher/grch38/decipher_cnv_syndromes_liftover.b38.bed

		- Example run for grch38 version:
		
		$ python processDB.py -a CONFIG/Analysis.xml -d decipher -r grch38 -l false
	
	* Processing dbVar
		- Download these files (grch37):
		 	(1) Clingen (nstd45):https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/vcf/nstd45.GRCh38.variant_call.vcf.gz
			(2) User curated (nstd51): https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/vcf/nstd51.GRCh38.variant_call.vcf.gz
			(3) Clinvar (nstd102): https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/vcf/nstd102.GRCh38.variant_call.vcf.gz
		- Put these files in: <resource-path>/demo/resources/dbvar/grch37/

		- Download these files (grch38):
		 	(1) Clingen (nstd45):https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/vcf/nstd45.GRCh38.variant_call.vcf.gz
			(2) User curated (nstd51): https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/vcf/nstd51.GRCh38.variant_call.vcf.gz
			(3) Clinvar (nstd102): https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/vcf/nstd102.GRCh38.variant_call.vcf.gz
		- Put these files in: <resource-path>/demo/resources/dbvar/grch38/
		- Example run of the command for grch38 version:
		
		$ python processDB.py -a CONFIG/Analysis.xml -d dbvar -r grch38 -l false
	
	* Processing  specific databases:
		- Promoter regions (Ensembl)
		- Download: (grch37) http://ftp.ensembl.org/pub/grch37/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz
		- Put it here: <resource-path>/demo/resources/promoter/grch37/

		- Download: (grch38) http://ftp.ensembl.org/pub/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh38.98.gtf.gz
		- Put it here: <resource-path>/demo/resources/promoter/grch38/
		
		- Example command for running grch38 version:
		
		$ python processDB.py -a CONFIG/Analysis.xml -d promoter -r grch38 -l false
	
		
		- Blacklist regions
		- Download: (grch37) http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg19-human/wgEncodeHg19ConsensusSignalArtifactRegions.bed.gz 
		- Put it here: <resource-path>/demo/resources/blacklist/grch37/

		- Download: (grch38): http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg38-human/hg38.blacklist.bed.gz
		- Put it here: <resource-path>/demo/resources/blacklist/grch37/
		
		- Example command for processing grch38 version:
		
		$ python processDB.py -a CONFIG/Analysis.xml -d blacklist -r grch38 -l false
	
		
		
#
# Running SV annotation pipeline
#

##### SV pipeline using Vcfanno ############
		
After creating the annotation source the input vcf files can be annotated using following command. It generates shell scripts. Generally for large cohort analysis it can be integrated into any HPC services such as SLURM/LSF/PBS. For a given trio family, the pipeline annoates the input SVs of proband. The main steps are:

	(1) Normalizing vcfs for multi-allelic sites.
	(2) Annotation using Vcfanno by incorporating above databases.
	(3) Extract relevant fields of annotation splitted per annotation source  (using bcftools).
	(4) Compute overlap of query SVs with annotation source:
		- Reciprocal overlap (default 70%) with gnomAD, Internal cohort(such as NGC)
		- Single base overlap with Ensembl, Refseq, Decipher, dbVar, promoter, blacklist regions.
	(5) Merge all the extracted fields
	(6) Incorporate the frequency and inheritance filtering on the anntated SVs.

N.B: Current implementation steps of the pipeline (step 1-6) has been optimized for NGC project WGS samples. However for demonstration purpose we provide example scripts that were used to annotate The International Genome Sample Resource (IGSR) trio SVs called by LUMPY. See demo/examples/samples/ directory

	- Example command to generate the scripts:
	Input requirement: The manifest file which holds pedigree information and their correpsonding SV vcf file.

	$ python processSV.py 	-m <manifest-file>
				-a <config-XML-file>
				-w <work-dir>
				-e <analysis-type:vcfanno_demo/vcfanno_ngc>
				-f <list-of-proband-id-to-be-analyzed>
				-r <ref-genome-version>
				-l <HPC cluster launch flag>
				-p <project-date>

	Example run:
	$ python processSV.py -m <resource-path>/demo/examples/manifest.txt 
			-a CONFIG/Analysis.xml 
			-w <resource-path>/demo/examples/ 
			-e vcfanno_demo 
			-f <resource-path>/demo/examples/family_id.txt 
			-r grch38 -l F -p 20211210
	
	N.B: For replicating the same for NGC-WGS-cohort, please use the analysis type: vcfanno_ngc
		

	- Output: This command will generate the shell scripts in <resource-path>demo/examples/20211210/tmp_binaries/NGC001_01.vcfanno_demo.sh
	- One can launch the script their own preferred HPC scheduler such as SLURM/LSF/PBS by adding relevant memory and CPU configurations.
	- With the provided demo/examples the SV trio data (IGSR) only the final merged annotated SVs (called by LUMPY) of the proband can be obtained. The filtering pipeline has been tested and customized for NGC-WGS-cohort samples which were eventually called vby Manta & Canvas.
	
#
# Contact details
# 

Please contact us for any quesries related to pipeline:
Ajay : aak@ebi.ac.uk
Courtney : cf458@cam.ac.uk



