#!/usr/bin/python

from xml.etree import ElementTree
from xml.dom import minidom
from xml.etree.ElementTree import Element, SubElement, Comment, tostring
import xmltodict,os
import argparse

######################################################################################
#
# Description: Generic class for creating XML based configuration files; 
#
# Dependencies: Please install the xmltodict package for python. Do this:
#       (1) Create the directory:  mkdir -p ~/.local/lib/python2.7/site-packages/
#       (2) Download the package: wget https://github.com/martinblech/xmltodict/archive/master.zip
#       (3) python setup.py install --prefix=$HOME/.local
#
#
######################################################################################

class CONFIG:

    def __init__(self, elements=[]):
        self.__elements={}
        for e in elements:
            self.__elemets[e]=1

    def str2bool(self,v):
        ''' Utility function to convert string/unicode objects ('yes','true','t','1')
            of boolean type (True/False) '''
        print v
        v = v[0]        
        return v.lower() in ('yes','true','t','1')

    def display(self):
        print "Inside CONFIG class. Creating a configuration file for the NGC-SNV pipeline"

    def processXMLArguments(self):
        ''' Subroutine to process user configuration file '''
        
        cmd_dict = {}
        parser = argparse.ArgumentParser(
                     description = 'python createAnalysisXML.py -u <user-config> \n'
                                           '-b <base-xml-file> -o <out-xml-file> \n'
                    )

        parser.add_argument('-u','--userConfig',help='User specified config file',
                                    action='store',dest='userConfig',required=True)
        parser.add_argument('-b','--baseXml',help='Base XML config file',
                                    action='store',dest='baseXML',required=True)
        parser.add_argument('-o','--outXml',help='Output XML file',
                                    action='store',dest='outXML',required=True)

        self.results = parser.parse_args()
        self.user_config = os.path.abspath(self.results.userConfig)
        self.base_xml = os.path.abspath(self.results.baseXML)
        self.out_xml = os.path.abspath(self.results.outXML)

        cmd_dict = {'userConfig':self.user_config,'baseXml':self.base_xml,'outXml':self.out_xml}

        return cmd_dict

    def parseSVCommandArgs(self):

        cmdDict = {}

        parser = argparse.ArgumentParser(
                 description = 'python processSV.py -m <manifest-file> \n' 
                                '-a <analysis.xml> -p <project-date> '
                                '-e <analysis-type> -w <work-dir> '
                                '-f <family-id-file> -r <genome-build> '
                                '-l <launch-flag> '
                )

        parser.add_argument('-m','--mainfest',help='Manifest file',
                            action='store',dest='manifestFile',required=True)
        parser.add_argument('-a','--analXML',help='Analysis XML file: Analysis.xml',
                            action='store',dest='xmlFile',required=True)
        parser.add_argument('-p','--projDate',help='Project date: 20100906',
                            action='store',dest='projDate',required=True)
        parser.add_argument('-e','--analType',help='Analysis type: vcfanno/annotSV',
                            action='store',dest='expType',required=True)
        parser.add_argument('-w','--workDir',help='Working directory', 
                            action='store',dest='workDir',required=True)
        parser.add_argument('-f','--famFile',help='Family IDs for analysis', 
                            action='store',dest='famFile',required=True)
        parser.add_argument('-r','--ref',help='Genome build version: GRCh37/GRCh38', 
                            action='store',dest='ref',required=True,default='v38')
        parser.add_argument('-l','--launch',help='Launch the SV pipeline:True/False', 
                            action='store',dest='launchFlag',required=True)
        parser.add_argument("-v","--version",help="show program's version and exit\n\n\n", 
                            action="version",version='%(prog)s 1.0')

        self.results = parser.parse_args()
        self.manifest_file = os.path.abspath(self.results.manifestFile)
        self.xml_file = os.path.abspath(self.results.xmlFile)
        self.proj_date = self.results.projDate
        self.anal_type = self.results.expType
        self.work_dir = os.path.abspath(self.results.workDir)
        self.fam_file = os.path.abspath(self.results.famFile)
        self.ref_built = self.results.ref
        self.launch_flag = self.str2bool(self.results.launchFlag)
        
        cmd_dict = {'manifest':self.manifest_file,'xml':self.xml_file,
                    'proj':self.proj_date,'analType':self.anal_type,
                    'workDir':self.work_dir,'famFile':self.fam_file,
                    'ref':self.ref_built,'launch':self.launch_flag
                   }
        print "\n\n"
        print 'Entered Manifest file is: ',''.join(self.manifest_file)
        print 'Entered Analysis file is: ',''.join(self.xml_file)
        print 'Entered Project date  is: ',''.join(self.proj_date)
        print 'Entered Analysis step is: ',''.join(self.anal_type)
        print 'Entered Work directory is: ',''.join(self.work_dir)
        print 'Entered family IDs file is: ',''.join(self.fam_file)
        print 'Entered Analysis genome build is: ',''.join(self.ref_built)
        print 'Entered SV pipeline launch flag: ', self.launch_flag

        print "\n\n"

        #return (self.manifest_file, self.xml_file, self.proj_date, self.anal_type, 
        #        self.work_dir, self.fam_file, self.ref_file)
        return cmd_dict
        

    def parseVCFCommandArgs(self):
        
        parser = argparse.ArgumentParser()
        parser.add_argument("-i","--input",help="input annotated vcf file",
                                action="store",dest="InputFile",required=True)
        #parser.add_argument("-o","--output",help="output annotated vcf file",action="store",dest="OutFile",required=True)
        parser.add_argument("-f","--familyID",help="familyID",action="store",
                                                    dest="famID",required=True)
        parser.add_argument("-a","--annotation",help="annotation databases used:\
                                gnomad/ngc38/ngclov",action="store",dest="dbType",
                                                                    required=True)
        parser.add_argument("-t","--task",help="execution task:overlap/freq",
                                        action="store",dest="task",required=True)
        parser.add_argument("-x","--config",help="Configuration XML file",
                                        action="store",dest="xml",required=True)
        parser.add_argument("-d","--debug",help="debugging the code. Enter the \
                              variant line number to start debugging. Default is line number =1",
                              action="append",dest="collection",default=[],required=False)
        parser.add_argument("-m","--manifest",help="Manifest File",action="store",
                                                    dest="manifest",required=True)
        parser.add_argument("-r","--ovFrac",help="Overlap Fraction:",action="store",
                                            dest="ovFrac",required=True,default=0.70)
        parser.add_argument("-v","--version",help="show program's version and exit\n",
                                                action="version",version='%(prog)s 1.0')
    
        self.results = parser.parse_args()
        self.inp_file = os.path.abspath(self.results.InputFile)
        self.fam_id = self.results.famID
        self.xml_file = os.path.abspath(self.results.xml)
        self.task     = self.results.task
        self.db_type  = self.results.dbType.lower()
        self.manifest_file = os.path.abspath(self.results.manifest)
        self.line_index = self.results.collection
        self.ovFrac = float(self.results.ovFrac)

        if self.line_index:
            self.debug_flag = True
            self.line_index = self.line_index[0]
        else:
            self.debug_flag=False
            self.line_index = 1

        return self.inp_file,self.fam_id, self.line_index,self.manifest_file, self.xml_file,self.db_type,self.task,self.debug_flag, self.ovFrac

    def parseTomlCommandArgs(self):
        ''' Method to parse the arguments for creating Toml file as list of annotation sources used '''

        parser = argparse.ArgumentParser()
        parser.add_argument('-a','--analXml',help='Analysis xml file',action='store',dest='xmlFile',required=True)
        parser.add_argument('-v','--version',help='show programs\'s version and exit\n',action='version',version='%(prog)s 1.0')

        self.results = parser.parse_args()
        self.xml_file = results.xmlFile

        return self.xml_file

    def parseIntOverlapCommandArgs(self):

        ''' Method to parse the argument for computing Internal Overlap withing family '''

        parser = argparse.ArgumentParser()
        parser.add_argument("-m","--manifest",help="Input manifest file",action="store",dest="manifestFile",required=True)
        parser.add_argument("-o","--output",help="output overlap file",action="store",dest="outFile",required=True)
        parser.add_argument("-f","--familyID",help="familyID",action="store",dest="famID",required=True)
        parser.add_argument("-w","--window",help="window length of base pairs (bp)",action="store",dest="winLen",required=True)
        parser.add_argument("-d","--debug",help="debugging the code. Enter the variant line number to start debugging. Default is line number =1",action="append",dest="collection",default=[],required=False)
        parser.add_argument("-v","--version",help="show program's version and exit\n",action="version",version='%(prog)s 1.0')
    
        results = parser.parse_args()
        m_file = results.manifestFile
        fam_id = results.famID
        out_file = results.outFile
        win_len = results.winLen
        line_index = results.collection

        if line_index:
            debug_flag = True
            line_index = line_index[0]
        else:
            debug_flag=False
            line_index = 1

        return [m_file, out_file, fam_id, win_len, debug_flag, line_index]

    def parseDBCommandArgs(self):
        ''' Method to parse the argument for processing Annotation DBs '''

        parser = argparse.ArgumentParser(
                 description = 'Process SV annotation sources'
                 #usage = 'python processSDB.py -m <manifest-file> \n' 
                 #              '-a <analysis.xml> -p <project-date> '
                 #              '-r <genome-build> -l <lift-over-flag:0/1;\
                 #               0=> False and 1=>True >' 
                 )

        parser.add_argument('-m','--manifest',help='manifest file',
                            action='append',dest='manifestFile',required=False)
        parser.add_argument('-a','--configxml',help='Analysis XML file',action='store',
                            dest='xmlFile',required=True)
        parser.add_argument('-d','--db',help='Database name: dbVar,NGC,GNOMAD,Ensembl\
                            DECIPHER',action='store',dest='dbType',required=True)
        parser.add_argument('-r','--reference',help='Reference genome: v37/v38',
                            action='store',dest='refGenome',required=True)
        parser.add_argument('-l','--liftOverFlag',help='Lift over flag:'+ 
                            'True/False ',action='append',dest='lovFlag',
                            required=False)
        parser.add_argument('-p','--projectDate',help='Project date',
                            action='append',dest='project',required=False)

        self.results = parser.parse_args()
        self.m_file = self.results.manifestFile
        self.xml_file = self.results.xmlFile
        self.db_name = self.results.dbType
        self.ref_file = self.results.refGenome
        self.lov_flag = self.str2bool(self.results.lovFlag)
        self.proj_date = self.results.project

        cmd_dict = {'manifest':self.m_file,'analysis':self.xml_file,
                    'db':self.db_name,'ref':self.ref_file,'lov':self.lov_flag,
                    'project':self.proj_date}

        print "\n\n"
        print 'Entered Manifest file is: ',self.m_file
        print 'Entered Analysis file is: ',''.join(self.xml_file)
        print 'Entered Project date  is: ',self.proj_date
        print 'Entered database is: ',''.join(self.db_name)
        print 'Entered refernce genome: ',''.join(self.ref_file)
        print 'Entered lift over flag is: ',self.lov_flag
        print "\n\n" 

        return cmd_dict
    
    def parseMergeAnnotCommandArgs(self):
        ''' Method to process command line arguments for merge annotation 
            steps '''
        parser = argparse.ArgumentParser(
                 description = 'python mergeAnnotVCF.py -m <manifest-file> \n' 
                               '-i <input-rest-file> -a <annot-dir>'
                               '-f <family-id>  -r <genome-build>'
                               '-o <overlap-frac> -s <analysis-type>' 
                 )
        parser.add_argument('-m','--mainifest',help='Input manifest file',
                            action='store',dest='manifestFile',required=True)
        parser.add_argument('-i','--rest',help='Input SV coordinates file',
                            action='store',dest='restFile',required=True)
        parser.add_argument('-a','--annoDir',help='Annotation Directory',
                            action='store',dest='annoDir',required=True)
        parser.add_argument('-f','--famID',help='Family ID',action='store',
                            dest='famID',required=True)
        parser.add_argument('-r','--ref',help='reference built',action='store',
                            dest='refFile',required=True)
        parser.add_argument('-o','--frac',help='Overlap fraction',
                             action='store',dest='overlapFrac',required=True)
        parser.add_argument('-s','--analType',help='Analysis Type: ngc/demo',
                             action='store',dest='analType',required=True)


        self.results = parser.parse_args()
        self.m_file = os.path.abspath(self.results.manifestFile)
        self.rest_file = os.path.abspath(self.results.restFile)
        self.anno_dir = os.path.abspath(self.results.annoDir)
        self.fam_id = self.results.famID
        self.ref_file = self.results.refFile
        self.overlap_frac = self.results.overlapFrac
        self.anal_type = self.results.analType

        cmd_dict = {'manifest':self.m_file,'rest':self.rest_file,'annoDir':self.anno_dir,
                    'fam':self.fam_id,'ref':self.ref_file,'overlap':self.overlap_frac,
                    'analType':self.anal_type
                   }

        print "\n\n"
        print 'Entered Manifest file is: ',self.m_file
        print 'Entered Rest annotation file is: ',self.rest_file
        print 'Entered Annotation directory is: ',self.anno_dir
        print 'Entered Family ID is: ',self.fam_id
        print 'Entered Reference genome built is: ',self.ref_file
        print 'Entered Overlap fraction is: ',self.overlap_frac
        print 'Entered Analysis type is: ',self.anal_type

        print "\n\n"

        return cmd_dict

    def parseFormatGTCommandArgs(self):
        ''' Method to process command line arguments for formatting GT 
            steps '''
        parser = argparse.ArgumentParser(
                 description = 'python formatGT.py -m <manifest-file> \n' 
                               '-i <input-rest-file>'
                               '-f <family-id>'
                               '-p <par-region-file>'
                               '-o <out-formatted-file> ' 
                 )
        parser.add_argument('-m','--mainifest',help='Input manifest file',
                            action='store',dest='manifestFile',required=True)
        parser.add_argument('-i','--inpFile',help='Input merge annotated file',
                            action='store',dest='inpFile',required=True)
        parser.add_argument('-o','--outFile',help='Output Formatted file',
                            action='store',dest='outFile',required=True)
        parser.add_argument('-f','--famID',help='Family ID',action='store',
                            dest='famID',required=True)
        parser.add_argument('-r','--ref',help='reference built',action='store',
                            dest='refFile',required=True)
        parser.add_argument('-p','--par',help='PAR region file',action='store',
                            dest='parFile',required=True)

        self.results = parser.parse_args()
        self.m_file = os.path.abspath(self.results.manifestFile)
        self.inp_file = os.path.abspath(self.results.inpFile)
        self.out_file = os.path.abspath(self.results.outFile)
        self.fam_id = self.results.famID
        self.par_file = self.results.parFile
        self.ref_file = self.results.refFile

        cmd_dict = {'manifest':self.m_file,'inp':self.inp_file,'out':self.out_file,
                    'fam':self.fam_id,'par':self.par_file,'ref':self.ref_file}

        print "\n\n"
        print 'Entered Manifest file is: ',self.m_file
        print 'Entered input file is: ',self.inp_file
        print 'Entered out file is: ',self.out_file
        print 'Entered Family ID is: ',self.fam_id
        print 'Entered PAR region file is: ',self.par_file
        print 'Entered Reference built is: ',self.ref_file
        print "\n\n"

        return cmd_dict

    def parseRefseqArgs(self):
        ''' Method to process command line arguments for formatting Refseq file 
        '''
        parser = argparse.ArgumentParser(
                 description = 'python convertRefseqAccession.py -i <refseq-file> \n' 
                               '-r <assembly-report-file>'
                               '-o <out-file>'
                 )
        parser.add_argument('-i','--inpFile',help='Input Refseq file',
                            action='store',dest='inpFile',required=True)
        parser.add_argument('-o','--outFile',help='Output Formatted file',
                            action='store',dest='outFile',required=True)
        parser.add_argument('-r','--report',help='Assembly Report File',action='store',
                            dest='report',required=True)

        self.results = parser.parse_args()
        self.inp_file = os.path.abspath(self.results.inpFile)
        self.out_file = os.path.abspath(self.results.outFile)
        self.report_file = self.results.report

        cmd_dict = {'inp':self.inp_file,'out':self.out_file,
                    'report':self.report_file}

        print "\n\n"
        print 'Entered input file is: ',self.inp_file
        print 'Entered out file is: ',self.out_file
        print 'Entered report file is: ',self.report_file
        print "\n\n"

        return cmd_dict


    def processInit(self,work_dir,cmd_args, proj_date):

        """ Creating temporary directories inside the work directory """
        
        print "Creating temporary directories inside the work directory \n"

        self.date_str = proj_date
        self.tmp_dir = os.path.abspath(self.work_dir+'/'+self.date_str)
        self.tmp_bin = os.path.abspath(self.tmp_dir+'/tmp_binaries')
        self.tmp_data = os.path.abspath(self.tmp_dir+'/tmp_data')
        self.tmp_status = os.path.abspath(self.tmp_dir+'/tmp_status')
        self.tmp_log = os.path.abspath(self.tmp_dir+'/tmp_log')
        self.sb_log = os.path.abspath(self.tmp_dir+'/sb_log')
    
        self.tmp_dict = {'tmpDir':self.tmp_dir,'tmpBin':self.tmp_bin,'tmpData':self.tmp_data,
                    'tmpStat':self.tmp_status,'tmpLog':self.tmp_log,'sbLog':self.sb_log
                    }

        if os.path.isdir(self.work_dir+'/'+self.date_str):

            print ' -- The entered project space already exists. Do you want to delete it?'
            self.text = raw_input('Enter y/n: ')
            print self.text

            if self.text.lower() == 'y':
                print ' -- Deleting existing project. Creating fresh workspace with name: ',self.date_str
                os.system('rm -rf '+self.work_dir+'/'+self.date_str)

                for keys,values in self.tmp_dict.items():
                    os.system('mkdir -p '+values)

            elif self.text.lower() == 'n':
                pass
        else:
            for keys,values in self.tmp_dict.items():
                os.system('mkdir -p '+values)

        wh = open(self.tmp_dir+"/Command_Entered.txt","a+")
        print >>wh,"python "+" ".join(cmd_args)
        wh.close()

        return self.tmp_dict


    def prettify(self,elem):

        self.rough_string = ElementTree.tostring(elem,'utf-8')
        self.reparsed = minidom.parseString(rough_string)
        
        return reparsed.toprettyxml(indent="  ")

    def getConfigTop(self, version_type):

        self.top = Element('config')
        self.comment = Comment("Configuration file for the NGC-SNV pipeline")
        self.top.append(self.comment)
        self.version = SubElement(self.top,"version")
        self.version.text = version_type

        return self.top

    def getGeneral(self):
        pass
  
    def getConfigDict(self, xml_file,module_name=[]):

        print '\nReading and processing configuration-analysis-xml file: ',xml_file
        with open(xml_file) as fd:
            self.doc = xmltodict.parse(fd.read())

        self.doc = self.doc["config"]
        print "\n -- The item classes present in Configuration files are: "+",".join(self.doc.keys())

        if module_name:
            self.config_dict = self.doc[module_name]
        else:
            self.config_dict = self.doc

        return self.config_dict

    def getUserConfigDict(self,userConfigFile):
        ''' Subroutine to read process the User defined configuration file and
            return as dictionary object 
        '''

        fh = open(userConfigFile)

        for lines in fh:
            lines = lines.strip()
            strs = re.split('\n',lines)
            
            if re.search('\=',strs[0]) and strs[0] !='':
                pass
            elif strs[0]:
                print 'The input user configuration file is not in correct format'+\
                        'Please read the online documentation '
                print 'The faulty line is: ',lines
                sys.exit()
        
        fh.close()
        
        userConfigDict = OrderedDict()

        fh = open(userConfigFile)
        for lines in fh:
            lines = lines.strip()
            strs = re.split('\=',lines);strs=[x.strip() for x in strs]
            if strs[0] and strs[1]:
                userConfigDict[strs[0]] = strs[1]
        
        fh.close()

        return userConfigDict

    def updateBaseXML(self,baseConfigDict,userConfigDict):
        ''' Subroutine to update the Base XML file with user specified 
            configuration file 
        '''

        for user_key in userConfigDict:
            if re.search('genomeBuild',user_key):
                try:
                    baseConfigDict['general'][user_key] = userConfigDict[user_key]
                except KeyError:
                    print 'The corresponding xml-field: ',user_key,' doesnot exist'
                    print 'Please enter correct xml-field tag: genomeBuild'
                
            if re.search('resourceDir',user_key):
                try:
                    baseConfigDict['general'][user_key] = userConfigDict[user_key]
                except KeyError:
                    print 'The corresponding xml-field: ',user_key,' doesnot exist'
                    print 'Please enter correct xml-field tag: resourceDir'
            if re.search('refFasta',user_key):
                try:
                    baseConfigDict['general'][user_key] = userConfigDict[user_key]
                except KeyError:
                    print 'The corresponding xml-field: ',user_key,' doesnot exist'
                    print 'Please enter correct xml-field tag: refFasta'
            
            if re.search('gnomad_g',user_key):
                try:
                    baseConfigDict['annotation']['cust-annot']['GNOMADg']['vcf'] = userConfigDict[user_key]
                except KeyError:
                    print 'The corresponding xml-field: ',user_key,' doesnot exist'
                    print 'Please enter correct xml-field tag: GNOMADg'
                 
            if re.search('gnomad_e',user_key):
                try:
                    baseConfigDict['annotation']['cust-annot']['GNOMADe']['vcf'] = userConfigDict[user_key]
                except KeyError:
                    print 'The corresponding xml-field: ',user_key,' doesnot exist'
                    print 'Please enter correct xml-field tag: GNOMADe'

            if re.search('dir_cache',user_key):
                try:
                    baseConfigDict['annotation']['defPar'][user_key] = userConfigDict[user_key]
                except KeyError:
                    print 'The corresponding xml-field: ',user_key,' doesnot exist'
                    print 'Please enter correct xml-field tag: <dir_cache>'

            if re.search('^exac$',user_key):
                try:
                    baseConfigDict['annotation']['cust-annot']['EXAC']['vcf'] = userConfigDict[user_key]
                except KeyError:
                    print 'The corresponding xml-field: ',user_key,' doesnot exist'
                    print 'Please enter correct xml-field tag: <EXAC>'
  
            if re.search('^exac_t$',user_key):
                try:
                    baseConfigDict['annotation']['cust-annot']['EXACt']['vcf'] = userConfigDict[user_key]
                except KeyError:
                    print 'The corresponding xml-field: ',user_key,' doesnot exist'
                    print 'Please enter correct xml-field tag: <exac_t>'
 
            if re.search('^clinvar$',user_key):
                try:
                    baseConfigDict['annotation']['cust-annot']['CLINVAR']['vcf'] = userConfigDict[user_key]
                except KeyError:
                    print 'The corresponding xml-field: ',user_key,' doesnot exist'
                    print 'Please enter correct xml-field tag: <CLINVAR>'
 
            if re.search('^hgmd$',user_key):
                try:
                    baseConfigDict['annotation']['cust-annot']['HGMD']['vcf'] = userConfigDict[user_key]
                except KeyError:
                    print 'The corresponding xml-field: ',user_key,' doesnot exist'
                    print 'Please enter correct xml-field tag: <HGMD>'
 
            if re.search('^cadd_snv$',user_key):
                try:
                    baseConfigDict['annotation']['plugIn']['CADD']['cadd_snv'] = userConfigDict[user_key]
                except KeyError:
                    print 'The corresponding xml-field: ',user_key,' doesnot exist'
                    print 'Please enter correct xml-field tag: <cadd_snv>'
 
            if re.search('^cadd_indel$',user_key):
                try:
                    baseConfigDict['annotation']['plugIn']['CADD']['cadd_indel'] = userConfigDict[user_key]
                except KeyError:
                    print 'The corresponding xml-field: ',user_key,' doesnot exist'
                    print 'Please enter correct xml-field tag: <cadd_indel>'
 
            if re.search('^exac_pli$',user_key):
                try:
                    baseConfigDict['annotation']['plugIn']['ExACpLI'] = userConfigDict[user_key]
                except KeyError:
                    print 'The corresponding xml-field: ',user_key,' doesnot exist'
                    print 'Please enter correct xml-field tag: <ExACpLI>'
 
            if re.search('^revel$',user_key):
                try:
                    baseConfigDict['annotation']['plugIn']['REVEL'] = userConfigDict[user_key]
                except KeyError:
                    print 'The corresponding xml-field: ',user_key,' doesnot exist'
                    print 'Please enter correct xml-field tag: <REVEL>'
 
            if re.search('^region_exons$',user_key):
                try:
                    baseConfigDict['regions']['exons'] = userConfigDict[user_key]
                except KeyError:
                    print 'The corresponding xml-field: ',user_key,' doesnot exist'
                    print 'Please enter correct xml-field tag: <exons>'
 
            if re.search('^region_par$',user_key):
                try:
                    baseConfigDict['regions']['PAR'] = userConfigDict[user_key]
                except KeyError:
                    print 'The corresponding xml-field: ',user_key,' doesnot exist'
                    print 'Please enter correct xml-field tag: <PAR>'

            if re.search('^ensembl$',user_key):
                try:
                    baseConfigDict['varPrior']['params']['l'] = userConfigDict[user_key]
                except KeyError:
                    print 'The corresponding xml-field: ',user_key,' doesnot exist'
                    print 'Please enter correct xml-field tag: <ensembl>'

            if re.search('^hpo$',user_key):
                try:
                    baseConfigDict['famFilter'][user_key] = userConfigDict[user_key]
                except KeyError:
                    print 'The corresponding xml-field: ',user_key,' doesnot exist'
                    print 'Please enter correct xml-field tag: <hpo>'

            if re.search('^phenotypes$',user_key):
                try:
                    baseConfigDict['famFilter'][user_key] = userConfigDict[user_key]
                except KeyError:
                    print 'The corresponding xml-field: ',user_key,' doesnot exist'
                    print 'Please enter correct xml-field tag: <phenotypes>'


        return baseConfigDict

    def convertTxt2XL(self,inpDir):
        ''' Utility function to convert list of Text files to Excel sheets '''

        #textfiles = [join(inpDir,f) for f in os.listdir(inpdir) 
        #                        if isfile(join(inpDir,f,'filtering')) and '.txt' in f ]

        wb = Workbook()

        for fam_id in os.listdir(inpDir):
            #fam_id = 'TwEx_EX1912835-TwEx_EX1912836-TwEx_EX1912837'
            try:
                txt_file = os.listdir(join(inpDir,fam_id,'filtering'))
                #print fam_id,'\t',fam_id[0]
                #print txt_file
                proband_id = re.split('\-',fam_id)[0]
                path_txt_file = join(inpDir,fam_id,'filtering',txt_file[0])
                #worksheet = workbook.add_sheet(str(proband_id))

                f = open(path_txt_file,'r+')
                data = f.readlines()
               
                proband_id = fam_id
                wb.create_sheet(proband_id)
                this_sheet = wb[proband_id]
                for row in data:
                    this_sheet.append(row.split('\t'))
                f.close()

            except IndexError:
                print fam_id

        wb.save(join(inpDir,'exeter_exomes_vars_famID.xlsx'))



