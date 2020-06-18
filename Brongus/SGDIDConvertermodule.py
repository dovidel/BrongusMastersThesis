import csv
from os import path

class SGD_gene(object):
    def __init__(self,SGDID,feature_type, status,feature_name,common_name,alt_names,chromosome,description):
        self.SGDID = SGDID.strip()
        self.feature_type= feature_type.strip()
        self.status = status.strip().upper()
        self.feature_name = feature_name.strip().upper()
        self.common_name = common_name
        self.alt_names = alt_names.split('|')
        self.chromosome = chromosome
        self.description = description
    
    def __repr__(self):
        #this function returns the desired gene name in a kind of tabular format that is easy for users to read. The \t makes a tab, and the \n makes a new line
        #the %s turns whatever is corresponding to that particular argument into a string
        #the two underscores before and after the name make this method "hidden", and it runs automatically when object is instantiated
        return 'SGD ID:\t%s\nFeature type:\t%s\nStatus: \t%s\nFeature name:\t%s\nCommon name:\t%s\nAlternative names: \t%s\nChromosome: \t%s\nDescription:\t%s\n\n' % (
        self.SGDID, self.feature_type, self.status, self.feature_name, self.common_name, self.alt_names, self.chromosome, self.description)
    
    
class idconverter(object):
    def __init__(self):
        current_directory = path.dirname(__file__)
        filename = 'SGD_features.tab'
        data_folder= 'Data'
        data_filepath= path.join(data_folder, filename)
        SGDfilepath= path.join(current_directory,data_filepath)
        f = open(SGDfilepath,'r')
        c = csv.reader(f, dialect="excel-tab")
        self.SGDdict = {}
        self.feature_name_map= {}
        self.common_name_map= {}
        self.alt_name_map= {}
        for row in c: 
            #the following parts of the tab file are parsed and passed into the SGD_gene object as arguments to build up the various parts
            idinf = SGD_gene(row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[15])
            feature_type = row[1].strip()
            if feature_type == "ORF" or feature_type=="tRNA_gene":
                self.SGDdict[idinf.SGDID] = idinf
                if idinf.feature_name != "":
                    self.feature_name_map[idinf.feature_name]= idinf.SGDID
                if idinf.common_name != "":
                    self.common_name_map[idinf.common_name]= idinf.SGDID
                for entry in idinf.alt_names:
                    if entry != "":
                        self.alt_name_map[entry]= idinf.SGDID
           
                
            
    def getgene(self, genid):
        genid = genid.strip().upper()
        gene = []
        if genid in self.SGDdict:
            gene = self.SGDdict[genid]
        elif genid in self.common_name_map:
            gene = self.SGDdict[self.common_name_map[genid]] 
        elif genid in self.feature_name_map:
            gene = self.SGDdict[self.feature_name_map[genid]] 
        elif genid in self.alt_name_map:
            gene = self.SGDdict[self.alt_name_map[genid]] 

        if type(gene)!=SGD_gene:
            raise Exception("The given gene {} is not present in SGD_features.tab. Sorry.".format(genid))
        else:
            return gene
            
