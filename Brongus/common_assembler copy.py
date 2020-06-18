###this submodule will look at multiple pickles at once and intersect
from os import path
import os
import pickle
import Brongus.SGDIDConvertermodule as sgdc
import numpy as np
from pandas import DataFrame
from warnings import warn
import Brongus.motifscorecompiler as motifscorecompiler
from tqdm import tqdm
import sys

current_directory = path.dirname(__file__)
idconverter= sgdc.idconverter()
motif= motifscorecompiler.MotifDict()

def common_tf_pickle_assembler(tf_id_list, filename):
    score_list_folder= "TFbinding"
    folder_location= path.join(current_directory,'Output', score_list_folder)
    folder_contents = os.listdir(folder_location)
    tf_csv_folder= "TF_CSV"
    tf_csv_folder_location= path.join(current_directory,'Output', tf_csv_folder)
    ColumnList=[]
    gene_id_dict= {}
    total_genes_score_dict= {}
    for tf_id in tf_id_list:
        try:
            tf= idconverter.getgene(tf_id).SGDID
        except Exception:
            invalid_tf= tf_id_list.index(tf_id)
            del tf_id_list[invalid_tf]
            warn ("The transcripton factors {} entered could not be located in the SGD database. Please check the name.".format(tf_id))
        ColumnList+=['{} Number of Hits'.format(tf_id),'{} Total Sum of Scores'.format(tf_id)]
        if 'tf_{}_score_list.pickle'.format(tf) in folder_contents: 
            tf_pickle_file = open(path.join(folder_location,'tf_{}_score_list.pickle'.format(tf)),'rb')
            tf_dict = pickle.load(tf_pickle_file)
            tf_pickle_file.close()
            gene_ids= []
            gene_score_dict= {}
            for gene_id in tf_dict.keys():
                total_sum_scores= sum(tf_dict[gene_id])
                number_binding_sites= len(np.atleast_1d(tf_dict[gene_id]))
                gene_score_dict[gene_id]= [number_binding_sites, total_sum_scores]
                gene_ids += [gene_id] 
            gene_id_dict[tf_id] = set(gene_ids)
            total_genes_score_dict[tf_id]= gene_score_dict      
    gene_sets= gene_id_dict.values()
    common_genes= set.intersection(*gene_sets)
    common_genes= list(common_genes)
    data= []
    for gene_id in tqdm(common_genes):
        gene_info= ()
        gene_description= idconverter.getgene(gene_id).description 
        gene_common_name= idconverter.getgene(gene_id).common_name
        gene_feature_name= idconverter.getgene(gene_id).feature_name
        gene_info += (gene_feature_name, gene_common_name)
        #df.set_value(gene_id,'Gene Feature Name',gene_feature_name)
        #df.set_value(gene_id,'Gene Common Name', gene_common_name)
        #df.set_value(gene_id,'Gene Description', gene_description)
        for tf_id in total_genes_score_dict.keys():
            if gene_id in total_genes_score_dict[tf_id].keys(): 
                gene_info += (total_genes_score_dict[tf_id][gene_id][0],)
                gene_info += (total_genes_score_dict[tf_id][gene_id][1],)
                #df.set_value(gene_id, '{} Number of Hits'.format(tf_id), total_genes_score_dict[tf_id][gene_id][0])
                #df.set_value(gene_id, '{} Total Sum of Scores'.format(tf_id), total_genes_score_dict[tf_id][gene_id][1])
        gene_info+= tuple(gene_description)
        data+= [gene_info]
    print("\n"+str(len(data)))
    df= DataFrame.from_records(data, columns= ['Gene Feature Name','Gene Common Name']+ColumnList+['Gene Description'])
    tf_csv_filename= filename
    tf_csv= path.join(tf_csv_folder_location, "{}.csv".format(tf_csv_filename))
    df.to_csv(tf_csv) 


def common_gene_pickle_assembler(gene_id_list, filename):
    score_list_folder= "TFbinding"
    folder_location= path.join(current_directory, score_list_folder)
    folder_contents = os.listdir(folder_location)
    tf_csv_folder= "Gene_CSV"
    tf_csv_folder_location= path.join(current_directory, tf_csv_folder)
    ColumnList=[]
    tf_id_dict= {}
    total_genes_score_dict= {}
    for gene_id in gene_id_list:
        gene= idconverter.getgene(gene_id).SGDID
        ColumnList+=['{} Number of Hits'.format(gene_id),'{} Total Sum of Scores'.format(gene_id)]
        if 'gene_{}_score_list.pickle'.format(gene) in folder_contents: 
            gene_pickle_file = open(path.join(folder_location,'gene_{}_score_list.pickle'.format(gene)),'rb')
            gene_dict = pickle.load(gene_pickle_file)
            gene_pickle_file.close()
            tf_ids= []
            tf_score_dict= {}
            for tf_id in gene_dict.keys():
                total_sum_scores= sum(gene_dict[tf_id])
                number_binding_sites= len(np.atleast_1d(gene_dict[tf_id]))
                tf_score_dict[tf_id]= [number_binding_sites, total_sum_scores]
                tf_ids += [tf_id] 
            tf_id_dict[gene_id] = set(tf_ids)
            total_genes_score_dict[gene_id]= tf_score_dict      
    tf_sets= tf_id_dict.values()
    common_tfs= set.intersection(*tf_sets)
    common_tfs= list(common_tfs)
    df= DataFrame(columns= ['TF Feature Name','TF Common Name']+ColumnList+['Medline']+['TF Description'])
    for tf_id in common_tfs:
        tf_description= idconverter.getgene(tf_id).description 
        tf_common_name= idconverter.getgene(tf_id).common_name
        tf_feature_name= idconverter.getgene(tf_id).feature_name
        tf_medline= motif[tf_id].medline
        tf_medline_url= 'www.pubmed.com/{}'.format(tf_medline)
        df.set_value(tf_id,'TF Feature Name',tf_feature_name)
        df.set_value(tf_id,'TF Common Name', tf_common_name)
        df.set_value(tf_id,'TF Description', tf_description)
        df.set_value(tf_id,'Medline', tf_medline_url)
        for gene_id in total_genes_score_dict.keys():
            if tf_id in total_genes_score_dict[gene_id].keys(): 
                df.set_value(tf_id, '{} Number of Hits'.format(gene_id), total_genes_score_dict[gene_id][tf_id][0])
                df.set_value(tf_id, '{} Total Sum of Scores'.format(gene_id), total_genes_score_dict[gene_id][tf_id][1])
    genes_csv_filename= filename
    genes_csv= path.join(tf_csv_folder_location, "{}.csv".format(genes_csv_filename))
    df.to_csv(genes_csv) 

#tf_id_list= ['ROX1', 'MSN2']
#test= common_tf_pickle_assembler(tf_id_list, 'pleasework')

gene_list= ['ZIP1', 'COX2']
test2= common_gene_pickle_assembler(gene_list, 'pleasework')