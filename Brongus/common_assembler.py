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

current_directory = path.dirname(__file__)
idconverter= sgdc.idconverter()
motif= motifscorecompiler.MotifDict()

def common_tf_pickle_assembler(tf_id_list, filename):
    pickle_dir= "Pickles"
    output_folder= 'Output'
    pickle_folder= path.join(output_folder, pickle_dir)
    if not os.path.isdir(pickle_folder):
        raise IOError
    folder_contents = os.listdir(pickle_folder)
    tf_csv_folder= "TF_CSV"
    tf_csv_folder_location= path.join(output_folder, tf_csv_folder)
    if not os.path.isdir(tf_csv_folder_location):
        raise IOError
    ColumnList=[]
    gene_id_dict= {}
    total_genes_score_dict= {}
    for tf_id in tf_id_list:
        try:
            tf= idconverter.getgene(tf_id).SGDID
        except Exception:
            invalid_tf= tf_id_list.index(tf_id)
            del tf_id_list[invalid_tf]
            warn("The transcripton factors {} entered could not be located in the SGD database. Please check the name.".format(tf_id))
        ColumnList+=['{} Number of Hits'.format(tf_id),'{} Total Sum of Scores'.format(tf_id)]
        if 'tf_{}_scores.pickle'.format(tf) in folder_contents: 
            tf_pickle_file = open(path.join(pickle_folder,'tf_{}_scores.pickle'.format(tf)),'rb')
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
        gene_info+=(gene_feature_name, gene_common_name)
        for tf_id in total_genes_score_dict.keys():
            if gene_id in total_genes_score_dict[tf_id].keys(): 
                gene_info+= (total_genes_score_dict[tf_id][gene_id][0],)
                gene_info+= (total_genes_score_dict[tf_id][gene_id][1],)
        gene_info+= (gene_description,)
        data+= [gene_info]
    df= DataFrame.from_records(data, index= common_genes, columns= ['Gene Feature Name','Gene Common Name']+ColumnList+['Gene Description'])
    tf_csv= path.join(tf_csv_folder_location, "{}.csv".format(filename))
    df.to_csv(tf_csv) 


def common_gene_pickle_assembler(gene_id_list, filename):
    pickle_dir= "Pickles"
    output_folder= 'Output'
    pickle_folder= path.join(output_folder, pickle_dir)
    if not path.isdir(output_folder):
        os.mkdir(output_folder)
    if not path.isdir(pickle_folder):
        raise IOError
    folder_contents = os.listdir(pickle_folder)
    gene_csv_dir= "Gene_CSV"
    gene_csv_folder= path.join(output_folder, gene_csv_dir)
    if not path.isdir(gene_csv_folder):
        os.mkdir(gene_csv_folder)
    ColumnList=[]
    tf_id_dict= {}
    total_genes_score_dict= {}
    for gene_id in gene_id_list:
        gene= idconverter.getgene(gene_id).SGDID
        ColumnList+=['{} Number of Hits'.format(gene_id),'{} Total Sum of Scores'.format(gene_id)]
        if 'gene_{}_scores.pickle'.format(gene) in folder_contents: 
            gene_pickle_file = open(path.join(pickle_folder,'gene_{}_scores.pickle'.format(gene)),'rb')
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
    #tf_sets= (set(tf_ids) for tf_ids in tf_id_dict.values())
    common_tfs= set.intersection(*tf_sets)
    common_tfs= list(common_tfs)
    data = []
    for tf_id in tqdm(common_tfs):
        tf_info= ()
        tf_description= idconverter.getgene(tf_id).description 
        tf_common_name= idconverter.getgene(tf_id).common_name
        tf_feature_name= idconverter.getgene(tf_id).feature_name
        tf_medline= motif[tf_id].medline
        tf_medline_url= 'www.pubmed.com/{}'.format(tf_medline)
        tf_info+= (tf_feature_name, tf_common_name)
        for gene_id in total_genes_score_dict.keys():
            if tf_id in total_genes_score_dict[gene_id].keys(): 
                tf_info+= (total_genes_score_dict[gene_id][tf_id][0],)
                tf_info+= (total_genes_score_dict[gene_id][tf_id][1],)
        tf_info+= (tf_medline_url, tf_description)
        data+= [tf_info]
    df= DataFrame.from_records(data, index= common_tfs, columns= ['TF Feature Name','TF Common Name']+ColumnList+['Medline']+['TF Description'])
    genes_csv= path.join(gene_csv_folder, "{}.csv".format(filename))
    df.to_csv(genes_csv) 
