###once the pickles of the input genes or transcription factors are built up, they must be accessed, then the results must be displayed for the user
###in the case of multiple inputs, the results must be combined together in order to see which hits are shared between the input genes or tfs

import matplotlib.pyplot as plt
import pickle
from os import path
import os
import numpy as np
import Brongus.SGDIDConvertermodule as sgdc
from pandas import DataFrame
from tqdm import tqdm
import Brongus.motifscorecompiler as motifscorecompiler


current_directory = path.dirname(__file__)
idconverter= sgdc.idconverter()

def tf_pickle_reader(tf_id):
    tf= idconverter.getgene(tf_id).SGDID
    output_folder= "Output"
    if not os.path.isdir(output_folder):
        raise IOError
    pickle_dir= "Pickles"
    tf_histogram_folder= "TF_Histograms"
    tf_csv_folder= "TF_CSV"
    pickle_folder= path.join(output_folder,pickle_dir)
    if not os.path.isdir(pickle_folder):
        raise IOError
    tf_histogram_folder_location= path.join(output_folder,tf_histogram_folder)
    if not os.path.isdir(tf_histogram_folder_location):
        os.mkdir(tf_histogram_folder_location)
    tf_csv_folder_location= path.join(output_folder, tf_csv_folder)
    if not os.path.isdir(tf_csv_folder_location):
        os.mkdir(tf_csv_folder_location)
    folder_contents = os.listdir(pickle_folder)
    scores_distribution= []
    if 'tf_{}_scores.pickle'.format(tf) in folder_contents:
        tf_pickle_file = open(path.join(pickle_folder,'tf_{}_scores.pickle'.format(tf)),'rb')
        tf_dict = pickle.load(tf_pickle_file)
        tf_pickle_file.close()
        df = DataFrame(columns=['Gene Feature Name','Common Name','Number of Hits','Total Sum of Scores','Avg. Score per Hit','List of Scores','Gene Description'])
        tf_common_name= idconverter.getgene(tf).common_name
        tf_csv= path.join(tf_csv_folder_location, "tf_{}.csv".format(tf_id))
        gene_ids= []
        for gene_id in tqdm(tf_dict.keys()):
            list_of_scores= []
            gene_description= idconverter.getgene(gene_id).description
            gene_common_name= idconverter.getgene(gene_id).common_name
            gene_feature_name= idconverter.getgene(gene_id).feature_name
            total_sum_scores= sum(tf_dict[gene_id])
            number_binding_sites= len(np.atleast_1d(tf_dict[gene_id]))
            average_score_per_binding_site= total_sum_scores/number_binding_sites
            for score in tf_dict[gene_id]:
                list_of_scores+= [score]
                scores_distribution += [score]
            df.loc[gene_id]= [gene_feature_name, gene_common_name, number_binding_sites, total_sum_scores, average_score_per_binding_site, list_of_scores, gene_description]
            gene_ids += [gene_id]
        df= df.sort_values(['Total Sum of Scores'], ascending=False)
        df.to_csv(tf_csv)
        tf_histogram_filename= path.join(tf_histogram_folder_location, tf+"_"+tf_id)
        plt.cla()
        plt.hist(scores_distribution, bins=30)
        plt.xlabel('Log-Odd Scores of {} Binding Sites'.format(tf_common_name))
        plt.ylabel('Number of Potential Binding Sites')
        plt.suptitle("Distribution of Binding Site Scores for the transcription factor {},{} \nagainst the promoter regions of yeast genes".format(tf_id,tf_common_name), fontsize= 12)
        plt.savefig(tf_histogram_filename)
    else:
        raise IOError


def gene_pickle_reader(gene_id):
    gene= idconverter.getgene(gene_id).SGDID
    output_folder= "Output"
    pickle_dir= "Pickles"
    gene_histogram_folder= "Gene_Histograms"
    gene_csv_folder= "Gene_CSV"
    if not os.path.isdir(output_folder):
        os.mkdir(output_folder)
    pickle_folder= path.join(output_folder, pickle_dir)
    if not os.path.isdir(pickle_folder):
        raise IOError
    gene_histogram_folder_location= path.join(output_folder, gene_histogram_folder)
    if not os.path.isdir(gene_histogram_folder_location):
        os.mkdir(gene_histogram_folder_location)
    gene_csv_folder_location= path.join(output_folder, gene_csv_folder)
    if not os.path.isdir(gene_csv_folder_location):
        os.mkdir(gene_csv_folder_location)
    pickle_folder_contents = os.listdir(pickle_folder)
    scores_distribution= []
    if 'gene_{}_scores.pickle'.format(gene) in pickle_folder_contents:
        df = DataFrame(columns=['TF Feature Name','TF Common Name','Number of Hits','Total Sum of Scores','Avg. Score per Hit','List of Scores','Medline','TF Description'])
        gene_pickle_file = open(path.join(pickle_folder,'gene_{}_scores.pickle'.format(gene)),'rb')
        gene_dict = pickle.load(gene_pickle_file)
        gene_pickle_file.close()
        gene_common_name= idconverter.getgene(gene).common_name
        motif= motifscorecompiler.MotifDict()
        gene_csv= path.join(gene_csv_folder_location, "gene_{}.csv".format(gene_id))
        for tf_id in tqdm(gene_dict.keys()):
            tf_feature_name= idconverter.getgene(tf_id).feature_name
            tf_common_name= idconverter.getgene(tf_id).common_name
            tf_description= idconverter.getgene(tf_id).description
            tf_medline= motif[tf_id].medline
            tf_medline_url= 'www.pubmed.com/{}'.format(tf_medline)
            total_sum_scores= sum(gene_dict[tf_id])
            number_binding_sites= len(np.atleast_1d(gene_dict[tf_id]))
            average_score_per_binding_site= total_sum_scores/number_binding_sites
            list_of_scores=[]
            for score in gene_dict[tf_id]:
                list_of_scores+= [score]
                scores_distribution += [score]
            df.loc[tf_id]= [tf_feature_name, tf_common_name, number_binding_sites, total_sum_scores, average_score_per_binding_site, list_of_scores,tf_medline_url,tf_description]
        df= df.sort_values(['Total Sum of Scores'], ascending=False) 
        df.to_csv(gene_csv)
        gene_histogram_filename= path.join(gene_histogram_folder_location,gene+"_"+gene_id)
        plt.cla()
        plt.hist(scores_distribution, bins=30)
        plt.xlabel('Log-Odd Scores of {} Binding Sites'.format(gene_common_name))
        plt.ylabel('Number of Potential Binding Sites')
        plt.suptitle("Distribution of Binding Site Scores for the gene {},{} \nagainst the transcription factors of yeast genes".format(gene_id, gene), fontsize= 12)
        plt.savefig(gene_histogram_filename)
    else:
        raise IOError
