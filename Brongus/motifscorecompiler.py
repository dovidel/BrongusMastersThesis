from Bio import SeqIO
from Bio.Alphabet.IUPAC import IUPACUnambiguousDNA
from Bio import motifs
import csv
import pickle
import numpy as np
import os
from os import path
import SGDIDConvertermodule as sgdc
from tqdm import tqdm
import change_keys_SGDID as cks

current_directory = path.dirname(__file__)
idconverter = sgdc.idconverter()
SGD_dict= idconverter.SGDdict

class PromoterSequences(object):
    def __init__(self, promoter_sequences_filename= None):
        if promoter_sequences_filename is None:
            promoter_sequences_filename= "promoter_sequences.fasta"   
            fullfile = path.join('Data', promoter_sequences_filename)
            self.promoter_sequences= SeqIO.to_dict(SeqIO.parse(fullfile, "fasta", alphabet=IUPACUnambiguousDNA()))
            self.promoter_sequences= cks.change_keys_SGDID(self.promoter_sequences)
        else:
            self.promoter_sequences= SeqIO.to_dict(SeqIO.parse(promoter_sequences_filename, "fasta", alphabet=IUPACUnambiguousDNA()))
            self.promoter_sequences= cks.change_keys_SGDID(self.promoter_sequences)

    def __getitem__(self, key):
        return self.promoter_sequences[key]



class MotifDict(object):
    def __init__(self,matrixfile=None,pfmdir=None):
        annotationmatrixfile= 'SaccCereAnnotation.txt'
        matrixfile='SaccCereTFMATRIX.txt'
        annotationmatrixfile= path.join('Data', annotationmatrixfile)
        matrixfile= path.join('Data', matrixfile)
        TFmatrixFile= path.join(current_directory, matrixfile)
        annotationfile= path.join(current_directory, annotationmatrixfile)
        medline= open(annotationfile, 'r')
        medlinematrix = csv.reader(medline, dialect='excel-tab')
        self.medline_dict= {}
        for row in medlinematrix:
                    if row[1] == 'medline':
                        self.medline_dict[row[0]]= row[2]  
        if pfmdir is None:
            pfmdir="SaccCerePFMFlatFileDir"
            pfm_folder= path.join(current_directory,'Data', pfmdir)
        else:
            pfm_folder=pfmdir
        TFmatrixReader= open(TFmatrixFile, 'r')
        TFmatrix = csv.reader(TFmatrixReader, dialect='excel-tab')
        self.motif_dict = dict()
        for row in TFmatrix:
            tf_accession= row[0]
            common_name= row[4]
            tf_pfm= "{}.{}.pfm".format(row[2],row[3])
            filename= path.join('Data',pfm_folder, tf_pfm)
            with open(filename, 'r') as handle:
                self.motif_dict[common_name] = motifs.read(handle,'pfm')
                try:
                    sgdid= idconverter.getgene(common_name).SGDID
                    self.motif_dict[sgdid]= self.motif_dict.pop(common_name)
                    self.motif_dict[sgdid].medline= self.medline_dict[tf_accession]
                except Exception:
                    self.motif_dict.pop(common_name, None)
                    continue 
                         
        
    def __getitem__(self, key):
        return self.motif_dict[key]
        


class TFCompiledScores(object):
    def __init__(self, motif_dict,promoter_sequences):
        self.motif_dict = MotifDict
        self.promoter_sequences= PromoterSequences

    def compile_scores(self, tf_id, threshold):
        tf_SGDID = idconverter.getgene(tf_id).SGDID
        output_folder= "Output"
        if not os.path.isdir(output_folder):
            os.mkdir(output_folder)
        pickle_dir= "Pickles"
        pickle_folder= path.join(output_folder,pickle_dir)
        if not os.path.isdir(pickle_folder):
            os.mkdir(pickle_folder)
        motif= MotifDict()
        promoter_sequences= PromoterSequences().promoter_sequences
        tf_score_dict = {}
        try:
            motif= motif.motif_dict[tf_SGDID]
        except KeyError:
            print("The transcription factor "+str(tf_id)+ " was not found in the current dictionary of yeast transcription factor motifs.")
        motif = motif.pssm
        reverse_motif= motif.reverse_complement()
        for gene_name in tqdm(promoter_sequences.keys()):
            promoter= promoter_sequences[gene_name]
            promoter_seq= promoter.seq
            gene_SGDID= idconverter.getgene(gene_name).SGDID
            try:
                forward_scores= motif.calculate(promoter_seq)
                reverse_scores= reverse_motif.calculate(promoter_seq)
            except MemoryError:
                raise "Gene error:{}".format(gene_name)
            forward_and_reverse_scores= np.concatenate([forward_scores,reverse_scores])
            forward_and_reverse_scores= np.sort(forward_and_reverse_scores)
            final_forward_and_reverse_scores= [score for score in forward_and_reverse_scores if score> threshold]
            if len(final_forward_and_reverse_scores)>0:
                tf_score_dict[gene_SGDID]= final_forward_and_reverse_scores
        pickle_file = 'tf_{}_scores.pickle'.format(tf_SGDID)
        folder_location= path.join(pickle_folder, pickle_file)
        pickle_file= open(folder_location, 'wb')
        pickle.dump(tf_score_dict, pickle_file)
        pickle_file.close()
        return(tf_score_dict)
    
    
class GeneCompiledScores(object):
    def __init__(self, motif_dict,promoter_sequences):
        self.motif_dict = MotifDict
        self.promoter_sequences= PromoterSequences
    
    def compile_scores(self,gene_id, threshold):
        gene_SGDID= idconverter.getgene(gene_id).SGDID
        pickle_dir= "Pickles"
        output_location= 'Output'
        if not path.isdir(output_location):
            os.mkdir(output_location)
        pickle_location= path.join(output_location, pickle_dir)
        if not path.isdir(pickle_location):
            os.mkdir(pickle_location)
        motif= MotifDict()
        motif_dict= motif.motif_dict
        promoter_sequences= PromoterSequences().promoter_sequences
        promoter= promoter_sequences[gene_SGDID]
        promoter_seq= promoter.seq
        gene_score_dict = {}
        for transcription_factor in tqdm(motif_dict.keys()):
            motif= motif_dict[transcription_factor].pssm
            reverse_motif = motif.reverse_complement()
            try:
                forward_scores= motif.calculate(promoter_seq)
                reverse_scores= reverse_motif.calculate(promoter_seq)
            except MemoryError:
                raise "Gene error:{}".format(transcription_factor)
            forward_and_reverse_scores= np.concatenate([forward_scores,reverse_scores])
            forward_and_reverse_scores= np.sort(forward_and_reverse_scores)
            final_forward_and_reverse_scores= [score for score in forward_and_reverse_scores if score> threshold]
            if len(final_forward_and_reverse_scores)>0:
                gene_score_dict[transcription_factor]= (final_forward_and_reverse_scores)
        pickle_file = 'gene_{}_scores.pickle'.format(gene_SGDID)
        folder_location= path.join(pickle_location, pickle_file)
        pickle_file= open(folder_location, 'wb')
        pickle.dump(gene_score_dict, pickle_file)
        pickle_file.close()
        return(gene_score_dict)