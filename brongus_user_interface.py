#!/usr/bin/env python
##user_interface
import argparse
import Brongus.motifscorecompiler as msc
import Brongus.motifscorereader as msr
import Brongus.common_assembler as ca
import sys
import Brongus.SGDIDConvertermodule as sgdc

if __name__ == "__main__":
    Brongusparser = argparse.ArgumentParser()
    Brongusparser.description= "After some genes or transcription factor names of S. cerevisae are entered in, data files containing the Position Weight Matrices of the transcription factor motifs will be scanned against the promoter sequences of the genes in order to generate a list of most likely gene targets or transcription factors will be generated."
    
    
    Brongusparser.add_argument("-fn","--filename", type= str, 
    help= "If you entered multiple genes or transcription factors, the please enter the name for the csv file you will generate.")
    
    Brongusparser.add_argument("-threshold", type= float, default= 0.0, 
    help= "Changing the threshold will change which log-odd scores are included in the array of hits, i.e. which associations between transcription factor and genetic sequence are deemed to be hits. By default, this is set to zero")
    
    
    Brongusparser.add_argument("-genetic_sequences", type= str, default= None,
    help= "This argument is not to be changed unless the user specifically wants to run the transcription factor motifs against a new assembly of genetic sequences (which must be in a .fasta format and in the Data folder of Brongus). Otherwise, the file promoter_sequences.fasta is always used as the assembly of genetic sequences.")
    
    ###only one of the arguments for -tf and -g can be run at the same. So they must be added to a mutually exclusive group
    tf_or_gene_group= Brongusparser.add_mutually_exclusive_group(required= True)
    
    tf_or_gene_group.add_argument("-tf",
    "--transcription_factors", 
    nargs= '*',
    type= str,
    help= "Please enter the transcription factor or factors whose top gene targets you wish to analyze. Please write the name of the transcription factor with no blanks and in capital letters.")
    
    tf_or_gene_group.add_argument("-g", 
    "--genes", 
    nargs= '*',
    type= str,
    help= "Please enter the gene or genes in order to analyze the transacription factors most likely to target them. Please write the genes with no blanks and in capital letters.")
    
    tf_or_gene_group.add_argument("-tf_fn",
    "--transcription_factor_filename",
    type= str,
    help= "It is also possible to upload a .txt file to the working directory that contains the list of transcription factors to be analyzed. Please format it so that there is one transcription factor per line.")
    
    
    tf_or_gene_group.add_argument("-g_fn",
    "--gene_filename",
    type= str,
    help= "It is also possible to upload a .txt file to the working directory that contains the list of genes to be analyzed. Please format it so that there is one gene per line.")
    
    
    args = Brongusparser.parse_args()
    
    ##run the proper functions and objects down below with arguments passed into them
    
    
    m= msc.MotifDict().motif_dict
    promoter_sequences= msc.PromoterSequences(args.genetic_sequences).promoter_sequences
    
    #tf pipeline
    if args.transcription_factors:
        tf_list= args.transcription_factors
        tf_list= list(tf_list)
        tf_ids= [sgdc.idconverter().getgene(tf) for tf in tf_list]
        if len(tf_list)==1:
            args.filename = None
            print("Now carrying out "+str(2)+" analyses for the transcription factor given... Please be patient.")
            sys.stdout.flush()
            msc.TFCompiledScores(m,promoter_sequences).compile_scores(tf_list[0], args.threshold)
            msr.tf_pickle_reader(tf_list[0])
        else:
            if args.filename is None:
                raise argparse.ArgumentTypeError("You must enter a filename if you provide more than one gene or transcription factor.")
            tf_processes= (len(tf_list)*2)+1
            print("Now carrying out "+str(tf_processes)+" analyses for the "+str(len(tf_list))+" transcription factors given... Please be patient.")
            sys.stdout.flush()
            for tf in tf_list:
                test3= msc.TFCompiledScores(m,promoter_sequences).compile_scores(tf, args.threshold)
                test4= msr.tf_pickle_reader(tf)
            test6= ca.common_tf_pickle_assembler(tf_list, args.filename)    
    
    #gene pipeline
    if args.genes:
        gene_list= args.genes
        gene_list= list(gene_list)
        if len(gene_list)== 1:
            args.filename = None
            print("Now carrying out "+str(2)+" analyses for the gene given... Please be patient.")
            sys.stdout.flush()
            test7= msc.GeneCompiledScores(m, promoter_sequences).compile_scores(gene_list[0], args.threshold)
            test8= msr.gene_pickle_reader(gene_list[0])
        else:
            if args.filename is None:
                raise argparse.ArgumentTypeError("You must enter a filename if you provide more than one gene or transcription factor.")
            gene_processes= (len(gene_list)*2)+1
            print("Now carrying out "+str(gene_processes)+" analyses for the "+str(len(gene_list))+" genes given... Please be patient.")
            sys.stdout.flush()
            for gene in gene_list:
                test9= msc.GeneCompiledScores(m,promoter_sequences).compile_scores(gene, args.threshold)
                test10= msr.gene_pickle_reader(gene)
            test11= ca.common_gene_pickle_assembler(gene_list, args.filename)        
    
    ##tf filename pipeline
    if args.transcription_factor_filename:
        tfs= open(args.transcription_factor_filename, 'r')
        tf_list= []
        for line in tfs:
            tf_list+=[line]
        if len(tf_list)==1:
            args.filename = None
            tf_processes= 2
            print("Now carrying out "+str(tf_processes)+" analyses for the "+str(len(tf_list))+" transcription factors given... Please be patient.")
            sys.stdout.flush()
            test3= msc.TFCompiledScores(m,promoter_sequences).compile_scores(tf_list[0], args.threshold)
            test4= msr.tf_pickle_reader(tf_list[0])
        else:
            if args.filename is None:
                raise argparse.ArgumentTypeError("You must enter a filename if you provide more than one gene or transcription factor.")
            tf_processes= (len(tf_list)*2)+1
            print("Now carrying out "+str(tf_processes)+" analyses for the "+str(len(tf_list))+" transcription factors given... Please be patient.")
            sys.stdout.flush()
            for tf in tf_list:
                test3= msc.TFCompiledScores(m,promoter_sequences).compile_scores(tf, args.threshold)
                test4= msr.tf_pickle_reader(tf)
            test6= ca.common_tf_pickle_assembler(tf_list, args.filename) 
    
    ##gene filename pipeline
    if args.gene_filename:
        genes= open(args.gene_filename, 'r')
        gene_list= []
        for line in genes:
            gene_list += [line]
        if len(gene_list)==1:
            args.filename = None
            gene_processes= 2
            print("Now carrying out "+str(gene_processes)+" analyses for the "+str(len(gene_list))+" genes given... Please be patient.")
            sys.stdout.flush()
            test20= msc.GeneCompiledScores(m,promoter_sequences).compile_scores(gene_list[0], args.threshold)
            test10= msr.gene_pickle_reader(gene_list[0])
        else:
            if args.filename is None:
                raise argparse.ArgumentTypeError("You must enter a filename if you provide more than one gene or transcription factor.")
            gene_processes= (len(gene_list)*2)+1
            print("Now carrying out "+str(gene_processes)+" analyses for the "+str(len(gene_list))+" genes given... Please be patient.")
            sys.stdout.flush()
            for gene in gene_list:
                test9= msc.GeneCompiledScores(m,promoter_sequences).compile_scores(gene, args.threshold)
                test10= msr.gene_pickle_reader(gene)
            test11= ca.common_gene_pickle_assembler(gene_list, args.filename) 
    
    
