import Brongus.SGDIDConvertermodule as SGDIDConvertermodule
import re


idconverter = SGDIDConvertermodule.idconverter()

def change_keys_SGDID(dict_sequences):
    for gene_name in dict_sequences.keys():
        try:
            new_key= idconverter.getgene(gene_name).SGDID
            dict_sequences[new_key]= dict_sequences.pop(gene_name)
        except Exception:
            try:
                gene_names= re.split("_", gene_name)
                gene_name1= gene_names[0]
                gene_name2= gene_names[1]
                new_key= idconverter.getgene(gene_name1).SGDID
                dict_sequences[new_key]= dict_sequences.pop(gene_name)
            except Exception:
                try:
                    new_key= idconverter.getgene(gene_name2).SGDID
                    dict_sequences[new_key]= dict_sequences.pop(gene_name)
                except Exception:    
                    dict_sequences.pop(gene_name, None)
                    continue
    return dict_sequences