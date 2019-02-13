# -*- coding: utf-8 -*-
"""
Created on Wed Nov pain 18:21:22 2018

@author: Catarina

This script takes the output file "OutputnotuniqueXspeciesuid.txt" from Trialer.py, parses the file, and generates a fasta file per gene 
The output files are named by geneID, which is necessary for the script trialpipes2trial.py

This script runs by the command "$python TrialFastaSplitoriginalcp.py"

Note the syntax is built for python 2
"""

import sys, json, os, csv, time, copy, math
from collections import Counter
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein





def UniqueProteinIDFilterFrequency(filename):
  ##this function reads data from the table file of the NOT unique pairs and sorts it into a dictionary which contains the ProteinID as the key
## and the Protein Sequence as the value    
    proteindict= {}
	
    with open(filename, 'r') as f :  
    
        reader = csv.reader(f, delimiter = ',')
        for row in reader :
            
                	#matching data with export of data from Navicat
                UID = row[0]
                SpeciesUID = row[1]
                Newick_Species = row[2]
                GeneID = row[3]
                ProteinID = row[4]
                CodingSequence = row[5]
			               
               
                	#making two sets, one of Protein Tags, one of Sequences
                if ProteinID not in proteindict:
                    proteindict[ProteinID] = [CodingSequence]
                    
                else:
                    proteindict[ProteinID].append(CodingSequence)
                    
        proteindict = {k: list(set(v)) for k, v in proteindict.items()}
         
        return proteindict   

print('done')  

def UniqueGeneIDFilterFrequency(filename):
  ##this function reads data from the table file and sorts it into two dictionaries: uniquedict, which contains
## the GeneIDs  (key) with only one protein ID (value) each, and notuniquedict, which contains GeneIDs with multiple protein IDs
    genedict= {}
    
    with open(filename, 'r') as f :  
    
        reader = csv.reader(f, delimiter = ',')
        for row in reader :
            
                #matching data with export of data from Navicat
                UID = row[0]
                SpeciesUID = row[1]
                Newick_Species = row[2]
                GeneID = row[3]
                ProteinID = row[4]
                CodingSequence = row[5]
               
               
                #making two lists, one of Genes, one of Protein tags
                if GeneID not in genedict:
                    genedict[GeneID] = [ProteinID]
                    
                else:
                    genedict[GeneID].append(ProteinID)
                    
        genedict = {k: list(set(v)) for k, v in genedict.items()}
    
    notuniquegenedict ={}
    uniquegenedict = {}
    
    for k,v in genedict.items() :
        values = list(v)
        keys = list(k)
        
        #filter for uniqueness
        if len(values) >1 :
            notuniquegenedict[k]= v
            
        else:
            uniquegenedict[k]=v
     
    return notuniquegenedict 

print('done')
   

def Unique_data_reproduction(filename):
    
    ##This function uses the dictionaries as a filter to sort all of the data into two files: Unique Gene IDs and not Unique gene IDs
    genedictionary = UniqueGeneIDFilterFrequency(filename)
    proteindictionary = UniqueProteinIDFilterFrequency(filename)

    print(genedictionary.keys())

    
    with open(filename, 'r') as f :  
    
        	reader = csv.reader(f, delimiter = ',')
        	for row in reader :
            
                	#matching data with export of data from Navicat
                    UID = row[0]
                    SpeciesUID = row[1]
                    Newick_Species = row[2]
                    GeneID = row[3]
                    ProteinID = row[4]
                    CodingSequence = row[5]
                    
                    if GeneID in genedictionary.keys():
                        
                            print(GeneID)
                            
                            output = [UID, SpeciesUID, Newick_Species, GeneID, ProteinID]
                            outputc = [CodingSequence]
                            
                            c = Seq(CodingSequence)
                            d = Seq.translate(c)
                            
                            with open("%s.txt"%GeneID, "a") as f: #do NOT switch the GeneID for UID, else will have a file per UID instead of per gene
                                f.write('>')
                                f.write(','.join(output[0:]) + '\n')
                                f.write(''.join(d[0:]) + '\n' + '\n')
                                
                            f.close()
                    
                    else:
                        print('uh oh %s'%GeneID)
								
					
					
                             
    
print('done')


Unique_data_reproduction("Outputnotuniquetrial87.txt")


###Test module to check different outputs
def tester(filename):

    c = Unique_data_reproduction(filename)
    d = UniqueGeneIDFilterFrequency(filename)
    e = UniqueProteinIDFilterFrequency(filename)

	#print(d.keys())
	#print(d.values())
	#print(d)

    print(e.keys())
	#print(e.values())
	#print(e)

