# -*- coding: utf-8 -*-
"""
Created on Tues Oct 23 16:21:22 2018

@author: Catarina


This script takes species data (WITH SEQUENCES) in csv file form and returns two files in identical csv formats except one file (listed as "unique") 
contains only the data with unique gene:protein pairs, while the other (listed as "notunique") contains the data with gene: [protein 1, protein2, ...protein n]
sets.

This script has two functions:
    UniqueGeneIDFilterFrequency - parses and sorts thte gene:protein pairings into two dictionaries
    Uniquedatareproduction - reparses and sorts the full sequence data according to the dictionaries into two csv files
    
The end bit:
    
    for i in range(n,m):

	print(i)

	Unique_data_reproduction('%s.txt' %i)

	print('mission success for species %s'%i)
    
    this part basically iterates the the functions through the exported MySQL csv files, labeled by species ID
    
This script does parse redundantly, but sorts a 150,00 kb file in about 3 seconds.
"""

import sys, json, os, csv, time, copy, math
from collections import Counter
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein





def UniqueGeneIDFilterFrequency(filename):
  ##this function reads data from the table file and sorts it into two dictionaries: uniquedict, which contains
## the GeneIDs  (key) with only one protein ID (value) each, and notuniquedict, which contains GeneIDs with multiple protein IDs
    
    dict= {}
    
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
                if GeneID not in dict:
                    dict[GeneID] = [ProteinID]
                    
                else:
                    dict[GeneID].append(ProteinID)
                    
        dict = {k: list(set(v)) for k, v in dict.items()}
    
    notuniquedict ={}
    uniquedict = {}
    
    for k,v in dict.items() :
        values = list(v)
        keys = list(k)
        
        #filter for uniqueness
        if len(values) >1 :
            notuniquedict[k]= v
            
        else:
            uniquedict[k]=v
     
    return uniquedict       
   

def Unique_data_reproduction(filename):
    
    ##This function uses the dictionaries as a filter to sort all of the data into two files: Unique Gene IDs and not Unique gene IDs
    
    dictionary = UniqueGeneIDFilterFrequency(filename)
    
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
               
                
                if GeneID not in dictionary.keys():
                    notuniqueoutput = [UID,SpeciesUID,Newick_Species,GeneID,ProteinID,CodingSequence]
                    #print(SpeciesUID,GeneID, ProteinID)
                    
                    with open("Outputnotuniquetrial%s.txt" %i, "a") as f2:#append instead of overwrite "w"
                        f2.write(','.join(notuniqueoutput[0:]) + '\n')
                
                #elif EnsemblGeneID not in dictionary.keys():
                else:
                    uniqueoutput = [UID,SpeciesUID,Newick_Species,GeneID,ProteinID,CodingSequence]
                    #print(UID, GeneID, ProteinID)
                    with open("Outputuniquetrial%s.txt" %i, "a") as f3:#append instead of overwrite "w"
                        f3.write(','.join(uniqueoutput[0:]) + '\n')
                        
        print('done')                
        #f2.close()
        #f3.close()
print('done')

for i in range(413,524):

	print(i)

	Unique_data_reproduction('%s.txt' %i)

	print('mission success for species %s'%i)

