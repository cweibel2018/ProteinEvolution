# -*- coding: utf-8 -*-
"""

Created on Wed Jan pain 18:21:22 2018

@author: Catarina


Note the syntax is built for python 2
"""

import sys, json, os, csv, time, copy, math
#from collections import Counter
#from pathlib import Path
#from Bio import SeqIO
from Bio.Seq import Seq
#from Bio.SeqRecord import SeqRecord
#from Bio.Alphabet import generic_protein


def errorgenelistbuilder(errorfilename):
    errorgenelist =[]
    
    with open(errorfilename, 'r') as f :  
    
        reader = csv.reader(f, delimiter = ',')
        for row in reader :
            
            GeneID = row[0]
            
            if GeneID not in errorgenelist:
                errorgenelist.append(GeneID)
            else:
                print('huh')
                
    return errorgenelist

def errorgenefilegenerator(errorfilename, file):
    
    g = errorgenelistbuilder(errorfilename)
    
    with open(file, 'r') as f :  
    
        reader = csv.reader(f, delimiter = ',')
        for row in reader :
            
            #matching data with export of data from Navicat
            UID = row[0]
            SpeciesUID = row[1]
            Newick_Species = row[2]
            GeneID = row[3]
            ProteinID = row[4]
            CodingSequence = row[5]
            
            if GeneID in g:
                print(GeneID)
                
                output = [UID, SpeciesUID, Newick_Species, GeneID, ProteinID]
                #outputc = [CodingSequence]
                            
                c = Seq(CodingSequence)
                d = Seq.translate(c)
                            
                with open("%s.txt"%GeneID, "a") as f:
                    f.write(','.join(output[0:]) + ',' + ''.join(d[0:]) )
                f.close()
                
            else:
                pass
            
    print('done')
    
def bestsequencebylength(filename):
    
    sequencelist = []
    with open(filename, 'r') as f :  
    
        reader = csv.reader(f, delimiter = ',')
        for row in reader :
            
            #matching data with export of data from Navicat
            UID = row[0]
            SpeciesUID = row[1]
            Newick_Species = row[2]
            GeneID = row[3]
            ProteinID = row[4]
            ProteinSequence = row[5]
            
            if ProteinSequence not in sequencelist:
                sequencelist.append(ProteinSequence)
                
        c = max(sequencelist, key=len)
        
    return c

def bestsequencefulldata(filename,outfile):
    
    bestsequence = bestsequencebylength(filename) 
    
    with open(filename, 'r') as f :  
    
        reader = csv.reader(f, delimiter = ',')
        for row in reader :
            
            #matching data with export of data from Navicat
            UID = row[0]
            SpeciesUID = row[1]
            Newick_Species = row[2]
            GeneID = row[3]
            ProteinID = row[4]
            ProteinSequence = row[5]
            
            if ProteinSequence == bestsequence:
                output = [UID,SpeciesUID,Newick_Species,GeneID,ProteinID, ProteinSequence]
                
                with open(outfile, 'a') as f1:
                    f1.write(','.join(output[0:])+ '\n')
                    
                f1.close()
            else:
                print('nope %s'%GeneID)


def garbagedisposal(file):

	os.system('rm %s'%file)
	
	print('garbage disposal success for %s'%file)

                    
                    
def longestsequencepicker(errorfilename,originalfile,speciesUID):
    genelist = errorgenelistbuilder(errorfilename)
    print(genelist)
    
    errorgenefilegenerator(errorfilename, originalfile)
    
    genelistsize = len(genelist)
    print(genelistsize)
    
    for i in range (0, genelistsize):
        print(genelist[i])
        bestsequencefulldata("%s.txt"%genelist[i],"%s_longest.txt"%speciesUID) 
	garbagedisposal("%s.txt"%genelist[i])
print('done')  
            
longestsequencepicker("errors_159.txt","Outputnotuniquetrial159.txt",159)
        

    
            