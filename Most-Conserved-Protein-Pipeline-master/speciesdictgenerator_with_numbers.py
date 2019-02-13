# -*- coding: utf-8 -*-
"""
Created on Monday Oct 8 18:34:30 2018
@author: Catarina


This script takes the SpeciesList.txt file and a Newick formatted tree (with the species from SpeciesList) to pair the most closely related species in a dictionary.
This dictionary will then be used in the trialpipes.py script for the sequences of one species to be Blastp -ed against the database of the sister species. 
"""
from ete3 import Tree
import sys, json, os, csv, time, copy, math
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
def TreeDataComparison(filename):
   
    Specieslist = []
    nodename = []
    Truesearchlist =[]
	
	#parse the SpeciesList.txt
    with open(filename, 'r') as f : 
   
        reader = csv.reader(f, delimiter = ',')
        for row in reader :
           
                #matching data with export of data from Navicat
                SpeciesUID = row[0]
                Species_name = row[1]
                Newick_Formatted_Species = row[2]
                Species_commom_name = row[3]
                Highest_gold_status = row[4]
                Study_ID = row[5]
                Gene_count = row[6]
                Ensemble_accession = row[7]
                Ensembl_DB = row[8]
                Gene_build_method = row[9]
                TaxonID = row[10]
                inBiomart = row[11]
               
                #building a list of the species names with no repeats, to search against the tree
                if Newick_Formatted_Species in Specieslist:
                        pass
                       
                else:
                    Specieslist.append(Newick_Formatted_Species)
                   
        #Species list from previous output is saved into a list in Newick Format            
        c= Specieslist
       
        #read in the tree
        t = Tree("specieslistfortree2realtree.nwk",format=1)
       
        #turn the node names into a species list (including leaf/branch names in the form of numbers)
        leaf=t.get_leaf_names(is_leaf_fn=None)
       
        for x in leaf:
            nodename.append(x)
           
        for x in c:
            if x in nodename:
                Truesearchlist.append(x)
           
            else:
                pass
       
        return Truesearchlist


   
#function to get closest relative on tree 
def Sisterspecies(filename):
   
    	sisterlist =[]
    	truelist=[]
    	t = Tree("specieslistfortree2realtree.nwk",format=1)
   
    	d= TreeDataComparison(filename)
   
	# dictionary of all sisterspecies pairs
    	sisterspeciesdict = {} 

	# dictionary of species:sisterspecies with 1:1 ratio
	cleansisterdict = {}
   
    	#making a list of all nodes that are terminal (AKA specieslist from tree)
    	for node in t.traverse():
        	#node
        		if node.is_leaf() :
            			speciesname = node.name
            			sisterlist.append(speciesname)
        		else:
				pass
	
   
#loop to make a dictionary of all species/sisters  
#t.traverse() reads through the input tree in a specific order
	
	for node in t.traverse():

       
		if node.name in d:
           
			firstspecies = node.name
			#print(firstspecies)
           
            			#the next node up from the leaf
            		parent = node.up
           
            		#get all descendants of the node
            		parentree= parent.get_descendants(strategy='levelorder',is_leaf_fn=None)
           		#print(parentree)
           
            		#prints the node name if it is listed in the species list
			sisters = []
            		for node in parentree:
				
				#tests that all is reading correctly
                		if node.name in sisterlist:
					sisters.append(node.name)
				
				#removes the key from the list of future values 
				if firstspecies in sisters:
					#print('welp')
					sisters.remove(firstspecies)
					#print(sisters)
				else:
					pass
					#print('huh')
			#print(sisters)	
			
			#NOTE: dictionary objects don't have '.append' attributes	
			sisterspeciesdict[firstspecies]= sisters
			cleansisterdict[firstspecies] = sisters[0]
		
	return cleansisterdict

#function to generate a dictionary with Newick_Species associated with SpeciesUID

def SpeciesUIDPairing(filename):

	UID_Speciesdict = {}

	with open(filename, 'r') as f : 
   
        	reader = csv.reader(f, delimiter = ',')
        	for row in reader :
           
                	#matching data with export of data from Navicat
                	SpeciesUID = row[0]
                	Species_name = row[1]
                	Newick_Formatted_Species = row[2]
                	Species_commom_name = row[3]
                	Highest_gold_status = row[4]
                	Study_ID = row[5]
                	Gene_count = row[6]
                	Ensemble_accession = row[7]
                	Ensembl_DB = row[8]
                	Gene_build_method = row[9]
                	TaxonID = row[10]
                	inBiomart = row[11]

			UID_Speciesdict[Newick_Formatted_Species] = SpeciesUID
               
	#print(UID_Speciesdict)

	return UID_Speciesdict



#function to tie everything together
def True_Blast_Pairing(prayer):
 

	cleansisterdict = Sisterspecies('SpeciesTable.txt')
	UID_Speciesdict = SpeciesUIDPairing('SpeciesTable.txt')

	cleanvalues = cleansisterdict.values()
	specified = UID_Speciesdict.keys()

	True_Blast_dict = cleansisterdict

	#print(True_Blast_dict)	

	#print(specified)
	#print(UID_Speciesdict)

	for key in cleansisterdict:
		
		if key in specified:
			c = key
			d = UID_Speciesdict[c]
			print(c)
			print(c,d)
			
			True_Blast_dict[d] = True_Blast_dict[c]
			del True_Blast_dict[c]
		else:
			print('nay', key)

	for value in cleansisterdict:
		
		if value in specified:
			e = value
			f = UID_Speciesdict[e]
			print(e)
			print(e,f)
			
			True_Blast_dict[d] = True_Blast_dict[c]
			del True_Blast_dict[c]
		else:
			print('nay', value)


	print(True_Blast_dict)

True_Blast_Pairing('prayer')
