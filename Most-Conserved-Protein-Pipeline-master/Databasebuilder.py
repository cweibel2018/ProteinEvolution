# -*- coding: utf-8 -*-
"""
Created on Thursday Nov 8th 16:24:56 2018

@author: Catarina

Note that this script must be run with python 3, and Blastx must me installed to call makeblastdb

I had to break up the steps as individual scripts because the sytax for step 1 only works in python3, and the syntax for steps 2-4 are written for python 2


"""

import mysql.connector
from itertools import chain
import sys, json, os, csv, time, copy, math


          
        		     
def databasebuilder(database):
	
    cnx = mysql.connector.connect(user = 'username',password= 'password' ,host='DNS' , database= 'database')

    mycursor = cnx.cursor(buffered = True)


	
#do this in batches, and change the range so the script can be run in parallel on multiple cores 
    for i in range(n,m):#substitute n,m for species UID
	
        extractionStatement= 'SELECT UID,CodingSeqTableUID,SpeciesUID,NewickSpecies,GeneID,ProteinID,ProteinSequence FROM Table WHERE SpeciesUID = %s'%i

        mycursor.execute(extractionStatement)
        
        #this exports the data from MySQL
	
        for entry in mycursor:
	
            string = ','.join(map(str,entry))
			#print(string)
			
            with open("%s.txt"%i, "a") as f1:
                f1.write(string + '\n')

        f1.close()
        print(i)

        print('done with step 1 of %s'%i)
        
        #this step turns the exported data into a fasta file

        with open("%s.txt"%i, 'r') as f :  
    
            reader = csv.reader(f, delimiter = ',')
            for row in reader :
            
                #matching data with export of data from Navicat
                    UID = row[0]
                    CodingSeqTableUID = row[1]
                    SpeciesUID = row[2]
                    Newick_Species = row[3]
                    GeneID = row[4]
                    ProteinID = row[5]
                    ProteinSequence = row[6]
               			
                    output = [UID,SpeciesUID, Newick_Species,GeneID,ProteinID]
                    outputc = [ProteinSequence]
			
			#turning the .txt file into a fasta format PER SPECIES
                    with open("fasta%s.txt"%i,"a") as f2:
                        f2.write('>')
                        f2.write(','.join(output[0:]) + '\n')
                        f2.write(''.join(outputc[0:]) +'\n' + '\n')
            f2.close()
			
        print('done with step 2 of %s'%i)        

        #this step generates the blast databases and appropriate indexes. NOTE THIS MUST HAVE THE -hash_index FOR BLASTp AGAINST THIS DATABASE TO WORK
        os.system("path-to/blastbinprograms/makeblastdb -in fasta%s.txt -dbtype 'prot' -input_type fasta -out %sdb -hash_index"%(i,i))

        print('done with step 3 of %s'%i) 
        
        #this step makes the path to a species-specific directory
        os.system("path-to/Databases/mkdir %s"%i)

        print('done with step 4 of %s'%i)
        
    	#this step moves the files to the database
        os.system("path-to/Databases/mv %s"%i)
    
        print('done building databases')
        

	
#run the program				
databasebuilder('NCBIGenomes_Protein_Complete')
