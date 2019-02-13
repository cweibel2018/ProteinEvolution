# -*- coding: utf-8 -*-
"""
Created on Tuesday October 9 16:24:56 2018

@author: Catarina

This script extracts the contents of MySQL tables for use in the main pipeline BY SPECIES UID.

Note: mysql.connector and python 3 dependencies
"""
import mysql.connector
from itertools import chain




def speciesextractor(database):
    
    #put in MySQL info below
	
	cnx = mysql.connector.connect(user = 'user',password= 'password' ,host='DNS' , database= 'database')

	mycursor = cnx.cursor(buffered = True)

#extractor
#*=entire row
	
#do this in batches, and change the range so the script can be run in parallel on multiple cores 
	for i in range(n,m):#n,m are the species UIDs
		print(i)
		extractionStatement= 'SELECT UID,SpeciesUID,NewickSpecies,GeneID,ProteinID,CodingSequence FROM table WHERE SpeciesUID = %s'%i

		mycursor.execute(extractionStatement)
        
        #writing the file
	
		for entry in mycursor:
	
			string = ','.join(map(str,entry))
			#print(string)
			
			with open("%s.txt"%i, "a") as f1:
				f1.write(string + '\n')

			f1.close()
		print(i)

		print('done with step 1 of %s'%i)



speciesextractor('NCBI')


