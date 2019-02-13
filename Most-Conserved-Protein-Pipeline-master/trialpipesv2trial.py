# -*- coding: utf-8 -*-
"""
Created on Thursday Nov 7th 16:21:22 2018

@author: Catarina

THIS IS THE MAIN FILTER SCRIPT. There is a flowchart with a clearer picture of how this script runs, but the main breakdown is as follows:

1.UniqueGeneIDFilterFrequency - Double check that everything that was tagged as "notunique" is in face so, and generate a list of Genes in this category
2.Find length of Genelist (index check)
3.Run loop through the full index of Genes
	a.Blaster- blastp of the gene.txt file against teh sisterspecies blast database, built by running DatabaseBuilder.py
	b.BlastParser- takes the output of the Blaster (the hits in the sister genome) and formats for use in ...(the next step)
	c.Reciprocal Blaster- takes the reformatted BlastParser output file  and re-blastp-s the hits in the sister genome against thte original species' database
	d.Top Hit-  takes the returned hits from Reciprocal Blaster, organizes by bitscore, and picks the protein with the highest bitscore per gene. The top hit is then added to the top_hits_file.
	e.Garbage disposal- deletes the files generated in the pipeline process
	
NOTE: There are cases where the Blastp steps will yield no output within the evalue threshold, due to insufficient sequence conservation between the species and the sister species. 
These genes are listed as errros in Blaster and Reciprocal Blaster, generate no files to pasre, and are thus are shunted to the errors_speciesuid.txt file. 
The error files are then used in another script to pick proteins by maximum length for problematic genes.

THIS SCRIPT IS SLOW- depending on the length of the protein sequences, the similarity between species and sister species, and the number of actual hits, the Blaster stes can take 4-10 seconds apiece PER GENE.
If the first filter returns 5000 not uniquely paired genes in a species (for our dataset, this is the average per species), this script will take abot 10 hours to run.
We have 450 species. On one core, this thing takes A LONG TIME.
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

    	#print(notuniquegenedict)

	Genelist = notuniquegenedict.keys()
 	
    	return Genelist


#this function will make a list of all the gene names (and hence, file names) per species that will go through the pipeline

def filesforblaster_no_uniqueness_check(filename):

	Genelist = []

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

			if GeneID not in Genelist:
				Genelist.append(GeneID)
			
			else:
				pass
	return Genelist


# This function will take an input file (file per gene) and blast against a given database, and yield an output blastp file

def Blaster(file,databasefile,database,outfile):

	#os.system('/home/catherineweibel/Projects/pipeline/blastbinprograms/blastp -query %s -db ../%s/%sdb -evalue 1e-3 -outfmt "6 qseqid qacc qseq sseqid sacc sseq evalue bitscore" -max_target_seqs 2 |sort -u>%s'%(file, database,database,outfile))

	os.system('/home/catherineweibel/Projects/pipeline/blastbinprograms/blastp -query %s -db ../%s/%sdb -evalue 1e-3 -outfmt "6 qseqid qacc qseq sseqid sacc sseq evalue bitscore"  -num_alignments 3 -out %s'%(file, databasefile,database,outfile))
	
	print('blaster success for %s'%file)

def BlastParser(filename,outfile,database):

# This function creates an alignment of the ouput from the BlastParser (should align sequences with the top 3 hits from BLASTp)  - note this takes the ALIGNED sequences from Blastp
	exists = os.path.isfile(filename)
	if exists:
		print('file %s exists'%filename)
	

		with open(filename, 'r') as f :  
    
        		reader = csv.reader(f, delimiter = '	')
       			for row in reader :
				
				query_sequence_id = row[0]
				query_accession = row[1]
				query_sequence_aligned = row[2]
				subject_sequence_id = row[3]
				subject_accession = row[4]
				subject_sequence_aligned = row[5]
				e_value = row[6]
				bitscore = row[7]
			
				outputq = [query_sequence_id]
				seq = [query_sequence_aligned]
				outputs = [subject_sequence_id]
				ses = [subject_sequence_aligned]
				
				#get sequences into a string
				ses_string = ''.join(ses[0:])
			
				#remove hyphens to avoid reciprocal blast errors
				ses_string_no_hyphen = ses_string.replace("-", "")
				
				#writing to a new fastaformat for new blastp to align
			
			
				with open(outfile, "a") as f1:
					#f.write('>')
					#f.write(','.join(outputq[0:]) + '\n')
					#f.write(''.join(seq[0:]) + '\n' + '\n')
					f1.write('>')
					f1.write(','.join(outputs[0:]) + '\n')
					f1.write( ses_string_no_hyphen + '\n' + '\n')

		print('parser success for %s'%filename)

	else:
		print('something went wrong for %s in BlastParser'%filename)
		with open("errors_%s.txt"%(database), "a") as f2:
			f2.write(' something went wrong for %s in BlastParser'%filename + '\n')
			f2.close()

	
	

def Reciprocal_Blast(file, databasefile, database, outfile):
	
#this function aligns the hits from the Blaster back to the original genome by taking the fasta files of the hits generated by the parser and re-blasting them

	 os.system('/home/catherineweibel/Projects/pipeline/blastbinprograms/blastp -query %s -db ../%s/%sdb -evalue 1e-3 -outfmt "6 qseqid qacc qseq sseqid sacc sseq evalue bitscore"  -num_alignments 1 -out %s'%(file, databasefile,database,outfile))

	#print('reciprocal blaster success for %s'%file )
	

def Top_Hit(filename,outfile,database):

	bitscore_list= []

	exists = os.path.isfile(filename)
	if exists:
		print('file %s exists'%filename)

	

		with open(filename, 'r') as f :  
    
        		reader = csv.reader(f, delimiter = '	')
       			for row in reader :
			
				query_sequence_id = row[0]
				query_accession = row[1]
				query_sequence_aligned = row[2]
				subject_sequence_id = row[3]
				subject_accession = row[4]
				subject_sequence_aligned = row[5]
				e_value = row[6]
				bitscore = row[7]
			
				bitscore_list.append(bitscore)
		f.close()
		c = max(bitscore_list)
		print(c)
	
		with open(filename, 'r') as f :  
    
        		reader = csv.reader(f, delimiter = '	')
       			for row in reader :
			
				query_sequence_id = row[0]
				query_accession = row[1]
				query_sequence_aligned = row[2]
				subject_sequence_id = row[3]
				subject_accession = row[4]
				subject_sequence_aligned = row[5]
				e_value = row[6]
				bitscore = row[7]
			
				if bitscore == c:
					output = [subject_sequence_id, bitscore]
					
					with open(outfile, "a") as f1:
						f1.write(','.join(output[0:]) + '\n')
				
			f1.close()	
		print('Top Hits success for %s'%filename)
	else:
		print('something went wrong for %s in Top Hits'%filename)
		with open("errors_%s.txt"%(database), "a") as f2:
			f2.write('something went wrong for %s in Top Hits'%filename + '\n')
			f2.close()


	
def clustalaligner(file,outfile):

	#--in is infile, --out is outfile, -v means verbose output

	os.system('/home/catherineweibel/Projects/pipeline/Databases/./clustalo --in %s --dealign --resno --out %s'%(file,outfile))

	#os.system('/home/catherineweibel/Projects/pipeline/Databases/./clustalo --in %s -v --out %s'%(file,outfile))
	
	print('aligner success for %s'%file)

def garbagedisposal(file):

	os.system('rm %s'%file)
	
	print('garbage disposal success for %s'%file)


def Pipeline(filename, database_one,database_two):

	
	#list of genes for species
	Genelist = UniqueGeneIDFilterFrequency(filename)
	print(Genelist)

	#find length of list so we have an idea of the number of iterations
	Genelistsize = len(Genelist)
	print(Genelistsize)
	

	#loop through the Genes- 
	for i in range (0, Genelistsize):
		
		print(Genelist[i])
		#NOTE THIS TAKES FOREVER- blast takes 4-5 seconds per gene. Multiply that by around 30,000 for runtime per species- so roughly 41 hrs per species for this step alone
		Blaster("%s.txt"%Genelist[i],database_one,database_one,"%s_blasted.txt"%Genelist[i])

		#takes barely any time- note that first file is infile, second is outfile
		BlastParser("%s_blasted.txt"%Genelist[i],"blast_to_reciprocal_blast_%s.txt"%Genelist[i],database_two)

		#do it again
		Reciprocal_Blast("blast_to_reciprocal_blast_%s.txt"%Genelist[i], database_two, database_two, "reciprocal_blast_to_parser_%s.txt"%Genelist[i])
		
		#Top Hit runner
		Top_Hit("reciprocal_blast_to_parser_%s.txt"%Genelist[i],"real_top_hit_%s.txt"%database_two,database_two)

		#reparse for clustal alignment
		#BlastParser("reciprocal_blast_to_parser_%s.txt"%Genelist[i],"parser_to_clustal_%s.txt"%Genelist[i],database_two)		

		#first file in input_sequence_file, second is alignment
		#clustalaligner(" parser_to_clustal_%s.txt"%Genelist[i], "clustal_alignment_%s.txt"%Genelist[i])
		
		#start deleting
		garbagedisposal("%s_blasted.txt"%Genelist[i])

		#keep deleting
		garbagedisposal("blast_to_reciprocal_blast_%s.txt"%Genelist[i])
		
		#and some more
		garbagedisposal(" reciprocal_blast_to_parser_%s.txt"%Genelist[i])

		#clean some more
		garbagedisposal("%s.txt"%Genelist[i])

#Unique_data_reproduction("Outputnotuniquetrial211.txt")
Pipeline("Outputnotuniquetrial87.txt",n,m)
#n = sisterpecies
#m- species whose genes were split up into individual files

#USE SISTERPECIESDICT SCRIPT TO PICK PAIRS OF SPECIES FOR THIS



 
