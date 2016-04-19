#!/usr/bin/env python3.4
# Created by: Lee Hong Kai
# Descript: Converts multiline FASTAs to single line FASTAs
#  
# Usage: DBupdate.py <sequences.fa>
# Example: DBupdate.py mySeqs.fa
#----------------------------------------------------------------------------------------
#===========================================================================================================
#Imports:
import os, re, sys, transaction
from ZODB import FileStorage, DB
#===========================================================================================================
# Functions:

# 1: Checks if in proper number of arguments are passed gives instructions on proper use.
def argsCheck(numArgs):
	if len(sys.argv) < numArgs or len(sys.argv) > numArgs:
		print("Converts multiline FASTAs to single line FASTAs")
		print("Usage: " + sys.argv[0] + " <sequences.fa>")
		print("Examples: " + sys.argv[0] + " mySeqs.fa")
		exit(1) # Aborts program. (exit(1) indicates that an error occurred)
#===========================================================================================================
# Main program code:
	
# House keeping...
argsCheck(2) # Checks if the number of arguments are correct.

# Stores file one for input checking.
inFile  = sys.argv[1]

print(">> Opening FASTA file...")
# Reads sequence file list and stores it as a string object. Safely closes file:
try:	
	with open(inFile,"r") as newFile:
		sequences = newFile.read()
		sequences = re.split("^>", sequences, flags=re.MULTILINE) # Only splits string at the start of a line.
		del sequences[0] # The first fasta in the file is split into an empty empty element and and the first fasta
						 # Del removes this empty element.
		newFile.close()
except IOError:
	print("Failed to open " + inFile)
	exit(1)

print(">> Converting FASTA file from multiline to single line and writing to file.")
# Conversts multiline fasta to single line. Writes new fasta to file.

try:
	os.unlink('FASTA.dat')
	os.unlink('FASTA.dat.index')
	os.unlink('FASTA.dat.lock')
	os.unlink('FASTA.dat.tmp')
except FileNotFoundError:
	None

try:
	storage = FileStorage.FileStorage('FASTA.dat')
	db = DB(storage)
	connection = db.open()
	root = connection.root()
	for fasta in sequences:
		try:
			header, sequence = fasta.split("\n", 1) # Split each fasta into header and sequence.
			sequence = sequence.replace("\n","")# Replace newlines in sequence
		except ValueError:
			print(fasta)
		root[header] = sequence
	transaction.commit()
	storage.close()
except IOError:
	print("Failed to open " + inFile)
	exit(1)

print(">> Done!")
