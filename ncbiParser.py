#!/usr/bin/env python3.4
__author__ = 'mdc_hk'
# Description: To retrive the fasta sequences from the database, according to the blastn results.out
# Usage: NCBIparser.py <blastn.out> <Sample ID> <Output Directory>
# Example: NCBIparser.py <blastn.out> <sample001> <PATH>

import time, datetime, logging, os, sys
import pandas as pd
from pandas import Series
import numpy as np
from ZODB import FileStorage, DB
logging.basicConfig(filename='Log.txt', level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

#===========================================================================================================
# Functions:

# 1: Checks if in proper number of arguments are passed gives instructions on proper use.
def argsCheck(numArgs):
	if len(sys.argv) < numArgs or len(sys.argv) > numArgs:
		print('To retrive the fasta sequences from the database, according to the blastn results output')
		print('Usage: ' + sys.argv[0] + ' <blastn.out>' + ' <Sample ID>' + ' <Output Directory>')
		print('Examples: ' + sys.argv[0] + ' blastn.out' + ' sample001', '<Directory>')
		exit(1) # Aborts program. (exit(1) indicates that an error occurred)
#===========================================================================================================
# Main program code:
startTime = time.time()

# Housekeeping.
argsCheck(4) # Checks if the number of arguments are correct.

# Stores file one for input checking.
inFile  = sys.argv[1]
sampleID = sys.argv[2]
outputDirectory = sys.argv[3]


programFolder = os.path.expanduser('~/FluSeq/')
startTime = time.time()
print('>> '+ datetime.datetime.now().strftime('%Y/%m/%d %H:%M:%S'), '- Parsing the blastn results for fasta sequences '
                                                                    'retriving from the database...')
try:
    df = pd.read_table(inFile, sep='\t', names=['query id', 'subject id', '% identity', 'alignment length',
                                                   'mismatches', 'gap opens', 'q. start', 'q. end', 's. start',
                                                   's. end', 'evalue', 'bit score'])
except OSError:
    logging.info('No blastn output recorded for sample ' + sampleID + '; Possible solution: Regenerate the pair-end '
                'fastq files using bcl2fastq2 program')
    errorFile = open('Error.txt', 'a')
    errorFile.write('No blastn output recorded for sample ' + sampleID + '; Possible solution: Regenerate the '
                        'pair-end fastq files using bcl2fastq2 program' + '\n')
    errorFile.close()
    exit(1)

df_tojoin = df['subject id'].str.split("|").apply(Series, 1)
if len(df.columns) == 2:
    df_tojoin[2] = np.NaN
    df_tojoin[3] = np.NaN
elif len(df.columns) == 3:
    df_tojoin[3] = np.NaN
df_tojoin.columns = ['GenBank Accession', 'Strain ID', 'Segment', 'Serotype']

try:
    df = pd.merge(df, df_tojoin, left_index=True, right_index=True)
except IndexError:
    logging.info('The subject ids for the blastn output are not recorded in desired format for sample ' + sampleID + 
                 '; Possible solution: Reformat it the subject id as such: '
                 '\'AJ123456|StrainID|SegmentNumber|Serotype\'')
    errorFile = open('Error.txt', 'a')
    errorFile.write('The subject ids for the blastn output are not recorded in desired format for sample ' + sampleID 
                       + '; Possible solution: Reformat it the subject id as such: '
                         '\'AJ123456|StrainID|SegmentNumber|Serotype\'' + '\n')
    errorFile.close()
    exit(1)

df_final=df.iloc[df.groupby('Segment')['alignment length'].agg(pd.Series.idxmax)]

for index, row in df_final.iterrows():
    if row['Segment'] == str(1):
        Segment1 = df.ix[index]['subject id']
    elif row['Segment'] == str(2):
        Segment2 = df.ix[index]['subject id']
    elif row['Segment'] == str(3):
        Segment3 = df.ix[index]['subject id']
    elif row['Segment'] == str(4):
        Segment4 = df.ix[index]['subject id']
    elif row['Segment'] == str(5):
        Segment5 = df.ix[index]['subject id']
    elif row['Segment'] == str(6):
        Segment6 = df.ix[index]['subject id']
    elif row['Segment'] == str(7):
        Segment7 = df.ix[index]['subject id']
    elif row['Segment'] == str(8):
        Segment8 = df.ix[index]['subject id']

print('>> '+ datetime.datetime.now().strftime('%Y/%m/%d %H:%M:%S'), '- Retrieving fasta sequences from the database...')

# To retrieve fasta sequences from database.
try:
    storage = FileStorage.FileStorage(programFolder + 'FASTADatabase/FASTA.dat')
    db = DB(storage)
    connection = db.open()
    root = connection.root()
except:
    logging.info('FASTA.dat database file was not found')
    errorFile = open('Error.txt', 'a')
    errorFile.write('FASTA.dat database file was not found. Suggested solution: Please make sure that FASTA.dat '
                       'is placed in the correct PATH' + '\n')
    errorFile.close()
    exit(1)

with open(outputDirectory + '/' + sampleID + '_Reference.fa', 'w') as fastaFile:
    try:
        segmentFasta1 = '>Segment1_' + Segment1 + '\n' + root[Segment1] + '\n'
        fastaFile.write(segmentFasta1)
    except NameError:
        print('No gene segment 1 found in blastn.out.')
        logging.info('No gene segment 1 found in blastn.out for sample ' + sampleID + '; Possible solution: '
                        'Use a pre-determined reference fasta for analyses')
        errorFile = open('Error.txt', 'a')
        errorFile.write('No gene segment 1 found in blastn.out for sample ' + sampleID + '; Possible solution: '
                        'Use a pre-determined reference fasta for analyses' + '\n')
        errorFile.close()

    try:
        segmentFasta2 = '>Segment2_' + Segment2 + '\n' + root[Segment2] + '\n'
        fastaFile.write(segmentFasta2)
    except NameError:
        print('No gene segment 2 found in blastn.out.')
        logging.info('No gene segment 2 found in blastn.out for sample ' + sampleID + '; Possible solution: '
                        'Use a pre-determined reference fasta for analyses')
        errorFile = open('Error.txt', 'a')
        errorFile.write('No gene segment 2 found in blastn.out for sample ' + sampleID + '; Possible solution: '
                        'Use a pre-determined reference fasta for analyses' + '\n')
        errorFile.close()

    try:
        segmentFasta3 = '>Segment3_' + Segment3 + '\n' + root[Segment3] + '\n'
        fastaFile.write(segmentFasta3)
    except NameError:
        print('No gene segment 3 found in blastn.out.')
        logging.info('No gene segment 3 found in blastn.out for sample ' + sampleID + '; Possible solution: '
                        'Use a pre-determined reference fasta for analyses')
        errorFile = open('Error.txt', 'a')
        errorFile.write('No gene segment 3 found in blastn.out for sample ' + sampleID + '; Possible solution: '
                        'Use a pre-determined reference fasta for analyses' + '\n')
        errorFile.close()

    try:
        segmentFasta4 = '>Segment4_' + Segment4 + '\n' + root[Segment4] + '\n'
        fastaFile.write(segmentFasta4)
    except NameError:
        print('No gene segment 4 found in blastn.out.')
        logging.info('No gene segment 4 found in blastn.out for sample ' + sampleID + '; Possible solution: '
                        'Use a pre-determined reference fasta for analyses')
        errorFile = open('Error.txt', 'a')
        errorFile.write('No gene segment 4 found in blastn.out for sample ' + sampleID + '; Possible solution: '
                        'Use a pre-determined reference fasta for analyses' + '\n')
        errorFile.close()


    try:
        segmentFasta5 = '>Segment5_' + Segment5 + '\n' + root[Segment5] + '\n'
        fastaFile.write(segmentFasta5)
    except NameError:
        print('No gene segment 5 found in blastn.out.')
        logging.info('No gene segment 5 found in blastn.out for sample ' + sampleID + '; Possible solution: '
                        'Use a pre-determined reference fasta for analyses')
        errorFile = open('Error.txt', 'a')
        errorFile.write('No gene segment 5 found in blastn.out for sample ' + sampleID + '; Possible solution: '
                        'Use a pre-determined reference fasta for analyses' + '\n')
        errorFile.close()


    try:
        segmentFasta6 = '>Segment6_' + Segment6 + '\n' + root[Segment6] + '\n'
        fastaFile.write(segmentFasta6)
    except NameError:
        print('No gene segment 6 found in blastn.out.')
        logging.info('No gene segment 6 found in blastn.out for sample ' + sampleID + '; Possible solution: '
                        'Use a pre-determined reference fasta for analyses')
        errorFile = open('Error.txt', 'a')
        errorFile.write('No gene segment 6 found in blastn.out for sample ' + sampleID + '; Possible solution: '
                        'Use a pre-determined reference fasta for analyses' + '\n')
        errorFile.close()


    try:
        segmentFasta7 = '>Segment7_' + Segment7 + '\n' + root[Segment7] + '\n'
        fastaFile.write(segmentFasta7)
    except NameError:
        print('No gene segment 7 found in blastn.out.')
        logging.info('No gene segment 7 found in blastn.out for sample ' + sampleID + '; Possible solution: '
                        'Use a pre-determined reference fasta for analyses')
        errorFile = open('Error.txt', 'a')
        errorFile.write('No gene segment 7 found in blastn.out for sample ' + sampleID + '; Possible solution: '
                        'Use a pre-determined reference fasta for analyses' + '\n')
        errorFile.close()


    try:
        segmentFasta8 = '>Segment8_' + Segment8 + '\n' + root[Segment8] + '\n'
        fastaFile.write(segmentFasta8)
    except NameError:
        print('No gene segment 8 found in blastn.out.')
        logging.info('No gene segment 8 found in blastn.out for sample ' + sampleID + '; Possible solution: '
                        'Use a pre-determined reference fasta for analyses')
        errorFile = open('Error.txt', 'a')
        errorFile.write('No gene segment 8 found in blastn.out for sample ' + sampleID + '; Possible solution: '
                        'Use a pre-determined reference fasta for analyses' + '\n')
        errorFile.close()

storage.close()


with open(outputDirectory + '/' + sampleID + '_Reference.fa','r') as fastaFile:
    content1 = fastaFile.read().splitlines()

geneName = []
for i in range(len(content1)):
    if content1.index(content1[i]) % 2 == 0 or content1.index(content1[i]) == 0:
        geneName.append(content1[i][1:])

geneCont = []
for i in range(len(content1)):
    if content1.index(content1[i]) % 2 != 0:

        # Replacing all degenerated sequences to ATCG...
        content1[i] = content1[i].replace('R', 'A')
        content1[i] = content1[i].replace('Y', 'C')
        content1[i] = content1[i].replace('S', 'G')
        content1[i] = content1[i].replace('W', 'A')
        content1[i] = content1[i].replace('K', 'G')
        content1[i] = content1[i].replace('M', 'A')
        content1[i] = content1[i].replace('B', 'C')
        content1[i] = content1[i].replace('D', 'A')
        content1[i] = content1[i].replace('H', 'A')
        content1[i] = content1[i].replace('V', 'A')
        content1[i] = content1[i].replace('N', 'A')

        # Annotating sequences...
        sequence = list(content1[i])
        for k in range(0, len(sequence), 150):
            if sequence[k] == 'A':
                sequence[k] = 'T'
            else:
                sequence[k] = 'A'
        joinsequence = ''.join(sequence)
        geneCont.append(joinsequence)

with open(outputDirectory + '/' + sampleID + '_Reference_annotated.fa', mode='w') as file1:
    for i in range(len(geneName)):
        file1.write('>' + geneName[i] + '\n' + geneCont[i] + '\n')

endTime = time.time()
print('Took {} seconds to complete.'.format(round(endTime-startTime, 2)))
