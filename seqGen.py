#!/usr/bin/env python3.4
__author__ = 'mdc_hk'
version = '1.0'

import logging, os, sys
import pandas as pd
from pandas import Series
import re
logging.basicConfig(filename='Log.txt', level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

try:
    variant_xls=sys.argv[1]
    referencefa=sys.argv[2]
    proceed = True
    # print('\nProgram: SeqGen Version: {} is running'.format(version))
except:
    print('\nProgram: SeqGen\nVersion: {}'.format(version))
    print('Usage:', sys.argv[0], '<Variant_xls>', '<reference.fa>')
    print('Example:', sys.argv[0], 'sample001_AllGenes.xls>', 'reference_annotated.fa')
    proceed = False

if proceed == True:
    def iupacConverter2n(nucleotide1, nucleotide2):
        try:
            if nucleotide1 + nucleotide2 == 1:
                iupac = 'A'
            elif nucleotide1 + nucleotide2 == 2:
                iupac = 'T'
            elif nucleotide1 + nucleotide2 == 3:
                iupac = 'W'
            elif nucleotide1 + nucleotide2 == 4:
                iupac = 'C'
            elif nucleotide1 + nucleotide2 == 5:
                iupac = 'M'
            elif nucleotide1 + nucleotide2 == 6:
                iupac = 'Y'
            elif nucleotide1 + nucleotide2 == 8:
                iupac = 'G'
            elif nucleotide1 + nucleotide2 == 9:
                iupac = 'R'
            elif nucleotide1 + nucleotide2 == 10:
                iupac = 'K'
            elif nucleotide1 + nucleotide2 == 12:
                iupac = 'S'
        except:
            iupac = 'X'
            logging.info('Deletion/Insertion found in ' + segment[:8] + ' sample ' + stripped_root[0][1][:-len(suffix)])
            errorFile = open('Error.txt', 'a')
            errorFile.write('Deletion/Insertion found in ' + segment[:8] + ' sample ' + stripped_root[0][1][:-len(suffix)] + '\n')
            errorFile.close()
        return iupac

    def iupacConverter3n(nucleotide1, nucleotide2, nucleotide3):
        try:
            if nucleotide1 + nucleotide2 + nucleotide3 == 1:
                iupac = 'A'
            elif nucleotide1 + nucleotide2 + nucleotide3 == 2:
                iupac = 'T'
            elif nucleotide1 + nucleotide2 + nucleotide3 == 3:
                iupac = 'W'
            elif nucleotide1 + nucleotide2 + nucleotide3 == 4:
                iupac = 'C'
            elif nucleotide1 + nucleotide2 + nucleotide3 == 5:
                iupac = 'M'
            elif nucleotide1 + nucleotide2 + nucleotide3 == 6:
                iupac = 'Y'
            elif nucleotide1 + nucleotide2 + nucleotide3 == 8:
                iupac = 'G'
            elif nucleotide1 + nucleotide2 + nucleotide3 == 9:
                iupac = 'R'
            elif nucleotide1 + nucleotide2 + nucleotide3 == 10:
                iupac = 'K'
            elif nucleotide1 + nucleotide2 + nucleotide3 == 12:
                iupac = 'S'
            elif nucleotide1 + nucleotide2 + nucleotide3 == 7:
                iupac = 'H'
            elif nucleotide1 + nucleotide2 + nucleotide3 == 11:
                iupac = 'D'
            elif nucleotide1 + nucleotide2 + nucleotide3 == 13:
                iupac = 'V'
            elif nucleotide1 + nucleotide2 + nucleotide3 == 14:
                iupac = 'B'
        except TypeError:
            iupac = 'X'
            logging.info('Deletion/Insertion found in ' + segment[:8] + ' sample ' + stripped_root[0][1][:-len(suffix)])
            errorFile = open('Error.txt', 'a')
            errorFile.write('Deletion/Insertion found in ' + segment[:8] + ' sample ' + stripped_root[0][1][:-len(suffix)] + '\n')
            errorFile.close()
        return iupac


    with open(referencefa, mode="r") as file1:
        content1 = file1.read().splitlines()

    geneName = []
    for i in range(len(content1)):
        if content1.index(content1[i]) % 2 == 0 or content1.index(content1[i]) == 0:
            geneName.append(content1[i][1:])
    geneCont = []
    for i in range(len(content1)):
        if content1.index(content1[i]) % 2 != 0:
            geneCont.append(content1[i])

    suffix = '_AllGenes'
    root, extension = os.path.splitext(variant_xls)
    regEx = re.compile(r'([\S]*/([^/]*).xls)')
    stripped_root = regEx.findall(variant_xls)

    dictRefGene = dict(zip(geneName, geneCont))
    df = pd.read_excel(variant_xls)
    df2 = df[1:]
    df2.loc[:, 'position'] = df2.index

    df2 = df2.set_index(['CHROM'])
    cummulativeseq = []

    def seqGenerator(gene, outputfile):
        global segment
        segment = gene
        try:
            dfGene = df2.ix[gene]
            proceed = True
        except:
            proceed = False
            print('No variant found for ' + gene + ' gene segment for sample: ' + stripped_root[0][1][:-len(suffix)])

        if proceed == True:

            dfGene = dfGene.reset_index()

            keyDelBack = dfGene.index.max()
            while dfGene.ix[keyDelBack, 'Call(R1,R2)'] == 'Invalid':
                keyDelBack -= 1
            if keyDelBack != dfGene.index.max():
                dfGene = dfGene.drop(dfGene.index[keyDelBack+1:])

            keyDelFront = 0
            while dfGene.ix[keyDelFront, 'Call(R1,R2)'] == 'Invalid':
                keyDelFront += 1

            if keyDelFront != 0:
                dfGene = dfGene.drop(dfGene.index[:keyDelFront])


            dfGene = dfGene.set_index(['position'])
            geneSeq = list(dictRefGene[gene])
            takeSequence = True
            indexList = dfGene.index.tolist()

            for pos in indexList:
                nucleotide12 = dfGene['Call(R1,R2)'].str.split(",").apply(Series, 1)
                if len(nucleotide12.columns) == 1:
                    if not nucleotide12.ix[pos, 0] == 'Invalid':
                        geneSeq[pos-1] = nucleotide12.ix[pos, 0]
                    else:
                        takeSequence = False
                        break
                        # geneSeq[pos-1] = 'X'
                elif len(nucleotide12.columns) == 2:
                    try:
                        if nucleotide12.ix[pos, 0] == 'Invalid':
                            takeSequence = False
                            break
                        else:
                            nucleotide12[1].fillna(0, inplace=True)
                            nucleotide12 = nucleotide12.replace(['A','T','C','G'], [1, 2, 4, 8])
                            geneSeq[pos-1] = iupacConverter2n(nucleotide12.ix[pos, 0], nucleotide12.ix[pos, 1])
                    except ValueError:
                            geneSeq[pos-1] = 'X'
                            logging.info('Deletion/Insertion found in ' + segment[:8] + ' sample ' + stripped_root[0][1][:-len(suffix)])
                            errorFile = open('Error.txt', 'a')
                            errorFile.write('Deletion/Insertion found in ' + segment[:8] + ' sample ' + stripped_root[0][1][:-len(suffix)] + '\n')
                            errorFile.close()

                elif len(nucleotide12.columns) == 3:
                    try:
                        if nucleotide12.ix[pos, 0] == 'Invalid':
                            takeSequence = False
                            break
                        else:
                            nucleotide12[1].fillna(0, inplace=True)
                            nucleotide12[2].fillna(0, inplace=True)
                            nucleotide12 = nucleotide12.replace(['A','T','C','G'], [1, 2, 4, 8])
                            geneSeq[pos-1] = iupacConverter3n(nucleotide12.ix[pos, 0], nucleotide12.ix[pos, 1],
                                                        nucleotide12.ix[pos, 2])
                    except ValueError:
                            geneSeq[pos-1] = 'X'
                            logging.info('Deletion/Insertion found in ' + segment[:8] + ' sample ' + stripped_root[0][1][:-len(suffix)])
                            errorFile = open('Error.txt', 'a')
                            errorFile.write('Deletion/Insertion found in ' + segment[:8] + ' sample ' + stripped_root[0][1][:-len(suffix)] + '\n')
                            errorFile.close()
                else:
                    geneSeq[pos-1] = 'N'
            if takeSequence == True:
                with open(outputfile, mode='a') as file1:
                    file1.write('>' + stripped_root[0][1][:-len(suffix)] + '_' + gene + '\n' + ''.join(geneSeq) + '\n')
            else:
                logging.info(gene[:8] + ' gene of ' + stripped_root[0][1][:-len(suffix)] + ' was not generated due to low read counts.')
                errorFile = open('Error.txt', mode='a')
                errorFile.write(gene[:8] + ' gene of ' + stripped_root[0][1][:-len(suffix)] + ' was not generated due to low read counts.'+'\n')
                errorFile.close()

    with open(referencefa, mode="r") as file1:
        listingName = file1.read()

    segmentRegex = re.compile(r'((Segment\d)_[\S]*)')
    segmentTolist = segmentRegex.findall(listingName)

    for seg in segmentTolist:
        seqGenerator(seg[0], 'MiseqSeq_' + seg[1] +'.txt')

    print('Genome sequences generated for sample ' + stripped_root[0][1][:-len(suffix)])