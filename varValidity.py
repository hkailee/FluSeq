#!/usr/bin/env python3.4
__author__ = 'mdc_hk'
version = '1.0'

import logging, os, re, sys
import pandas as pd
from pandas import Series
from bs4 import BeautifulSoup
import numpy as np
from scipy import stats
logging.basicConfig(filename='Log.txt', level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

#===========================================================================================================

try:
    filein=sys.argv[1]
    xmlResequencingRunStatistics=sys.argv[2]
    reference_annotated_fa = sys.argv[3]
    confidenceValue = sys.argv[4]
    proceed = True
except:
    print('\nProgram: varValidity\nVersion: {}\n'.format(version) + 'Usage: ' + sys.argv[0] + ' <input_txt_file> '
            '<ResequencingRunStatistics_xml> <reference_annotated_fa> <Confidence value>')
    print('Examples:', sys.argv[0], 'sample001_table.txt', 'ResequencingRunStatistics.xml', 'reference_annotated.fa',
          '0.9500')
    proceed = False

if proceed == True:

    programFolder = os.path.expanduser('~/FluSeq/')
    dataFromTable = pd.read_table(filein, sep='\t')
    columnID = dataFromTable.columns[7]

    try:
        # dfStat = pd.read_csv('ms2SysCtrl.db', sep=',')
        dfStat = pd.read_csv(programFolder + 'ms2SysCtrl.db', sep=',')
        dfStat = dfStat.set_index(['GeneSegment', 'POS'])

    except IOError:
        errorFile = open('Error.txt', 'a')
        errorFile.write('Please make sure the MS2SysCtrl.db is placed in the correct PATH' + '\n')
        errorFile.close()
        exit(1)

    # The statistics for contamination starts here.

    def expectedConR(x2, gene, pos):
        dfStat_tmp = dfStat.ix[gene, pos]
        x = dfStat_tmp['PFreads_inMillion']
        y = dfStat_tmp['Total_Depth']

        # Modeling with Numpy
        p, cov = np.polyfit(x,y,1,cov=True)           # parameters and covariance from of the fit
        y_model = np.polyval(p, x)                    # model using the fit parameters; NOTE: parameters here are coefficients

        # Statistics
        n = y.size                                    # number of observations
        m = p.size                                    # number of parameters
        DF = n - m                                    # degrees of freedom
        t = stats.t.ppf(float(confidenceValue), n - m)       # used for CI and PI bands

        # Estimates of Error in Data/Model
        resid = y - y_model
        s_err = np.sqrt(np.sum(resid**2)/(DF))        # standard deviation of the error

        y2 = np.polyval(p, x2)

        # Prediction Interval
        PI = t*s_err*np.sqrt(1+1/n+(x2-np.mean(x))**2/np.sum((x-np.mean(x))**2))

        uppConfLimit_ContRead = y2 + PI

        return uppConfLimit_ContRead

### To import the main table result
    dataFromTable['SAMPLE'] = columnID[:-3]
    dataFromTable = dataFromTable.set_index(['SAMPLE'])

### to create table/library to get the percentage of reads identified for each sample

    try:
        bsObj = BeautifulSoup(open(xmlResequencingRunStatistics), 'lxml-xml')
        xmlOfSamples = bsObj.findAll(re.compile(r'SummarizedSampleStatis(t)*ics'))

    except IOError:
        errorFile = open('Error.txt', 'a')
        errorFile.write('Please make sure the run-specific ResequencingRunStatistics.xml is placed in the run folder'
                        + '\n')
        errorFile.close()
        exit(1)

    listOfSamples = []
    for name in xmlOfSamples:
        listOfSamples.append(name.SampleName.get_text())

    listOfPFReads = []
    for name in xmlOfSamples:
        listOfPFReads.append(round(int(name.NumberOfClustersPF.get_text())/1000000, 2))
    if len(listOfSamples) == len(listOfPFReads):
        dictOfSamplesAndPFreads = dict(zip(listOfSamples, listOfPFReads))

    dataFromTable['PFreads_inMillion'] = dictOfSamplesAndPFreads[columnID[:-3]]

### to set the index as 'CHROM' for easy gene data separation
    dataFromTable = dataFromTable.set_index(['CHROM'], drop=False)
    geneList = []

    def validityTester(gene):
        validityOfWildtype = []
        for index, row in gene.iterrows():
            if row['refCopy'] > row['UppPredLimit_ContRead']:
                validityOfWildtype.append(1)
            else:
                validityOfWildtype.append(0)

        validityOfAlt1 = []
        for index, row in gene.iterrows():
            if row['alt1Copy'] > row['UppPredLimit_ContRead']:
                validityOfAlt1.append(2)
            else:
                validityOfAlt1.append(0)

        validityOfAlt2 = []
        for index, row in gene.iterrows():
            if row['alt2Copy'] > row['UppPredLimit_ContRead']:
                validityOfAlt2.append(4)
            else:
                validityOfAlt2.append(0)

        scoring = [x + y + z for x, y, z in zip(validityOfWildtype, validityOfAlt1, validityOfAlt2)]
        Ref = list(gene['REF'])
        Alt12 = gene['ALT'].str.split(",").apply(Series, 1)
        if len(Alt12.columns) == 1:
            Alt12[1] = ''
        Alt1 = list(Alt12[0])
        Alt2 = list(Alt12[1])
        refCopyNumber = list(gene['refCopy'])
        alt1CopyNumber = list(gene['alt1Copy'])
        alt2CopyNumber = list(gene['alt2Copy'])

        call = []
        percentage = []
        for i in range(len(scoring)):
            if scoring[i] == 0:
                call.append('Invalid')
                percentage.append(np.NaN)
            elif scoring[i] == 1:
                call.append(Ref[i])
                percentage.append(1.00)
            elif scoring[i] == 2:
                call.append(Alt1[i])
                percentage.append(1.00)
            elif scoring[i] == 3:
                call.append(Ref[i]+','+Alt1[i])
                percentage.append(alt1CopyNumber[i]/(alt1CopyNumber[i]+refCopyNumber[i]))
            elif scoring[i] == 4:
                call.append(Alt2[i])
                percentage.append(1.00)
            elif scoring[i] == 5:
                call.append(Ref[i]+','+Alt2[i])
                percentage.append(alt2CopyNumber/(alt2CopyNumber[i]+refCopyNumber[i]))
            elif scoring[i] == 6:
                call.append(Alt1[i]+','+Alt2[i])
                percentage.append(alt2CopyNumber[i]/(alt1CopyNumber[i]+alt2CopyNumber[i]))
            elif scoring[i] == 7:
                call.append(Ref[i]+','+Alt1[i]+','+Alt2[i])
                percentage.append(alt2CopyNumber[i]/(alt1CopyNumber[i]+alt2CopyNumber[i]))

        gene['Call(R1,R2)'] = call
        gene['Percentage(R2/R1+R2)'] = [ '%.3f' % elem for elem in percentage ]

        return

    def geneAnalysis(Chrom, Gene):
        GeneData = dataFromTable.ix[Chrom]
        GeneData = GeneData.reset_index(drop=True)
        UppPredLimit_ContRead = []
        for i in GeneData.index:
            if GeneData.ix[i]['POS'] <= 225:
                UppPredLimit_ContRead.append(expectedConR(GeneData.ix[i]['PFreads_inMillion'], Gene, 151))
            elif 226 <= GeneData.ix[i]['POS'] <= 450:
                UppPredLimit_ContRead.append(expectedConR(GeneData.ix[i]['PFreads_inMillion'], Gene, 301))
            elif 451 <= GeneData.ix[i]['POS'] <= 525:
                UppPredLimit_ContRead.append(expectedConR(GeneData.ix[i]['PFreads_inMillion'], Gene, 451))
            elif 526 <= GeneData.ix[i]['POS'] <= 675:
                UppPredLimit_ContRead.append(expectedConR(GeneData.ix[i]['PFreads_inMillion'], Gene, 601))
            elif 676 <= GeneData.ix[i]['POS'] <= 825:
                UppPredLimit_ContRead.append(expectedConR(GeneData.ix[i]['PFreads_inMillion'], Gene, 751))
            elif 826 <= GeneData.ix[i]['POS'] <= 975:
                try:
                    UppPredLimit_ContRead.append(expectedConR(GeneData.ix[i]['PFreads_inMillion'], Gene, 901))
                except:
                    UppPredLimit_ContRead.append(expectedConR(GeneData.ix[i]['PFreads_inMillion'], Gene, 751))
            elif 976 <= GeneData.ix[i]['POS'] <= 1125:
                try:
                    UppPredLimit_ContRead.append(expectedConR(GeneData.ix[i]['PFreads_inMillion'], Gene, 1051))
                except:
                    UppPredLimit_ContRead.append(expectedConR(GeneData.ix[i]['PFreads_inMillion'], Gene, 901))
            elif 1126 <= GeneData.ix[i]['POS'] <= 1275:
                try:
                    UppPredLimit_ContRead.append(expectedConR(GeneData.ix[i]['PFreads_inMillion'], Gene, 1201))
                except:
                    UppPredLimit_ContRead.append(expectedConR(GeneData.ix[i]['PFreads_inMillion'], Gene, 1051))
            elif 1276 <= GeneData.ix[i]['POS'] <= 1425:
                try:
                    UppPredLimit_ContRead.append(expectedConR(GeneData.ix[i]['PFreads_inMillion'], Gene, 1351))
                except:
                    UppPredLimit_ContRead.append(expectedConR(GeneData.ix[i]['PFreads_inMillion'], Gene, 1201))
            elif 1426 <= GeneData.ix[i]['POS'] <= 1575:
                try:
                    UppPredLimit_ContRead.append(expectedConR(GeneData.ix[i]['PFreads_inMillion'], Gene, 1501))
                except:
                    UppPredLimit_ContRead.append(expectedConR(GeneData.ix[i]['PFreads_inMillion'], Gene, 1351))
            elif 1576 <= GeneData.ix[i]['POS'] <= 1725:
                try:
                    UppPredLimit_ContRead.append(expectedConR(GeneData.ix[i]['PFreads_inMillion'], Gene, 1651))
                except:
                    UppPredLimit_ContRead.append(expectedConR(GeneData.ix[i]['PFreads_inMillion'], Gene, 1501))
            elif 1726 <= GeneData.ix[i]['POS'] <= 1801:
                try:
                    UppPredLimit_ContRead.append(expectedConR(GeneData.ix[i]['PFreads_inMillion'], Gene, 1801))
                except:
                    UppPredLimit_ContRead.append(expectedConR(GeneData.ix[i]['PFreads_inMillion'], Gene, 1651))
            elif 1802 <= GeneData.ix[i]['POS'] <= 1951:
                try:
                    UppPredLimit_ContRead.append(expectedConR(GeneData.ix[i]['PFreads_inMillion'], Gene, 1951))
                except:
                    UppPredLimit_ContRead.append(expectedConR(GeneData.ix[i]['PFreads_inMillion'], Gene, 1801))
            elif 1952 <= GeneData.ix[i]['POS'] <= 2101:
                try:
                    UppPredLimit_ContRead.append(expectedConR(GeneData.ix[i]['PFreads_inMillion'], Gene, 2101))
                except:
                    UppPredLimit_ContRead.append(expectedConR(GeneData.ix[i]['PFreads_inMillion'], Gene, 1951))
            else:
                try:
                    UppPredLimit_ContRead.append(expectedConR(GeneData.ix[i]['PFreads_inMillion'], Gene, 2251))
                except:
                    UppPredLimit_ContRead.append(expectedConR(GeneData.ix[i]['PFreads_inMillion'], Gene, 2101))

        GeneData['UppPredLimit_ContRead'] = UppPredLimit_ContRead
        GeneData_tojoin = GeneData[(columnID[:-3]+'.AD')].str.split(",").apply(Series, 1)
        if len(GeneData_tojoin.columns) == 2:
            GeneData_tojoin[2] = np.NaN
        GeneData_tojoin.columns = ['refCopy', 'alt1Copy', 'alt2Copy']
        GeneData = pd.merge(GeneData, GeneData_tojoin, left_index=True, right_index=True)
        GeneData.refCopy = GeneData.refCopy.astype(float)
        GeneData.alt1Copy = GeneData.alt1Copy.astype(float)
        GeneData.alt2Copy = GeneData.alt2Copy.astype(float)
        GeneData = GeneData.set_index('POS')
        validityTester(GeneData)
        geneList.append(GeneData)

        return

    with open(reference_annotated_fa,'r') as refFile:
        sequences = refFile.read()

    segmentRegex = re.compile(r'((Segment\d)_[\S]*)')
    segment = segmentRegex.findall(sequences)
    for seg in segment:

        try:
            geneAnalysis(seg[0], seg[1])
        except:
            logging.info(seg[1] + ' of ' + columnID[:-3] + ' was not analysed')
            errorFile = open('Error.txt', 'a')
            errorFile.write(seg[1] + ' of ' + columnID[:-3] + ' was not analysed. Suggested solution: Something wrong '
                                                              'with the formating of the AllGenes.xls file' + '\n')
            errorFile.close()

    try:
        pd.concat(geneList, axis=0).to_excel(columnID[:-3] + "_AllGenes" + ".xls")
    except ValueError:
        print('All segments of sample ' + columnID[:-3] + ' were not analysed')
        logging.info('All segments of sample ' + columnID[:-3] + ' were not analysed. '
                    'Suggested solution: please use bwa-0.6.2 but not the latest version')
        errorFile = open('Error.txt', 'a')
        errorFile.write('All gene segments of sample ' + columnID[:-3] + ' were not analysed. '
                    'Suggested solution: please use bwa-0.6.2 but not the latest version' + '\n')
        errorFile.close()