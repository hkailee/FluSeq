#!/usr/bin/env python3.4
__author__ = 'mdc_hk'
version = '1.0'
#===========================================================================================================
import datetime, logging, multiprocessing, os, re, subprocess, sys, time
import pandas as pd
from pandas import Series
from bs4 import BeautifulSoup
#===========================================================================================================
# Functions:

# 1: Checks if in proper number of arguments are passed gives instructions on proper use.
def argsCheck(numArgs):
	if len(sys.argv) < numArgs or len(sys.argv) > numArgs:
		print('To learn about contamination rate from system control')
		print('Usage:', sys.argv[0], '<FolderInput>',' <DateOfRunInYYMMDD>')
		print('Examples:', sys.argv[0], 'FolderMS2_001', '151225')
		exit(1) # Aborts program. (exit(1) indicates that an error occurred)
#===========================================================================================================
# Housekeeping.
argsCheck(3) # Checks if the number of arguments are correct.

# Stores file one for input checking.
inFolder  = sys.argv[1]
dateOfRun = sys.argv[2]

# Setting up working and fluSeq directories...
workingFolder_tmp = '~/FluSeq/' + inFolder
os.chdir(os.path.expanduser(workingFolder_tmp))
workingFolder = os.getcwd()
fluSeqFolder = os.path.expanduser('~/FluSeq/')

# Logging events...
logging.basicConfig(filename=workingFolder + '/Log.txt', level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')
startTime = time.time()
logging.info('Command-line invocation: ' + sys.argv[0] + ' ' + sys.argv[1])
logging.info('Runfolder path: ' + workingFolder)

# determining number of logical CPUs availalbe in your PC...
numberOfProcessors = multiprocessing.cpu_count()
print('>> '+ datetime.datetime.now().strftime('%Y/%m/%d %H:%M:%S'), '- The program will be optimized to use '
      + str(numberOfProcessors) + ' logical CPUs available in your PC.')
logging.info('Detected logical CPUs: ' + str(numberOfProcessors))

# Filing number of unique samples found in the working folder...
workingFilesR1 = [f for f in os.listdir(workingFolder) if re.match(r'[\S]+S\d+_L001_R1_001\.fastq\.gz', f)]
fastqFileNameR1 = re.compile(r'(([\S]+)_S\d+_L001_R)1_001\.fastq\.gz')
fastqFilesR1 = []
for file in workingFilesR1:
    fastqFilesR1.append(fastqFileNameR1.findall(file))

# Starting the analyses...
for file in fastqFilesR1:
    # bwa aln...
    print('>> '+ datetime.datetime.now().strftime('%Y/%m/%d %H:%M:%S'), '- Performing bwa aln for paired-end reads '
                                                    'of sample ' + file[0][1])
    proc1 = subprocess.Popen(['bwa', 'aln', '-t', str(numberOfProcessors), '-q', '15', '-f',
                                  workingFolder + '/' + file[0][0] + '1_001.sai',
                                  fluSeqFolder + 'GENOMERepository/H3N2Genome_annotated.fa',
                                  workingFolder + '/' + file[0][0]+'1_001.fastq.gz'],
                                 stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    out, err = proc1.communicate()
    logging.info('bwa aln -t ' + str(numberOfProcessors) + ' -q 15 -f ' + workingFolder + '/'
                    + file[0][0] + '1_001.sai ' + fluSeqFolder + 'GENOMERepository/H3N2Genome_annotated.fa '
                    + workingFolder + '/' + file[0][0]+'1_001.fastq.gz')

    proc2 = subprocess.Popen(['bwa', 'aln', '-t', str(numberOfProcessors), '-q', '15',
                              '-f', workingFolder + '/' + file[0][0] + '2_001.sai',
                              fluSeqFolder + 'GENOMERepository/H3N2Genome_annotated.fa',
                              workingFolder + '/' + file[0][0]+'2_001.fastq.gz'],
                             stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    out, err = proc2.communicate()
    logging.info('bwa aln -t ' + str(numberOfProcessors) + ' -q 15 -f ' + workingFolder + '/'
                     + file[0][0] + '2_001.sai ' + fluSeqFolder + 'GENOMERepository/H3N2Genome_annotated.fa '
                     + workingFolder + '/' + file[0][0]+'2_001.fastq.gz')

#     bwa sampe...
    print('>> '+ datetime.datetime.now().strftime('%Y/%m/%d %H:%M:%S'), '- Performing bwa sampe for sample ' +
              file[0][1])
    proc3 = subprocess.Popen(['bwa', 'sampe', '-r' + '@RG\tID:'+file[0][1]+'\tPL:ILLUMINA\tSM:'+file[0][1],
                               '-f', workingFolder + '/' + file[0][1] + '.uncompressed.bam',
                               fluSeqFolder + 'GENOMERepository/H3N2Genome_annotated.fa',
                               workingFolder + '/' + file[0][0] + '1_001.sai',
                               workingFolder + '/' + file[0][0] + '2_001.sai',
                               workingFolder + '/' + file[0][0] + '1_001.fastq.gz',
                               workingFolder + '/' + file[0][0] + '2_001.fastq.gz'],
                              stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    out, err = proc3.communicate()
    logging.info('bwa sampe -r @RG\tID:'+file[0][1]+'\tPL:ILLUMINA\tSM:'+file[0][1] +
                     ' -f '+ workingFolder + '/' + file[0][1] + '.uncompressed.bam ' +
                     fluSeqFolder + 'GENOMERepository/H3N2Genome_annotated.fa ' +
                     workingFolder + '/' + file[0][0] + '1_001.sai ' +
                     workingFolder + '/' + file[0][0] + '2_001.sai ' +
                     workingFolder + '/' + file[0][0] + '1_001.fastq.gz ' +
                     workingFolder + '/' + file[0][0] + '2_001.fastq.gz')

    # Performing bam sorting using samtools...
    print('>> '+ datetime.datetime.now().strftime('%Y/%m/%d %H:%M:%S'), '- Performing bam sorting for sample '
              + file[0][1] + ' using samtools sort module')
    proc4 = subprocess.Popen(['samtools', 'sort', workingFolder + '/' + file[0][1] +
                               '.uncompressed.bam', workingFolder + '/' + file[0][1]],
                                stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    out, err = proc4.communicate()
    logging.info('samtools sort ' + workingFolder + '/' + file[0][1] +
                     '.uncompressed.bam ' + workingFolder + '/' + file[0][1])

    # Performing bam indexing using samtools...
    print('>> '+ datetime.datetime.now().strftime('%Y/%m/%d %H:%M:%S'), '- Performing samtools indexing for sample '
              + file[0][1] + ' using samtools index module')
    proc5 = subprocess.Popen(['samtools', 'index', workingFolder + '/' + file[0][1]+'.bam'],
                              stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    out, err = proc5.communicate()
    logging.info('samtools index ' + workingFolder + '/' + file[0][1]+'.bam')

    print('>> '+ datetime.datetime.now().strftime('%Y/%m/%d %H:%M:%S'), '- Analysing Depth of Coverage for each '
            'gene segments of sample ' + file[0][1] + ' using GATK DepthOfCoverage')
    proc6 = subprocess.Popen(['java', '-Xmx4g', '-jar', '/home/hklee/Software/GenomeAnalysisTK.jar', '-T',
                              'DepthOfCoverage', '-R', fluSeqFolder + 'GENOMERepository/H3N2Genome_annotated.fa',
                              '-o', workingFolder + '/' + 'MS2SysCtrl_base',
                              '-I', workingFolder + '/' + file[0][1]+'.bam', '-omitIntervals',
                              '-omitSampleSummary'],
                              stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    out, err = proc6.communicate()
    logging.info('java -Xmx4g -jar /home/hklee/Software/GenomeAnalysisTK.jar -T DepthOfCoverage -R'
                 + fluSeqFolder + 'GENOMERepository/H3N2Genome_annotated.fa -o ' + workingFolder + '/'
                 + 'MS2SysCtrl_base -I ' + workingFolder + '/' + file[0][1]
                 + '.bam -omitIntervals -omitSampleSummary')

    # Housekeeping...
    os.unlink(workingFolder + '/' + file[0][0] + '1_001.sai')
    os.unlink(workingFolder + '/' + file[0][0] + '2_001.sai')
    os.unlink(workingFolder + '/' + file[0][1] + '.uncompressed.bam')
    os.unlink(workingFolder + '/' + file[0][1]+'.bam')
    os.unlink(workingFolder + '/' + file[0][1]+'.bam.bai')

dataFromTable = pd.read_table(workingFolder + '/MS2SysCtrl_base', sep='\t')
columns_tojoin = dataFromTable[('Locus')].str.split(":").apply(Series, 1)
columns_tojoin.columns = ['CHROM', 'POS']
dataFromTable = pd.merge(dataFromTable, columns_tojoin, left_index=True, right_index=True)
dataFromTable = dataFromTable.set_index(['CHROM'], drop=False)

try:
    bsObj = BeautifulSoup(open(workingFolder + '/ResequencingRunStatistics.xml'), 'lxml-xml')

except IOError:
    errorFile = open('Error.txt', 'a')
    errorFile.write('Please make sure the run-specific ResequencingRunStatistics.xml is placed in the run folder'
                    + '\n')
    errorFile.close()
    exit(1)

xmlOfSamples = bsObj.findAll('SummarizedSampleStatistics')
if not xmlOfSamples:
    xmlOfSamples = bsObj.findAll('SummarizedSampleStatisics')

listOfSamples = []
for name in xmlOfSamples:
    listOfSamples.append(name.SampleName.get_text())
listOfPFReads = []
for name in xmlOfSamples:
    listOfPFReads.append(round(int(name.NumberOfClustersPF.get_text())/1000000, 2))
if len(listOfSamples) == len(listOfPFReads):
    dictOfSamplesAndPFreads = dict(zip(listOfSamples, listOfPFReads))

geneList = []

def geneAnalysis(Chrom, Gene):
    GeneData = dataFromTable.ix[Chrom]
    GeneData.is_copy = False
    GeneData['PFreads_inMillion'] = dictOfSamplesAndPFreads['MS2SysCtrl']
    GeneData['DateOfRun'] = dateOfRun
    GeneData['GeneSegment'] = Gene
    GeneData = GeneData.set_index(['POS'], drop=False)
    geneSeq = list(dictOfGeneSegment[Chrom])
    for k in range(0, len(geneSeq), 150):
        GeneData_a = GeneData.ix[[k],['DateOfRun', 'CHROM', 'GeneSegment', 'POS', 'Total_Depth', 'PFreads_inMillion']]
        geneList.append(GeneData_a)
    return


with open(fluSeqFolder + 'GENOMERepository/H3N2Genome_annotated.fa' ,'r') as refFile:
    content1 = refFile.read().splitlines()


geneName = []
for i in range(len(content1)):
    if content1.index(content1[i]) % 2 == 0 or content1.index(content1[i]) == 0:
        geneName.append(content1[i][1:])

geneCont = []
for i in range(len(content1)):
    if content1.index(content1[i]) % 2 != 0:
        geneCont.append(content1[i])

if len(geneName) == len(geneCont):
    dictOfGeneSegment = dict(zip(geneName, geneCont))

with open(fluSeqFolder + 'GENOMERepository/H3N2Genome_annotated.fa','r') as refFile:
    sequences = refFile.read()

segmentRegex = re.compile(r'((Segment\d)_[\S]*)')
segment = segmentRegex.findall(sequences)

for seg in segment:
    try:
        geneAnalysis(seg[0], seg[1])
    except:
        errorFile = open('Error.txt', 'a')
        errorFile.write(seg[1] + ' of MS2 System Control was not analysed. Suggested solution: Something is wrong '
                                 'with the formatting of the AllGenes.xls file' + '\n')
        errorFile.close()

try:
    print('>> '+ datetime.datetime.now().strftime('%Y/%m/%d %H:%M:%S'), '- Updating existing ms2SysCtrl.db')
    dataBase = pd.read_csv(fluSeqFolder + 'ms2SysCtrl.db', sep=',')
    dataBase_tmp = pd.concat(geneList, axis=0)
    with open(fluSeqFolder + 'ms2SysCtrl.db', 'a') as f:
        dataBase_tmp.to_csv(f, header=False)

except OSError:
    print('Warning: The MS2 contamination database file is not found in the designated directory. '
          'A new database will be created. Please note that it is important to have sufficient data in the database '
          'to provide a confident contamination statistics for data analyses')
    errorFile = open('Error.txt', 'a')
    errorFile.write('Warning: The MS2 contamination database file is not found in the designated directory. A new '
                    'database will be created. Please note that it is important to have sufficient data in the '
                    'database to provide a confident contamination statistics for data analyses.\n')
    errorFile.close()
    pd.concat(geneList, axis=0).to_csv(fluSeqFolder + 'ms2SysCtrl.db')

# Housekeeping...
try:
    os.unlink(workingFolder + '/MS2SysCtrl_base.sample_cumulative_coverage_counts')
except OSError:
    print(workingFolder + '/MS2SysCtrl_base.sample_cumulative_coverage_counts is not found')
try:
    os.unlink(workingFolder + '/MS2SysCtrl_base.sample_cumulative_coverage_proportions')
except OSError:
    print(workingFolder + '/MS2SysCtrl_base.sample_cumulative_coverage_proportions is not found')

print('>> '+ datetime.datetime.now().strftime('%Y/%m/%d %H:%M:%S'), '- MS2 System contamination learning done ')


