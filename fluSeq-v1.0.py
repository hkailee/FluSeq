#!/usr/bin/env python3.4
__author__ = 'mdc_hk'
version = '1.0'

import datetime, multiprocessing, os, re, shutil, sys, subprocess, time, logging

#===========================================================================================================

try:
    folderIn=sys.argv[1]
    proceed = True
    print('Program: FluSeq - Version: {} is starting now.\n'.format(version))
except:
    print('\nProgram: FluSeq\nVersion: {}\n'.format(version) + 'Usage: ' + sys.argv[0] + ' <input_folder> ')
    proceed = False

if proceed == True:

    # Setting up working and fluSeq directories...
    workingFolder_tmp = '~/FluSeq/' + folderIn
    os.chdir(os.path.expanduser(workingFolder_tmp))
    workingFolder = os.getcwd()
    fluSeqFolder = os.path.expanduser('~/FluSeq/')

    # Logging events...
    logging.basicConfig(filename= workingFolder + '/Log.txt', level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s')

    startTime = time.time()
    logging.info('Command-line invocation: ' + sys.argv[0] + ' ' + sys.argv[1])
    logging.info('Runfolder path: ' + workingFolder)

    # determining number of logical CPUs availalbe in your PC...
    numberOfProcessors = multiprocessing.cpu_count()
    print('The program will be optimized to use ' + str(numberOfProcessors) + ' logical CPUs available in your PC.')
    logging.info('Detected logical CPUs: ' + str(numberOfProcessors))

    # Filing number of unique samples found in the working folder...
    workingFilesR1 = [f for f in os.listdir(workingFolder) if re.match(r'[\w]+S\d+_L001_R1_001\.fastq\.gz', f)]
    fastqFileNameR1 = re.compile(r'(([\w]+)_S\d+_L001_R)1_001\.fastq\.gz')
    fastqFilesR1 = []
    for file in workingFilesR1:
        fastqFilesR1.append(fastqFileNameR1.findall(file))

    def delete_files(pth):
        for sub in os.listdir(pth):
            if os.path.isdir(pth + '/' + sub):
                None
            else:
                os.remove(pth + '/' + sub)

    # Creating folder for each sample and starting the analyses for each of them...
    for file in fastqFilesR1:
        # create if non-exists
        if not os.path.exists(workingFolder + '/' + file[0][1]):
            os.makedirs(workingFolder + '/' + file[0][1])
        else:
            # archive existing analysis if exists
            if not os.path.exists(workingFolder + '/' + file[0][1] + '/LastAnalysis'):
                os.makedirs(workingFolder + '/' + file[0][1] + '/LastAnalysis')
            suffixToCopy = ['_Reference_annotated.fa', '_Reference.fa', '.vcf', '_blastn.out', '_table.txt']
            try:
                for f in suffixToCopy:
                    shutil.copy(workingFolder + '/' + file[0][1] + '/' + file[0][1] + f, workingFolder + '/'
                                + file[0][1] + '/LastAnalysis')
                os.rename(workingFolder + '/' + file[0][1] + '_AllGenes.xls', workingFolder + '/'
                                + file[0][1] + '/LastAnalysis/' + file[0][1] + '_AllGenes.xls')
            except FileNotFoundError:
                None

            # remove existing analysis if exists
            delete_files(workingFolder + '/' + file[0][1])


        # copy fastq.gz to sample folder
        suffixFastqgz = ['1_001.fastq.gz', '2_001.fastq.gz']
        for i in suffixFastqgz:
            shutil.copy(file[0][0] + i, workingFolder + '/' + file[0][1])

        for f in suffixFastqgz:
            # unzip fastq.gz
            proc1 = subprocess.Popen(['gunzip', workingFolder + '/' + file[0][1] + '/' + file[0][0] + f],
                                     stdin=subprocess.PIPE, stdout=subprocess.PIPE)
            out, err = proc1.communicate()
            logging.info('gunzip ' + workingFolder + '/' + file[0][1] + '/' + file[0][0] + f)

        # Performing de novo assembly using VICUNA...
        print('>> '+ datetime.datetime.now().strftime('%Y/%m/%d %H:%M:%S'), '- Performing de-novo assembly for sample '
              + file[0][1])
        # copy vicuna_config.txt to sample folder
        shutil.copy(fluSeqFolder+'vicuna_config.txt', workingFolder+'/'+file[0][1])
        with open(workingFolder+'/'+file[0][1]+'/vicuna_config.txt', 'r') as configFile:
            configFileRead = configFile.read()
        # configure vicuna_config.txt
        newData = configFileRead.replace('/Folder001', workingFolder + '/' + file[0][1])
        with open(workingFolder+'/'+file[0][1]+'/vicuna_config.txt', 'w') as configFile:
            configFile.write(newData)
        # perform vicuna analysis
        os.environ['OMP_NUM_THREADS'] = '{:d}'.format(numberOfProcessors)
        proc2 = subprocess.Popen(['/home/hklee/Software/VICUNA_v1.3/bin/vicuna-omp-v1.0',
                                  workingFolder+'/'+file[0][1]+'/vicuna_config.txt'], stdin=subprocess.PIPE,
                                 stdout=subprocess.PIPE)
        out, err = proc2.communicate()
        logging.info('OMP_NUM_THREADS=' + str(numberOfProcessors)
                     + ' /home/hklee/Software/VICUNA_v1.3/bin/vicuna-omp-v1.0 '
                     + workingFolder+'/'+file[0][1]+'/vicuna_config.txt')

        # Performing blastn for each sample
        print('>> '+ datetime.datetime.now().strftime('%Y/%m/%d %H:%M:%S'), '- Performing blastn for sample '
              + file[0][1])
        proc3 = subprocess.Popen(['blastn', '-db', 'FASTA.fa', '-task', 'megablast',
                                  '-query', workingFolder+ '/' + file[0][1] + '/contig.fasta',
                                  '-out', workingFolder+ '/' + file[0][1] + '/' + file[0][1] + '_blastn.out',
                                  '-num_threads', str(numberOfProcessors), '-outfmt', '6', '-max_target_seqs', '15',
                                  '-use_index', 'True'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        out, err = proc3.communicate()
        logging.info('blastn -db FASTA.fa -task megablast -query ' + workingFolder + '/' + file[0][1] + '/contig.fasta '
                        '-out ' + workingFolder + '/' + file[0][1] + '/' + file[0][1] + '_blastn.out -num_threads ' +
                     str(numberOfProcessors) + ' -outfmt 6 -max_target_seqs 15 -use_index True')

        # Parsing blastn results to retrieve fasta sequences from the database for each sample; Annotating the fasta
        # sequences
        print('>> '+ datetime.datetime.now().strftime('%Y/%m/%d %H:%M:%S'), '- Retriving the fasta sequences for '
                'sample ' + file[0][1] + ' from the database, according to the blastn results')
        proc4 = subprocess.Popen(['python3.4', fluSeqFolder + 'ncbiParser.py',
                                  workingFolder+ '/' + file[0][1] + '/' + file[0][1] + '_blastn.out', file[0][1],
                                  workingFolder+ '/' + file[0][1]], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        out, err = proc4.communicate()
        logging.info('python3.4 ' + fluSeqFolder + 'ncbiParser.py ' + workingFolder+ '/' + file[0][1] + '/' + file[0][1]
                         + '_blastn.out ' + file[0][1] + ' ' + workingFolder+ '/' + file[0][1])

        # Generating bwa index, fasta file index, sequence dictionary for reference_annotated.fa of each sample
        print('>> '+ datetime.datetime.now().strftime('%Y/%m/%d %H:%M:%S'), '- Generating the bwa index, fasta file '
                'index, and sequence dictionary for reference.fa of sample ' + file[0][1])
        proc5 = subprocess.Popen(['bwa', 'index', workingFolder + '/' + file[0][1] + '/' + file[0][1] +
                                  '_Reference_annotated.fa'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        out, err = proc5.communicate()
        logging.info('bwa index ' + workingFolder + '/' + file[0][1] + '/' + file[0][1] + '_Reference_annotated.fa')

        proc6 = subprocess.Popen(['samtools', 'faidx', workingFolder + '/' + file[0][1] + '/' + file[0][1] +
                                  '_Reference_annotated.fa'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        out, err = proc6.communicate()
        logging.info('samtools faidx ' + workingFolder + '/' + file[0][1] + '/' + file[0][1]
                     + '_Reference_annotated.fa')

        proc7 = subprocess.Popen(['java', '-jar', '/home/hklee/Software/picard-tools-1.138/picard.jar',
                                  'CreateSequenceDictionary', 'REFERENCE='+ workingFolder + '/' + file[0][1] + '/'
                                  + file[0][1] + '_Reference_annotated.fa',
                                  'OUTPUT=' + workingFolder + '/' + file[0][1] + '/' + file[0][1] +
                                  '_Reference_annotated.dict'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        out, err = proc7.communicate()
        logging.info('java -jar /home/hklee/Software/picard-tools-1.138/picard.jar CreateSequenceDictionary REFERENCE='
                     + workingFolder + '/' + file[0][1] + '/' + file[0][1] + '_Reference_annotated.fa OUTPUT='
                     + workingFolder + '/' + file[0][1] + '/' + file[0][1] +'_Reference_annotated.dict')


        # bwa aln...
        print('>> '+ datetime.datetime.now().strftime('%Y/%m/%d %H:%M:%S'), '- Performing bwa aln for paired-end reads '
                                                                            'of sample ' + file[0][1])

        proc8 = subprocess.Popen(['bwa', 'aln', '-t', str(numberOfProcessors), '-q', '15', '-f',
                                  workingFolder + '/' + file[0][1] + '/' + file[0][0] + '1_001.sai',
                                  workingFolder + '/' + file[0][1] + '/' + file[0][1] + '_Reference_annotated.fa',
                                  workingFolder + '/' + file[0][1] + '/' + file[0][0]+'1_001.fastq'],
                                 stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        out, err = proc8.communicate()
        logging.info('bwa aln -t ' + str(numberOfProcessors) + ' -q 15 -f ' + workingFolder + '/' + file[0][1] + '/'
                     + file[0][0] + '1_001.sai ' + workingFolder + '/' + file[0][1] + '/' + file[0][1]
                     + '_Reference_annotated.fa ' + workingFolder + '/' + file[0][1] + '/' + file[0][0]+'1_001.fastq')

        proc9 = subprocess.Popen(['bwa', 'aln', '-t', str(numberOfProcessors), '-q', '15', '-f',
                                  workingFolder + '/' + file[0][1] + '/' + file[0][0] + '2_001.sai',
                                  workingFolder + '/' + file[0][1] + '/' + file[0][1] + '_Reference_annotated.fa',
                                  workingFolder + '/' + file[0][1] + '/' + file[0][0]+'2_001.fastq'],
                                 stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        out, err = proc9.communicate()
        logging.info('bwa aln -t ' + str(numberOfProcessors) + ' -q 15 -f ' + workingFolder + '/' + file[0][1] + '/'
                     + file[0][0] + '2_001.sai ' + workingFolder + '/' + file[0][1] + '/' + file[0][1]
                     + '_Reference_annotated.fa ' + workingFolder + '/' + file[0][1] + '/' + file[0][0]+'2_001.fastq')

        # bwa sampe...
        print('>> '+ datetime.datetime.now().strftime('%Y/%m/%d %H:%M:%S'), '- Performing bwa sampe for sample ' +
              file[0][1])

        proc10 = subprocess.Popen(['bwa', 'sampe', '-r' + '@RG\tID:'+file[0][1]+'\tPL:ILLUMINA\tSM:'+file[0][1],
                                  '-f', workingFolder + '/' + file[0][1] + '/' + file[0][1] + '.uncompressed.bam',
                                  workingFolder + '/' + file[0][1] + '/' + file[0][1] + '_Reference_annotated.fa',
                                  workingFolder + '/' + file[0][1] + '/' + file[0][0] + '1_001.sai',
                                  workingFolder + '/' + file[0][1] + '/' + file[0][0] + '2_001.sai',
                                  workingFolder + '/' + file[0][1] + '/' + file[0][0] + '1_001.fastq',
                                  workingFolder + '/' + file[0][1] + '/' + file[0][0] + '2_001.fastq'],
                                 stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        out, err = proc10.communicate()
        logging.info('bwa sampe -r @RG\tID:'+file[0][1]+'\tPL:ILLUMINA\tSM:'+file[0][1] +
                     ' -f '+ workingFolder + '/' + file[0][1] + '/' + file[0][1] + '.uncompressed.bam ' +
                     workingFolder + '/' + file[0][1] + '/' + file[0][1] + '_Reference_annotated.fa ' +
                     workingFolder + '/' + file[0][1] + '/' + file[0][0] + '1_001.sai ' +
                     workingFolder + '/' + file[0][1] + '/' + file[0][0] + '2_001.sai ' +
                     workingFolder + '/' + file[0][1] + '/' + file[0][0] + '1_001.fastq ' +
                     workingFolder + '/' + file[0][1] + '/' + file[0][0] + '2_001.fastq')

        # Performing bam sorting using samtools...
        print('>> '+ datetime.datetime.now().strftime('%Y/%m/%d %H:%M:%S'), '- Performing bam sorting for sample '
              + file[0][1] + ' using samtools sort module')
        proc11 = subprocess.Popen(['samtools', 'sort', workingFolder + '/' + file[0][1] + '/' + file[0][1] +
                                  '.uncompressed.bam', workingFolder + '/' + file[0][1] + '/' + file[0][1]],
                                 stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        out, err = proc11.communicate()
        logging.info('samtools sort ' + workingFolder + '/' + file[0][1] + '/' + file[0][1] +
                     '.uncompressed.bam ' + workingFolder + '/' + file[0][1] + '/' + file[0][1])

        # Performing bam indexing using samtools...
        print('>> '+ datetime.datetime.now().strftime('%Y/%m/%d %H:%M:%S'), '- Performing samtools indexing for sample '
              + file[0][1] + ' using samtools index module')
        proc12 = subprocess.Popen(['samtools', 'index', workingFolder + '/' + file[0][1] + '/' + file[0][1]+'.bam'],
                                 stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        out, err = proc12.communicate()
        logging.info('samtools index ' + workingFolder + '/' + file[0][1] + '/' + file[0][1]+'.bam')

        # Housekeeping...
        os.unlink(workingFolder + '/' + file[0][1] + '/' + file[0][1] + '.uncompressed.bam')

        # creating list of gene segments via reference sequences..
        with open(workingFolder + '/' + file[0][1] + '/' + file[0][1] + '_Reference_annotated.fa','r') as refFile:
            sequences = refFile.read()
        segmentRegex = re.compile(r'((Segment\d)_[\S]*)')
        segment = segmentRegex.findall(sequences)

        # GATK UnifiedGenotyper to generate vcf for each genes (variant calling)...
        print('>> '+ datetime.datetime.now().strftime('%Y/%m/%d %H:%M:%S'), '- Performing variant calling for each '
                'gene segments of sample ' + file[0][1] + ' using GATK UnifiedGenotyper')
        def proc13(geneAndLength, geneSuffix):
            subproc = subprocess.Popen(['java', '-Xmx4g', '-jar', '/home/hklee/Software/GenomeAnalysisTK.jar', '-T',
                                        'UnifiedGenotyper', '-glm', 'BOTH', '-ploidy', '3', '-R', workingFolder + '/'
                                        + file[0][1] + '/' + file[0][1] + '_Reference_annotated.fa',
                                        '-I', workingFolder + '/' + file[0][1] + '/' + file[0][1]+'.bam',
                                        '-dcov', '10000000', '-l', str('INFO'), '-stand_emit_conf', '10',
                                        '-stand_call_conf', '10', '-nt', str(numberOfProcessors),
                                        '-L', str(geneAndLength),
                                        '-o', workingFolder + '/' + file[0][1] + '/' + 'tmp_'+file[0][1]+geneSuffix],
                                       stdin=subprocess.PIPE, stdout=subprocess.PIPE)
            out, err = subproc.communicate()
            logging.info('java -Xmx4g -jar /home/hklee/Software/GenomeAnalysisTK.jar -T UnifiedGenotyper -glm BOTH '
                         '-ploidy 3 -R ' + workingFolder + '/' + file[0][1] + '/' + file[0][1] +
                         '_Reference_annotated.fa ' + '-I ' + workingFolder + '/' + file[0][1] + '/' + file[0][1] +
                         '.bam -dcov 10000000 -l ' + 'INFO' + ' -stand_emit_conf 10 -stand_call_conf 10 -nt '
                         + str(numberOfProcessors) + ' -L ' + str(geneAndLength) + ' -o ' + workingFolder + '/'
                         + file[0][1] + '/' + 'tmp_' + file[0][1] +geneSuffix)

        for seg in segment:
            proc13(seg[0], '_' + seg[1]+'.vcf')

        # Listing tmp_Segment.vcf files for CatVariants tool
        tmpVCF = []
        for seg in segment:
            list1 = ['-V', workingFolder + '/' + file[0][1] + '/tmp_' + file[0][1] + '_' + seg[1] + '.vcf']
            tmpVCF = tmpVCF + list1

        proc14 = subprocess.Popen(['java', '-cp', '/home/hklee/Software/GenomeAnalysisTK.jar',
                                  'org.broadinstitute.gatk.tools.CatVariants',
                                   '-R', workingFolder + '/'+ file[0][1] + '/' + file[0][1] + '_Reference_annotated.fa']
                                  + tmpVCF + ['-out', workingFolder + '/' + file[0][1] + '/' + file[0][1]+'.vcf',
                                    '-assumeSorted'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        out, err = proc14.communicate()
        logging.info('java -cp /home/hklee/Software/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants -R '
                     + workingFolder + '/' + file[0][1] + '/' + file[0][1]+'_Reference_annotated.fa ' + ' '.join(tmpVCF)
                     + ' -out ' + workingFolder + '/' + file[0][1] + '/' + file[0][1]+'.vcf -assumeSorted')

        proc15 = subprocess.Popen(['java', '-jar', '/home/hklee/Software/GenomeAnalysisTK.jar',
                                   '-R', workingFolder + '/' + file[0][1] +'/'+ file[0][1] + '_Reference_annotated.fa',
                                   '-T', 'VariantsToTable',
                                   '-V', workingFolder + '/' + file[0][1] + '/' + file[0][1]+'.vcf',
                                   '-F', 'CHROM','-F', 'POS', '-F', 'REF', '-F', 'ALT', '-F', 'QUAL', '-F', 'FILTER',
                                   '-F', 'DP', '-GF', 'AD',
                                   '-o', workingFolder + '/' + file[0][1] + '/' + file[0][1]+'_table.txt'],
                                  stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        out, err = proc15.communicate()
        logging.info('java -jar /home/hklee/Software/GenomeAnalysisTK.jar -R '+ workingFolder + '/' + file[0][1] + '/'
                     + file[0][1] + '_Reference_annotated.fa -T VariantsToTable -V '+ workingFolder + '/' + file[0][1]
                     + '/' + file[0][1]+'.vcf -F CHROM -F POS -F REF -F ALT -F QUAL -F FILTER -F DP -GF AD -o '
                     + workingFolder + '/' + file[0][1] + '/' + file[0][1]+'_table.txt')

        # Housekeeping...
        suffixGeneVcfUnlink = ['_Segment1.vcf', '_Segment2.vcf', '_Segment3.vcf', '_Segment4.vcf', '_Segment5.vcf',
                               '_Segment6.vcf', '_Segment7.vcf', '_Segment8.vcf', '_Segment1.vcf.idx',
                               '_Segment2.vcf.idx', '_Segment3.vcf.idx', '_Segment4.vcf.idx', '_Segment5.vcf.idx',
                               '_Segment6.vcf.idx', '_Segment7.vcf.idx', '_Segment8.vcf.idx']

        for i in suffixGeneVcfUnlink:
            try:
                os.unlink(workingFolder + '/' + file[0][1] + '/' + 'tmp_' + file[0][1] + i)
            except FileNotFoundError:
                print(workingFolder + '/' + file[0][1] + '/' + 'tmp_' + file[0][1] + i + 'is not found.')

        # Verifying the validity of variants found via GATK UnifiedGenotyper (i.e. excluding contamination/background
        # noise)...
        print('>> '+ datetime.datetime.now().strftime('%Y/%m/%d %H:%M:%S'), '- Sample:' + file[0][1] + '- Verifying the'
                ' validity of variants found via GATK UnifiedGenotyper (i.e. excluding contamination/background noise)')
        proc16 = subprocess.Popen(['python3.4', fluSeqFolder + 'varValidity.py',
                                   workingFolder + '/' + file[0][1] + '/' + file[0][1]+'_table.txt',
                                   workingFolder + '/' + 'ResequencingRunStatistics.xml',
                                   workingFolder + '/' + file[0][1] + '/' + file[0][1] + '_Reference_annotated.fa',
                                   str(0.9999)], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        out, err = proc16.communicate()
        print(out.decode('utf-8'))
        logging.info('python3.4 ' + fluSeqFolder + 'varValidity.py ' +
                     workingFolder + '/' + file[0][1] + '/' + file[0][1]+'_table.txt ' +
                     workingFolder + '/' + 'ResequencingRunStatistics.xml ' +
                     workingFolder + '/' + file[0][1] + '/' + file[0][1] + '_Reference_annotated.fa 0.9900')

        print('>> '+ datetime.datetime.now().strftime('%Y/%m/%d %H:%M:%S'), '- Sample:' + file[0][1] + '- Generating '
                'consensus sequences from verified variants')
        proc17 = subprocess.Popen(['python3.4', fluSeqFolder + 'seqGen.py',
                                   workingFolder + '/' + file[0][1]+'_AllGenes.xls',
                                   workingFolder + '/' + file[0][1] + '/' + file[0][1] + '_Reference_annotated.fa'],
                                  stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        out, err = proc17.communicate()
        logging.info('python3.4 ' + fluSeqFolder + 'seqGen.py ' +
                     workingFolder + '/' + file[0][1]+'_AllGenes.xls ' +
                     workingFolder + '/' + file[0][1] + '/' + file[0][1] + '_Reference_annotated.fa')

    endTime = time.time()
    print('>> '+ datetime.datetime.now().strftime('%Y/%m/%d %H:%M:%S'), '- Folder run completed: Took {} seconds to '
                                                                        'complete.'.format(endTime-startTime))
