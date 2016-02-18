# FluSeq-v1.0

Version 1.0 for Illumina Miseq system.
For Research Use only. Not for use in diagnostic procedures.

FluSeq is a bioinformatics program developed and tested for an influenza A genome-wide amplicon-based high-throughput sequencing (HTS) method, using Illumina Miseq system [1]. It includes contigs assembly, blastn, filters for in-run contamination reads, and consensus sequence calling. Users are advised to refer the end-to-end laboratory protocol described by Lee, et al. 2016 [1]. The users are encouraged to modify the FluSeq program written in python scripts for other clinical virus amplicon-based HTS methods.

## INSTALLATION

Installation instructions are available for LINUX/UNIX or MacOSX. The entire analytic workflow was implemented and tested by the author on CentOS-6.6/RedHat Linux and MacOSX, but not on other operating systems. The operating system should contain an updated java.
	
### Installation of Python 3.4.3 
(https://www.python.org/downloads/)
- Upon installation, log in as root from Terminal (Linux and MacOSX): 
```
# pip install pandas
# pip install numpy
# pip install ZODB
# pip install lxml 
# pip install xlrd
# pip install xlwt
# pip install beautifulsoup4
# pip install scipy
# pip install transaction
```

###	Installation of VICUNA-v1.3 
(http://www.broadinstitute.org/scientific-community/science/projects/viral-genomics/vicuna) 
- The vicuna_config.txt used by this analytic workflow is stored in FluSeq Folder. Compilation of vicunAnalysis is not required during the installation. 
- path to the compiled vicuna executive file needs to be changed accordingly in the FluSeq-v1.0.py later, i.e. Line 100: First argument of the  subprocess.Popen, to “YourPathInstalled/VICUNA_v1.3/bin/vicuna-omp-v1.0”
- Optional: Line 104 - logging.info: Change program path to “YourPathInstalled/VICUNA_v1.3/bin/vicuna-omp-v1.0”.

### Installation of BLAST 2.2.31+ software 
(https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

### Installation of BWA-0.6.2 software 
(http://sourceforge.net/projects/bio-bwa/files/)
- Users are advised to install BWA-0.6.2 software but not the latest version, as bwa-aln/sampe in this version was found to be more stable in reads alignment for this study.

### Installation of samtools-1.2 
(http://sourceforge.net/projects/samtools/files/samtools/1.2/)

### Installation of GATK-3.4-46 software
(https://www.broadinstitute.org/gatk/download/)
- Path to the GenomeAnalysisTK.jar executive file needs to be changed accordingly in the FluSeq-v1.0.py later, i.e. Lines 233, 259,  and 269: Third or Fourth argument of the  subprocess.Popen, to “YourPathInstalled/GenomeAnalysisTK.jar”.
- Optional: Lines 243, 265,  and 278: logging.info: Change program path to “YourPathInstalled/GenomeAnalysisTK.jar”.
	
### Installation of picard-tools-1.138 
(https://github.com/broadinstitute/picard/releases/tag/1.138)
- Path to the GenomeAnalysisTK.jar executive file needs to be changed accordingly in the FluSeq-v1.0.py later, i.e. Line 146: Third argument of the  subprocess.Popen, to “YourPathInstalled/picard-tools-1.138/picard.jar”.
- Optional: Line 152: logging.info: Change program path to “YourPathInstalled/picard-tools-1.138/picard.jar”.

### Installation of FluSeq
a.	Download the FluSeq.tar, decompress, and mv the FluSeq Folder to ~.
```
$ tar –xvf  FluSeq.tar
$ mv FluSeq ~ 
```
b.	Change working directory to FluSeq 
```
$ cd ~/FluSeq
```	

