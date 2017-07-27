README FOR AutomART 
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

CONTENTS OF THIS FILE
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
- Introduction to AutomART
- Dependencies (Pre-installation)
- Installation
- Configuration
- Running Auto mART
- Interpreting the results 
- Acknowledgments/Contact Info

INTRODUCTION TO AUTO mART
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
AutomART works to automate the process of filtering steps for remote homologues of known mono-ADP-ribosylating toxins. Currently, Auto mART utilizes "SignalP", "TMHMM" and "OB-score" to iterate through a directory containing several FASTA files of remote mART homologue sequences. AutomART then outputs a list of putative mART toxins for further analysis and their corresponding sequences in a FASTA file format. 

DEPENDANCIES (PRE-INSTALLATION)
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
1) Local Installation
* Linux or UNIX environment
* Perl 5.6 or higher
* SignalP
	- SignalP is available for download for free for academic users from: http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?signalp
	- Please follow SignalP README for installation instructions
* TMHMM
	- TMHMM is available for download for free for academic users from: http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?tmhmm
	- Please follow TMHMM README for installation instructions
* OB-score 
	- OB-score is available for download for free for academic users from:
	http://www.compbio.dundee.ac.uk/obscore/OB_1.0.tar.gz
	
	- Please follow OB-score README for installation instructions

INSTALLATION
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Local Installation 
1) Git Hub
- git clone https://github.com/latours/AutomART.git

2) Download Zip File
- 


CONFIGURATION
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
1) Open automart.sh (nano automart.sh)
2) Change "/full/path/signalp/" and "/full/path/tmhmm" to entire paths for the commands "signalp" and "tmhmm" respectively.
3) Change "/full/directory_path/OB to the path to OB-score's directory (Not the command itself).
4) Save respective changes

RUNNING AUTO mART 
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
*COMMAND OVERVIEW:
sh automart.sh -1 -2
-1 = Input directory path
-2 = Output directory path

*INPUT FILE FORMAT: 
AutomART takes a directory of FASTA files as input. 
The program takes the full path to the directory or subdirectory containing your list of FASTA files.
Input example: /home/user/FASTA_files (Do not add an additional "/" at the end of the directory name)
NOTE: Please ensure all FASTA files are in correct format before placing in directory. 

*OUTPUT FILE FORMAT:
Auto mART allows you to output all files generated by the program to a seperate directory from that containing your FASTA files. 
Auto mART takes the full path to the directory where you wish to save your output files. 
Command line example: /home/user/OUTPUT_files (Do not add a "/" at the end of the directory name)
NOTE: Please ensure the output directory exists before attempting to run Auto mART.
Make directory example: mkdir OUTPUT_files 

*EXAMPLE COMMAND FOR AUTO mART:

sh automart.sh /home/user/FASTA_files/ /home/user/OUTPUT_files/

(If there are permission restrictions use sudo before the command)

AutomART OUTPUT FORMAT
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Directories Produced:
1) GramPositive_Output
	- Description: This directory will contain all fasta and output files pertaining to tests conducted for gram positive 	      organisms
	- Files Produced:
2) GramNegative_Output
	- Description: This directory will contain all fasta and output files pertaining to tests conducted for gram positive 		organisms
	- Files Produced:
3) Secreted
	- Description: Each GramPositive/GramNegative Output directory contains a secreted subdirectory which holds the output 		files for all secreted (classical pathway) IDs as determined by the 		program SignalP
	- Files Produced:
4) Not_Secreted
	- Description:Each GramPositive/GramNegative Output directory contains a non secreted subdirectory. This directory 	   will contain a list of IDs that are not classically secreted as
	determined by the program SignalP. The output is also provided for user interpretation.
	- Files Produced:
5) Final_IDs
	- Description: Within each of the Secreted/Not_Secreted directories will be a subdirectory called Final IDs which will 		contain the final ID list of organisms and the corresponding sequences.
	- Files Produced:

Final Output File:
1) Final_Report.txt
	- Description of File Contents: This file contains a summary for the number of organism IDs processed by AutomART from 		the provided input directory. The file also contains a count for the number 		of IDs that passed the 		filters by AutomART analysis for both Gram Positive and Gram Negative organisms.

AKNOWLEDGMENTS/CONTACT INFO
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Developed by: Sara Latour M.BINF University of Guelph
Contact Information: saralatour@outlook.com
Publish Date: August 15th 2017 (08/15/17)
Acknowledgements:
Thanks to my research advisors Dr. Dan Ashlock and Dr. Rod Merrill for their guidance in developing this script.
Special thanks to Olivier Tremblay for his guidance as well

