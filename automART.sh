#!/bin/bash
###CREATOR: SARA LATOUR
###CONTACT: saralatour@outlook.com
###UPDATED:31/07/17
################################################################################
################################CUSTOMIZATION###################################

#####               Please Customize Program Pathways:                     #####
#### Path to signalp command:
signalp=/home/latours/Documents/BINF6999/signalp-4.1/signalp
#### Path to TMHMM command:
tmhmm=/home/latours/Documents/BINF6999/tmhmm-2.0c/bin/tmhmm
#### OB-Score Directory (Full Directory not command):
ob=/home/latours/Documents/BINF6999/OB/



################################################################################
#####             NOTHING BELOW THIS LINE SHOULD REQUIRE CHANGING:         #####
################################################################################

#Command Line Arguments
#$1 is the location of the FASTA files (Include the full path in command line)
#2 is the lcoation where files will be output (Include the full path in the command line)
inputdir=$1
outputdir=$2

cat << "EOF"


   _         _                     _    ____ _____
   / \  _   _| |_ ___  _ __ ___    / \  |  _ \_   _|
  / _ \| | | | __/ _ \| '_ ` _ \  / _ \ | |_) || |
 / ___ \ |_| | || (_) | | | | | |/ ___ \|  _ < | |
/_/   \_\__,_|\__\___/|_| |_| |_/_/   \_\_| \_\|_|

Created by: Sara Latour
Date: 26/07/17

EOF

if [ $# -ne 2 ]
  then
    echo "ERROR MESSAGE:\nPlease ensure that the full paths for the directories are provided to AutomART and in the correct order:\n(1) Directory of FASTA files \n(2) Output Directory\nNOTE: See README for a more detailed explanation on how to run AutomART"
	exit 1
fi

echo "Checking if Input directory is valid..."
if   [ -d "${inputdir}" ]
then echo "✓"
else echo "Input is not a valid directory... exiting AutomART now\n";
     exit 1
fi

echo "Checking if Output directory is valid..."
if   [ -d "${outputdir}" ]
then echo "✓\nRunning AutomART...\n"
else echo "Output is not a valid directory... exiting AutomART now\n";
     exit 1
fi



################################DIRECTORIES#####################################
################################################################################
#Gram Positive
rm -Rf $outputdir/GramPositive_Output; mkdir -p $outputdir/GramPositive_Output
rm -Rf $outputdir/GramPositive_Output/Secreted; mkdir -p $outputdir/GramPositive_Output/Secreted
rm -Rf $outputdir/GramPositive_Output/Secreted; mkdir -p $outputdir/GramPositive_Output/Secreted/Final_IDs
rm -Rf $outputdir/GramPositive_Output/Not_Secreted; mkdir -p $outputdir/GramPositive_Output/Not_Secreted

# Gram Negative
rm -Rf $outputdir/GramNegative_Output; mkdir -p $outputdir/GramNegative_Output
rm -Rf $outputdir/GramNegative_Output/Secreted; mkdir -p $outputdir/GramNegative_Output/Secreted
rm -Rf $outputdir/GramNegative_Output/Secreted; mkdir -p $outputdir/GramNegative_Output/Secreted/Final_IDs
rm -Rf $outputdir/GramNegative_Output/Not_Secreted; mkdir -p $outputdir/GramNegative_Output/Not_Secreted



##ID Count:

echo "Sequences Detected from FASTA Directory:"
orgcount=$(grep ">" -R $inputdir/ | wc -l)
echo "$orgcount\n"
cd $inputdir
echo "Now Processing FASTA files...\n"


## Make Fasta files into one line & Sort on mART motif:
echo "Now filtering sequences for mART template...\nNote:This step takes some time.\n"
for file in *.fasta #only selects fasta files
do
sed -e 's/\(^>.*$\)/#\1#/' $file | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' | egrep -B1 '[YFL]R.{27,60}[YF].ST[SQT].{32,78}[QE].E' | egrep -v '^--' >> $outputdir/new.fasta
done

lineonefasta=$outputdir/new.fasta
seq_template=$(wc -l $lineonefasta)

################################GRAM[+]TESTING##################################
################################################################################

echo  "Now Processing SignalP (Gram +)...\n"

$signalp -t gram+ -f short $lineonefasta > $outputdir/GramPositive_Output/signalp_pos.short_out 

#SignalP

cd $outputdir/GramPositive_Output/ #move to Gram Positive Directory

for file in *.short_out #for every signalp file do the following:
do
head -n2 "$file" | tail -n1 > Secreted/temp_"$file".output
tail -n+2 "$file" | awk -F" " '{if ($10=="Y") print $1}' >> Secreted/secreted_IDs.output 
tail -n+2 "$file" | awk -F" " '{if ($10=="N") print $1}' >> Not_Secreted/not_secreted_IDs.output #keep (All no from signal p = not secreted)
done
grampos_not_secreted=$(wc -l $outputdir/GramPositive_Output/Not_Secreted/not_secreted_IDs.output)
rm Secreted/temp_"$file".output
sort Secreted/secreted_IDs.output | uniq > Secreted/pos_secreted_IDs.output
grampos_secreted=$(wc -l $outputdir/GramPositive_Output/Secreted/pos_secreted_IDs.output)
grep -A1 -f Secreted/pos_secreted_IDs.output $lineonefasta > Secreted/grampositive_secreted.fasta #grep for all the secreted IDs to find their corresponding fasta sequence and save to file with all gram positive and secreted sequences.


#TMHMM

echo "Now running TMHMM for Gram Positive Secreted Sequences...\n"
$tmhmm Secreted/grampositive_secreted.fasta > Secreted/Final_IDs/TMHMM_grampositive_secreted.output #All gram + secreted sequences processed by TMHMM
grep "Number of predicted TMHs:" Secreted/Final_IDs/TMHMM_grampositive_secreted.output > Secreted/Final_IDs/temp.output #search for line with just the predicted number of TM domains
awk -F " " '{if ($7==0)print $2}' Secreted/Final_IDs/temp.output > Secreted/Final_IDs/grampositive_secreted_TMHMM_IDS.output
grampos_tmhmm=$(wc -l $outputdir/GramPositive_Output/Secreted/Final_IDs/grampositive_secreted_TMHMM_IDS.output)
rm Secreted/Final_IDs/temp.output #remove
grep -A1 -f Secreted/Final_IDs/grampositive_secreted_TMHMM_IDS.output $lineonefasta > Secreted/Final_IDs/grampositive_secreted_TMHMM_IDS.fasta #This gets the corresponding sequences to all IDs with no TM domains that are secreted and gram positive

#OB-SCORE

echo "Now running OB-Score for Gram Positive Secreted Sequences containing no Transmembrane Domains...\n"
cp Secreted/Final_IDs/grampositive_secreted_TMHMM_IDS.fasta $ob
cd $ob
perl OB.pl -i grampositive_secreted_TMHMM_IDS.fasta -o obscores.output -n -p /home/latours/Documents/BINF6999/OB/Hydrophobicity_scores.dat -m /home/latours/Documents/BINF6999/OB/zmat.dat
cp obscores.output $outputdir/GramPositive_Output/Secreted/Final_IDs/
cd $outputdir/GramPositive_Output/Secreted/Final_IDs/
awk '$2>=1.5' obscores.output | cut -f1  > pos_final_mART_IDs.txt
pos_finalcount=$(wc -l pos_final_mART_IDs.txt)
rm obscores.output 

#Create Final Fasta File

cd $outputdir/GramPositive_Output/Secreted/Final_IDs/
grep -A1 -f pos_final_mART_IDs.txt $lineonefasta > $outputdir/GramPositive_Output/Secreted/Final_IDs/temppos_final_mART_sequences.fasta
egrep -v '^--' $outputdir/GramPositive_Output/Secreted/Final_IDs/temppos_final_mART_sequences.fasta > $outputdir/GramPositive_Output/Secreted/Final_IDs/pos_final_mART_sequences.fasta
#Create NCBI Links File:
cat pos_final_mART_IDs.txt | while read line; do echo "https://www.ncbi.nlm.nih.gov/Structure/cdd/wrpsb.cgi?INPUT_TYPE=live&SEQUENCE=$line"; done > NCBI_pos_final_mART_links.txt
#Create CSV NCBI Links:
cat pos_final_mART_IDs.txt | while read line; do echo "https://www.ncbi.nlm.nih.gov/Structure/cdd/wrpsb.cgi?INPUT_TYPE=live&SEQUENCE=$line"; done > NCBI_pos_final_mART_links.csv

################################GRAM[-]TESTING##################################
################################################################################

echo  "Now Processing SignalP (Gram -)...\n"
$signalp -t gram- -f short $lineonefasta > $outputdir/GramNegative_Output/signalp_neg.short_out 

#SignalP
cd $outputdir/GramNegative_Output/ #move to Gram Negative Directory

for file in *.short_out #for every signalp file do the following:
do
head -n2 "$file" | tail -n1 > Secreted/temp_"$file".output
tail -n+2 "$file" | awk -F" " '{if ($10=="Y") print $1}' >> Secreted/secreted_IDs.output #keep (All yes from signal p = secreted)
tail -n+2 "$file" | awk -F" " '{if ($10=="N") print $1}' >> Not_Secreted/not_secreted_IDs.output #keep (All no from signal p = not secreted)
done
gramneg_not_secreted=$(wc -l $outputdir/GramNegative_Output/Not_Secreted/not_secreted_IDs.output)
rm Secreted/temp_"$file".output
sort Secreted/secreted_IDs.output | uniq > Secreted/neg_secreted_IDs.output
gramneg_secreted=$(wc -l $outputdir/GramNegative_Output/Secreted/neg_secreted_IDs.output)
grep -A1 -f Secreted/neg_secreted_IDs.output $lineonefasta > Secreted/gramnegative_secreted.fasta

#TMHMM#
echo "Now running TMHMM for Gram Negative Secreted Sequences...\n"
$tmhmm Secreted/gramnegative_secreted.fasta > Secreted/Final_IDs/TMHMM_gramnegative_secreted.output #All gram - secreted sequences processed by TMHMM
grep "Number of predicted TMHs:" Secreted/Final_IDs/TMHMM_gramnegative_secreted.output > Secreted/Final_IDs/temp.output #search for line with just the predicted number of TM domains
awk -F " " '{if ($7==0)print $2}' Secreted/Final_IDs/temp.output > Secreted/Final_IDs/gramnegative_secreted_TMHMM_IDS.output 
gramneg_tmhmm=$(wc -l $outputdir/GramNegative_Output/Secreted/Final_IDs/gramnegative_secreted_TMHMM_IDS.output)
rm Secreted/Final_IDs/temp.output
grep -A1 -f Secreted/Final_IDs/gramnegative_secreted_TMHMM_IDS.output $lineonefasta > Secreted/Final_IDs/gramnegative_secreted_TMHMM_IDS.fasta

#OB-SCORE#
echo "Now running OB-Score for Gram Negative Secreted Sequences containing no Transmembrane Domains...\n"
cp Secreted/Final_IDs/gramnegative_secreted_TMHMM_IDS.fasta $ob
cd $ob
perl OB.pl -i gramnegative_secreted_TMHMM_IDS.fasta -o obscores.output -n -p /home/latours/Documents/BINF6999/OB/Hydrophobicity_scores.dat -m /home/latours/Documents/BINF6999/OB/zmat.dat
cp obscores.output $outputdir/GramNegative_Output/Secreted/Final_IDs/
cd $outputdir/GramNegative_Output/Secreted/Final_IDs/
awk '$2>=1.5' obscores.output | cut -f1 > neg_final_mART_IDs.txt
neg_finalcount=$(wc -l neg_final_mART_IDs.txt)
rm obscores.output 

#Create Final Fasta File#
cd $outputdir/GramNegative_Output/Secreted/Final_IDs/
grep -A1 -f neg_final_mART_IDs.txt $lineonefasta > $outputdir/GramNegative_Output/Secreted/Final_IDs/tempneg_final_mART_sequences.fasta
egrep -v '^--'$outputdir/GramNegative_Output/Secreted/Final_IDs/neg_final_mART_sequences.fasta > $outputdir/GramNegative_Output/Secreted/Final_IDs/neg_final_mART_sequences.fasta

#Create NCBI Links File:
cat neg_final_mART_IDs.txt | while read line; do echo "https://www.ncbi.nlm.nih.gov/Structure/cdd/wrpsb.cgi?INPUT_TYPE=live&SEQUENCE=$line"; done > NCBI_neg_final_mART_links.txt
#Create CSV NCBI Links:
cat neg_final_mART_IDs.txt | while read line; do echo "https://www.ncbi.nlm.nih.gov/Structure/cdd/wrpsb.cgi?INPUT_TYPE=live&SEQUENCE=$line"; done > NCBI_neg_final_mART_links.csv

################################################################################
				#REMOVE TEMP FILES#
################################################################################
cd $outputdir
rm new.fasta
cd $outputdir/GramPositive_Output/
rm *.short_out
rm -R TMHMM*/
cd $outputdir/GramNegative_Output/
rm *.short_out
rm -R TMHMM*/

rm $outputdir/GramPositive_Output/Secreted/Final_IDs/grampositive_secreted_TMHMM_IDS.output
rm $outputdir/GramNegative_Output/Secreted/Final_IDs/gramnegative_secreted_TMHMM_IDS.output
rm $outputdir/GramPositive_Output/Secreted/grampositive_secreted.fasta
rm $outputdir/GramNegative_Output/Secreted/gramnegative_secreted.fasta
rm $outputdir/GramPositive_Output/Secreted/Final_IDs/TMHMM_grampositive_secreted.output 
rm $outputdir/GramNegative_Output/Secreted/Final_IDs/TMHMM_gramnegative_secreted.output 
rm $outputdir/GramPositive_Output/Secreted/Final_IDs/grampositive_secreted_TMHMM_IDS.fasta 
rm $outputdir/GramNegative_Output/Secreted/Final_IDs/gramnegative_secreted_TMHMM_IDS.fasta
rm $outputdir/GramPositive_Output/Secreted/secreted_IDs.output
rm $outputdir/GramNegative_Output/Secreted/secreted_IDs.output
rm $outputdir/GramNegative_Output/Secreted/Final_IDs/tempneg_final_mART_sequences.fasta
rm $outputdir/GramPositive_Output/Secreted/Final_IDs/temppos_final_mART_sequences.fasta

################################################################################
				#Create Final Report#
################################################################################
cd $outputdir
touch Final_Report.txt
echo "Number of Sequences Submitted to AutomART:$orgcount" > Final_Report.txt
echo "AutomART Results:" >> Final_Report.txt
echo "-----------------------------------------------\n" >> Final_Report.txt
echo "Number of Sequences that passed the mART sequence template:$seq_template">> Final_Report.txt
echo "Gram Positive # of Secreted IDs:$grampos_secreted" >> Final_Report.txt
echo "Gram Positive # of IDs Not Classically Secreted:$grampos_not_secreted" >> Final_Report.txt
echo "Gram Negative # of Secreted IDs:$gramneg_secreted" >> Final_Report.txt
echo "Gram Negative # of IDs Not Classically Secreted:$gramneg_not_secreted" >> Final_Report.txt
echo "Gram Positive # of Secreted IDs without TM Domains:$grampos_tmhmm" >> Final_Report.txt
echo "Gram Negative # of Secreted IDs without TM Domains:$gramneg_tmhmm" >> Final_Report.txt
echo "Gram Positive # of Final IDs:$pos_finalcount" >> Final_Report.txt
echo "\nGram Negative # of Final IDs:$neg_finalcount" >> Final_Report.txt

cat << "EOF"

  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-
 / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \
`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'

		      AutomARt Run Complete!

  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-
 / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \
`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'

EOF
