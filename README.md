# pepMID_simul

## Descriptions
This package contains the following:

"peptide_mid.m" simulates the peptide mass isotopologue distribution (MID) for a given peptide in string format (e.g., 'AKRG'), and outputs the m/z and abundance for each isotopologue. The natural abundance for 13C, 15N, 2H, 17O, 18O, 32S, 33S and 36S are considered. The results have been tested to be fully consistent with the online tool https://prospector.ucsf.edu/prospector/cgi-bin/msform.cgi?form=msisotope 

"Main_fitting" determines the peptide turnover rate (r) in batch mode, given the measured MID of newly synthesized peptide and MID for 2D labeled amino acide. 
The measured isotope distribution is fitted to the linear equation below by least square fitting:

    MID_measured = MID_unlabeled * r + MID_labeled * (1 - r)
    
MID_labeled is calculated by the convolution of MID_unlabeled with the MID of 2D labeled peptide of its labeled amino acid composition.

other subroutines: 
"pep2mass.m" calculates the monoisotopic mass of a peptide.
"mergeM.m" finds the MID of the peptide from the measured MID of labeled amino acid by convolution
"formula2mass.m"  finds the monoisotopic mass from a given formula

## Usage
Main_fitting requires two input files in csv formats, the first is the measured amino acid MIDs, the second is the measured peptide MIDs.
See example input files in the data folder. An output file will be generated named "ex_*.csv". example code:
>> Main_fitting('..\data\1153aa-uncorrected.csv','..\data\xz1153peptide.csv','2H')

for the simulation of individual peptide MID, use the following example code: 
>> [out,tb] =peptide_mid('AKRG')

THe output stores the detailed MID information in structured array.




