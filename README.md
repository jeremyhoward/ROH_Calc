# ROH_Calc

C++ code that I have utilized to determine whether a region of the genome for an individual is in a ROH of a given length. The program reads in a parameter file that provides the location of necessary files along with quality control parameters. The program utilizes a sliding window approach and searches for stretches of continuous homozygous genotypes. The program output the proportion of the genome in an ROH, if a window or SNP is in a ROH across the genome for each individual. An updated version that generated errors if input parameters are not formated correctly was uploaded on May 5, 2018.

## To Compile and Run
Compile: g++ ROH_CAL.cpp -o ROH_CALC

Run: ./ROH_CALC param.txt

## Parameter File
The parameter file contains all the user-defined variables. The program reads the parameter file by searching for keywords (i.e. in bold below) that are capitalized and then followed by a colon. Therefore any phrase that does not meet the search criteria is ignored when initializing parameters within the program.  All parameters are required for the program to run. A list of them and some information about them is outlined below:

- **GENOTYPE_FILE:** The name of the genotype file and has to be a single word. The format is ID and followed by genotype string. The string of genotypes should be the same length across individuals. The delimiter is a space. The genotypes can’t be missing and they have to be in the format 0 (homozygote), 1 (heterozygote) and 2 (other homozygote). The row of the map should correspond the the location within the genotype string.
- **MAP_FILE:** The name of the map file and has to be a single word. The format of the map file is chromosome for column 1 and nucleotide position for column 2. The map file has to be ordered by chromosome and nucleotide position within chromosome. The delimiter is a space.
- **ROH_CUTOFF:** Refers to how long your ROH stretches have to be in Megabases in order to be called a ROH stretch. For example a value of “5” refers to a cutoff of 5 Mb.
- **ROH_THRESHOLD:** This threshold is given to remove regions that have a small number of SNP and therefore could give rise to spurious results. The value refers to how many standard deviations from the mean number of SNP that are removed. For instance a value of 2 refers to any window that contains less than 2 * Average SNP will be removed.
- **REMOVE_SNP:** Refers to whether you want the SNP that weren’t contained in a long enough ROH window to be removed. The two options are 'yes' or 'no'.
- **OUT_FILE:** Name of output files.

## Output Files
- Autozygosity: Autozygosity file. The first row of the file contains the ID and then the location of each SNP (chromosome_position). Any row after the file is then an indicator of whether that SNP is in an ROH. Not in an ROH is 0 and in an ROH is 1.
- ROH: ROH Window data. The first row of the file contains the ID and then location of each window (chromosome_ start position_end position). Any row after the file is then an indicator of whether that ROH window is a ROH or not. Not in an ROH is 0 and in an ROH is 1.
- Summary: The first column is the animalID, second is homozygosity and third is proportion of genome in ROH based on cutoff you chose in parameter file.
