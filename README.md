# ROH_Calc

The ROH_Calc is a program that determines whether a region of the genome for an individual is in a ROH of a given length. The program reads in a parameter file that provides the location of necessary files along with quality control parameters. The program utilizes a sliding window approach and searches for stretches of continuous homozygous genotypes. The program output the proportion of the genome in an ROH, if a window or SNP is in a ROH across the genome for each individual.

## Parameter File
The parameter file contains all the user-defined variables. A list of them and some information about them is outlined below:

- **GENOTYPE_FILE:** The name of the genotype file and has to be a single word. The format is ID and followed by genotype string. The string of genotypes should be the same length across individuals. The delimiter is a space. The genotypes can’t be missing and they have to be in the format 0 (homozygote), 1 (heterozygote) and 2 (other homozygote). The row of the map should correspond the the location within the genotype string.
- **MAP_FILE:** The name of the map file and has to be a single word. The format of the map file is chromosome for column 1 and nucleotide position for column 2. The map file has to be ordered by chromosome and nucleotide position within chromosome. The delimiter is a space.
- **ROH_CUTOFF:** Refers to how long your ROH stretches have to be in nucleotides in order to be called a ROH stretch. For example a value of “5000000” refers to a cutoff of 5 Mb.

