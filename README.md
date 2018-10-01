# crispRdesignR
Software used to design guide RNA sequences for CRISPR/Cas9 genome editing

This software aims to provide all scientifically pertinent information when designing guide RNA sequences for Cas9 genome editing. When provided a target DNA sequence for editing, a genome to check for off-targets in, and a genome annotation file (.gtf) to provide addition information about off-target matches it will out put information for two separate data tables. The first table contains all information on the generated sgRNA themselves (sgRNA sequence, PAM, Direction, Start, End, GC content, Presence of Homopolymers, Self Complementarity, Effciency Score (Doench 2014), and Genomic Matches). The second table contains all information on the found off-target sequences (Original sgRNA Sequence, Chromosome, Start, End, Number of Mismatches, Direction, CFD Scores, Matched Sequence, Gene ID, Gene Name, Sequence Type, and Exon Number). Additionally, a user may provide their own DNA libraries to search for off targets in and use a genome annotation file of their preference.

For a version of this software with a User interface through the Shiny package, see: https://github.com/dylanbeeber/Cas9-Guide-Designer.

# Requirements
This package requires three supplemental files:

Doench_Model_Weights_Singleonly.csv and Doench_Model_Weights_Doubleonly.csv - Two data tables used to assist with efficiency scoring. These must be put in the working directory when using the sgRNA_design function.

CFD_Scoring.csv - A data table that contains the information used to calculate the off-target effects of off-target sequences.

crispRdesignR dependencies:
The stringr package: `install.packages("stringr", repos='http://cran.us.r-project.org')`

The BioStrings and BSgenome packages through Bioconductor: `source("https://bioconductor.org/biocLite.R")`
`biocLite("Biostrings", "BSgenome", "BSgenome.Hsapiens.UCSC.hg19", "AnnotationHub")` The user may choose to download other genomes and use them with the crispRdesignR package.

# Usage
