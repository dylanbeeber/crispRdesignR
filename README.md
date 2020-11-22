# crispRdesignR
<!-- badges: start -->
[![Build Status](https://travis-ci.com/dylanbeeber/crispRdesignR.svg?branch=master)](https://travis-ci.com/dylanbeeber/crispRdesignR)
[![](https://www.r-pkg.org/badges/version/crispRdesignR?color=blue)](https://cran.r-project.org/package=crispRdesignR)
[![](https://img.shields.io/badge/doi-10.7150/jgen.41196-yellow.svg)](https://doi.org/10.7150/jgen.41196)
<!-- badges: end -->

Software used to design guide RNA sequences for CRISPR/Cas9 genome editing

This package aims to provide scientifically pertinent information when designing guide RNA sequences for Cas9 genome editing. When provided a target DNA sequence for editing, a genome to check for off-targets in, and a genome annotation file (.gtf), it will output information for two separate data tables. The first table contains information on the generated sgRNA (sgRNA sequence, PAM, Direction, Start, End, GC content, Presence of Homopolymers, Regions of Self Complementarity, Effciency Score (Doench 2016), Number of Potential Off-Target sequences, and Notes on the sgRNA). The second table contains information on the found off-target sequences (Original sgRNA Sequence, Chromosome, Start, End, Number of Mismatches, Direction, CFD Scores, Matched Sequence, Gene ID, Gene Name, Sequence Type, and Exon Number). This data can be generated through the command line or through crispRdesignR's GUI. Additionally, a user may provide their own DNA libraries and genome annotation files when searching for off-targets.

For more information, please see our article in the Journal of Genomics: http://www.jgenomics.com/v08p0062.htm

![crispRdesignRscreenshot4](https://user-images.githubusercontent.com/38253997/76813540-b40e9580-67ce-11ea-93c5-58f939ed6161.PNG)

## Installation and dependencies (tested in R version 3.6.2):

##### Install crispRdesignR v1.1.5 through CRAN

`install.packages("crispRdesignR")`

## Quick start with the GUI

Example Data is located in /inst/ folder.

The DAK1.fasta and DAK1_short.txt file contains a DNA sequence native to the DAK1 gene that can be copied and pasted into crispRdesignR or uploaded as a file (in the GUI version).

The "Saccharomyces_cerevisiae.R64-1-1.92.gtf.gz" file is the compressed genome annotation file for Saccharomyces cerevisiae. Both compressed and uncompressed .gtf files can be used.

Using the GUI version:

- `library(crispRdesignR)`

- `crispRdesignRUI()`

- Click on the “Use FASTA or txt file as target sequence” button and choose the DAK1.fasta or DAK1_short.txt file, or copy and paste the sequence in the box.

- select the Saccharomyces cerevisiae genome

- browse to choose the .gtf file Saccharomyces_cerevisiae.R64-1-1.92.gtf.gz

- click on the Find sgRNA button

![crispRdesignRscreenshot5](https://user-images.githubusercontent.com/38253997/76813686-20899480-67cf-11ea-89e3-447eef26d21a.PNG)

Additional Genome annotation files can be found here: https://useast.ensembl.org/info/data/ftp/index.html

Note: Even though it might be possible to select them in the GUI, Genomes must be installed (with `install.packages(BSgenome.yourgenome)`) before they can be used in the shiny app.

## Command-line version

- Load crispRdesignR with `library(crispRdesignR)` and ensure that your genome of interest is installed. The following commands can be used to output the same data tables as in the GUI version.

##### sgRNA Design

All data can be generated without the graphic interface by using a single function: `sgRNA_design(userseq, genomename, gtfname, userPAM, calloffs = TRUE, annotateoffs = TRUE)`

`userseq`: The target sequence to generate sgRNA guides for. Can either be a character sequence containing DNA bases or the name of a fasta/text file in the working directory.

`genomename`: The name of a genome (in BSgenome format) to check for off-targets in. These genomes can be downloaded through BSgenome or compiled by the user.

`gtfname`: The name of a genome annotation file (.gtf) in the working directory to check off-target sequences against.

`userPAM`: An optional argument used to set a custom PAM for the sgRNA. If not set, the function will default to the "NGG" PAM. Warning: Doench efficieny scores are only accurate for the "NGG" PAM.

`calloffs`: If TRUE, the function will search for off-targets in the genome chosen specified by the genomename argument. If FALSE, off-target calling will be skipped.

`annotateoffs`: If TRUE, the function will provide annotations for the off-targets called using the genome annotation file specified by the gtfname argument. If FALSE, off-target annotation will be skipped.

##### sgRNA Data Retrieval

The data on the generated sgRNA sequences can be retrieved with: `getsgRNAdata(x)`

`x`: The raw data generated by `sgRNA_design()`

##### Off-Target Data Retrieval

The additional off-target data can be retrieved with `getofftargetdata(x)`

`x`: The raw data generated by `sgRNA_design()`

###### Example:

`testseq <- "GGCAGAGCTTCGTATGTCGGCGATTCATCTCAAGTAGAAGATCCTGGTGCAGTAGG"
usergenome <- BSgenome.Scerevisiae.UCSC.sacCer2::BSgenome.Scerevisiae.UCSC.sacCer2
gtfname <- "Saccharomyces_cerevisiae.R64-1-1.92.gtf.gz"
annotation_file <- system.file("example_data", gtfname, package = "crispRdesignR")
alldata <- sgRNA_design(testseq, usergenome, annotation_file)`

`sgRNAdata <- getsgRNAdata(alldata)`

`offtargetdata <- getofftargetdata(alldata)`

###### Example:
`exampledata <- sgRNA_design("DAK1.fasta", BSgenome.Scerevisiae.UCSC.sacCer2, "Saccharomyces_cerevisiae.R64-1-1.92.gtf.gz", "NAG", calloffs = TRUE, annotateoffs = FALSE)`

## Example Data

Example Data is located in /inst/ folder. To use the DAK1.fasta file, place it in the working directory and refer to it in crispRdesignR. The DAK1_short.txt file contains a short DNA sequence that can be copied and pasted into crispRdesignR. Both sequences are native to the DAK1 gene in Saccharomyces cerevisiae. The "Saccharomyces_cerevisiae.R64-1-1.92.gtf.gz" file is a genome annotation file for Saccharomyces cerevisiae and must also be placed in the working directory (when using the command line version).
