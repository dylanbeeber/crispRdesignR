## This function needs a genome (from BSgenome to check for off-targets in and
## a gtf file to annotate the offtargets)

sgRNA_design <- function(userseq, genomename, gtfname, calloffs = TRUE, annotateoffs = TRUE){
  if (isTRUE(str_detect(userseq, ".fasta")) || (isTRUE(str_detect(userseq, ".txt")))) {
    if (isTRUE(class(try(import(userseq, format = "fasta"))) == "DNAStringSet")) {
      Biostrings_sequence <- import(userseq, format = "fasta")
      sequence <- as.character(Biostrings_sequence)
    } else {
      sequence <- read.table(userseq)
      sequence <- paste(sequence[1:nrow(sequence), 1], collapse = "")
      Biostrings_sequence <- DNAString(sequence)
    }
  } else {
    sequence <- paste(userseq, collapse = "")
    sequence <- str_replace_all(sequence, fixed(" "), "")
    Biostrings_sequence <- DNAString(sequence)
  }
  ## Creates a character string that contains the
  ## complementary sequence (Both in the reverse
  ## direction and with substituted nucleotides)
  Biostrings_rev_seq <- reverseComplement(Biostrings_sequence)
  rev_seq <- as.character(Biostrings_rev_seq)
  ## Create an empty list for the forward sgRNA to go
  sgRNA_list_f <- c()
  ## Create an empty list for the reverse sgRNA to go
  sgRNA_list_r <- c()
  ## Create four empty lists for forward and reverse start and end positions
  ## Check to make the start and end numbers are correct
  sgRNA_f_start <- c()
  sgRNA_f_end <- c()
  sgRNA_r_start <- c()
  sgRNA_r_end <- c()
  ## Sets number that helps determine when to stop looking
  ## for possible sgRNA (This prevents it from choosing
  ## an incomplete sgRNA at the very end of the sequence)
  num_char_in_seq <- nchar(sequence) - 29
  ## Sets the PAM sequence and determines what the program
  ## will describe as a possible sgRNA. Even though most sgRNA is
  ## only 20 nucleotides long, nucleotides surrounding the sgRNA
  ## are used for study-based scoring
  PAM <- ".........................GG..."
  ## Sets n to zero, which is used to incrementally increase
  ## and search for all possible sgRNA in the sequence
  n <- 0
  ## Searches all 23 nt streches in the sequence for
  ## possible matches to the PAM, then puts entire 30 nt
  ## matches into a list (including the PAM)
  for (x in 0:num_char_in_seq){
    poss_sgRNA <- substr(sequence, 1+n, 30+n)
    if (str_detect(poss_sgRNA, PAM) == TRUE){
      sgRNA_list_f[[length(sgRNA_list_f)+1]] <- poss_sgRNA
      sgRNA_f_start[[length(sgRNA_f_start)+1]] <- n+5
      sgRNA_f_end[[length(sgRNA_f_end)+1]] <- n+27
    }
    n <- n+1
  }
  ## Same as above but with the reverse sequence
  n <- 0
  for (x in 0:num_char_in_seq){
    poss_sgRNA <- substr(rev_seq, 1+n, 30+n)
    if (str_detect(poss_sgRNA, PAM) == TRUE){
      sgRNA_list_r[[length(sgRNA_list_r)+1]] <- poss_sgRNA
      sgRNA_r_start[[length(sgRNA_r_start)+1]] <- nchar(rev_seq)-n+5
      sgRNA_r_end[[length(sgRNA_r_end)+1]] <- nchar(rev_seq)-n+27
    }
    n <- n+1
  }
  ## Removes any sgRNA that contain degerate bases
  sgRNA_list_f <- sgRNA_list_f[grepl("[UWSMKRYBDHVNZ]", sgRNA_list_f) == FALSE]
  sgRNA_f_start <- sgRNA_f_start[grepl("[UWSMKRYBDHVNZ]", sgRNA_list_f) == FALSE]
  sgRNA_f_end <- sgRNA_f_end[grepl("[UWSMKRYBDHVNZ]", sgRNA_list_f) == FALSE]
  sgRNA_list_r <- sgRNA_list_r[grepl("[UWSMKRYBDHVNZ]", sgRNA_list_r) == FALSE]
  sgRNA_r_start <- sgRNA_r_start[grepl("[UWSMKRYBDHVNZ]", sgRNA_list_r) == FALSE]
  sgRNA_r_end <- sgRNA_r_end[grepl("[UWSMKRYBDHVNZ]", sgRNA_list_r) == FALSE]
  ## Creates list for all sgRNA to go
  sgRNA_list <- c(sgRNA_list_f, sgRNA_list_r)
  if (is.null(sgRNA_list) == FALSE) {
    sgRNA_start <- c(sgRNA_f_start, sgRNA_r_start)
    sgRNA_end <- c(sgRNA_f_end, sgRNA_r_end)
    ## Creates a list with only the sgRNA sequences (no PAMs or
    ## flanking sequence)
    breakseq <- function(seqlist){
      str_sub(seqlist, 5, 24)
    }
    sgRNA_seq <- sapply(sgRNA_list, breakseq)
    ## Creates a list with only the PAM sequences
    breakPAM <- function(seqlist){
      str_sub(seqlist, 25, 27)
    }
    sgRNA_PAM <- sapply(sgRNA_list, breakPAM)
    ## Makes a list of all of the sgRNA sequences with their PAM
    sgRNA_with_PAM <- paste(sgRNA_seq, sgRNA_PAM, sep = "")
    ## Makes a list of whether sgRNA are forward or reverse
    sgRNA_fow <- rep("+", each = length(sgRNA_list_f))
    sgRNA_rev <- rep("-", each = length(sgRNA_list_r))
    sgRNA_fow_or_rev <- c(sgRNA_fow, sgRNA_rev)
    ## Find GC percentage for each sgRNA and puts that data into
    ## a list called "GCinstance"
    FindGC <- function(seqlist){
      ((str_count(seqlist, "G") + str_count(seqlist, "C")) / 20)
    }
    GCinstance <- sapply(sgRNA_seq, FindGC)
    ## Creates a list that determines if the GC percentage
    ## is within the user-defined (WIP) threshold
    EvalGC <- function(GC){
      isTRUE(30<=GC&GC<= 70)
    }
    GClist <- sapply(GCinstance, EvalGC)
    ## Find homopolmers
    Findhomopolymer <- function(seqlist){
      str_detect(seqlist, "TTTT|AAAA|GGGG|CCCC")
    }
    Homopolymerdetect <- sapply(sgRNA_seq, Findhomopolymer)
    Homopolymerdetect
    ## Find TTTT homopolymers
    FindTTTThomopolymer <- function(seqlist){
      str_detect(seqlist, "TTTT")
    }
    TTTTHomopolymerdetect <- sapply(sgRNA_seq, FindTTTThomopolymer)
    ## Detect Self complementarity
    self_comp_list <- c()
    backbone_area <- DNAString("AGGCTAGTCCGT")
    revcomp_backbone <- reverseComplement(backbone_area)
    SpeFindGC <- function(seqlist){
      ((str_count(seqlist, "G") + str_count(seqlist, "C")) / 4)
    }
    for (j in 1:length(sgRNA_seq)){
      testDNA <- DNAString(sgRNA_seq[j])
      revcompDNA <- reverseComplement(testDNA)
      individ_comp_list <- c()
      for (r in 1:17) {
        test_region <- substr(testDNA, r, r+3)
        if (SpeFindGC(test_region) >= .5) {
          if (r <= 10) {
            compcount <- countPattern(test_region, reverseComplement(DNAString(substr(testDNA, r+7, length(testDNA)))))
            individ_comp_list[[length(individ_comp_list)+1]]  <- compcount
          }
          compcount <- countPattern(test_region, revcomp_backbone)
          individ_comp_list[[length(individ_comp_list)+1]]  <- compcount
        } else {
          individ_comp_list[[length(individ_comp_list)+1]] <- 0
        }
      }
      self_comp_list[[length(self_comp_list)+1]] <- sum(individ_comp_list)
    } ## Self comp checking ends here
    ## Assign a study-based efficiency score
    ## Following two lines retrieve the penalty constants (one for single nucleotides, the other for paired nucleotides)
    Doench_model_weights_singleonly <- read.csv("Doench_Model_Weights_Singleonly.csv", header = FALSE)
    Doench_model_weights_doubleonly <- read.csv("Doench_Model_Weights_Doubleonly.csv", header = FALSE)
    ## Creates an empty list for Doench study-based scores to go into
    Doench_Score <- c()
    ## "g" allows this script to go through the sgRNA_list, one by one
    g <- 0
    for (h in 1:length(sgRNA_list)){
      ## Splits the sgRNA into individual nucleotides
      split_sgRNA <- str_split(sgRNA_list[1+g], "", simplify = TRUE)
      ## Creates the sgRNA_model_weight list and adds the intercept to it
      sgRNA_model_weights <- c(0.597636154)
      ## Finds the model weights for each individual nucleotide in the sgRNA and adds them to sgRNA_model_weights
      n <- 0
      for (t in 1:39){
        if ((split_sgRNA[,Doench_model_weights_singleonly[1+n, 2]] == Doench_model_weights_singleonly[1+n, 1]) == TRUE){
          sgRNA_model_weights[[length(sgRNA_model_weights)+1]] <- Doench_model_weights_singleonly[1+n, 3]
        }
        n <- n+1
      }
      ## Splits the sgRNA into pieces two nucleotides long
      ## Serves the same purpose as: "split_sgRNA <- str_split(sgRNA_list[1+g], "", simplify = TRUE)"
      double <- ".."
      split_double_sgRNA <- c()
      n <- 0
      for (xp in 1:29){
        double_sgRNA <- substr(sgRNA_list[1+g], 1+n, 2+n)
        split_double_sgRNA[[length(split_double_sgRNA)+1]] <- double_sgRNA
        n <- n+1
      }
      ## Finds the model weights for each double nucleotide in the sgRNA and adds them to sgRNA_model_weights
      n <- 0
      for (t in 1:31){
        if ((split_double_sgRNA[Doench_model_weights_doubleonly[1+n, 2]] == Doench_model_weights_doubleonly[1+n, 1]) == TRUE){
          sgRNA_model_weights[[length(sgRNA_model_weights)+1]] <- Doench_model_weights_doubleonly[1+n, 3]
        }
        n <- n+1
      }
      GC <- (((str_count(sgRNA_seq[1+g], "G") + str_count(sgRNA_seq[1+g], "C")) / 20) * 100)
      if (isTRUE(GC>70 & 80>=GC) == TRUE) {
        sgRNA_model_weights[[length(sgRNA_model_weights)+1]] <- -0.202625894258
      }
      if (isTRUE(GC<30 & GC>=20) == TRUE) {
        sgRNA_model_weights[[length(sgRNA_model_weights)+1]] <- -0.202625894258
      }
      if (isTRUE(GC<20 | GC>80) == TRUE) {
        sgRNA_model_weights[[length(sgRNA_model_weights)+1]] <- -0.166587751983
      }
      ## Completes the calculation described in Doench et al. (2014) and adds it to the Doench Score list
      Doench_Score[[length(Doench_Score)+1]] <- round(1/(1+exp(-(sum(sgRNA_model_weights)))), digits = 3)
      g <- g+1
    }
    if (calloffs == FALSE) {
      mm0_list <- rep("NA", each = length(sgRNA_list))
      mm1_list <- rep("NA", each = length(sgRNA_list))
      mm2_list <- rep("NA", each = length(sgRNA_list))
      mm3_list <- rep("NA", each = length(sgRNA_list))
      mm4_list <- rep("NA", each = length(sgRNA_list))
      ## Creates data table with all available sgRNA data
      sgRNA_data <- data.frame(sgRNA_seq, sgRNA_PAM, sgRNA_fow_or_rev, sgRNA_start, sgRNA_end, GCinstance, TTTTHomopolymerdetect, Homopolymerdetect, self_comp_list, Doench_Score, mm0_list, mm1_list, mm2_list, mm3_list, mm4_list)
      ## Set the names of each column
      colnames(sgRNA_data) <- c("sgRNA sequence", "PAM sequence", "Direction", "Start", "End", "GC content", "TTTT Homopolymer", "Homopolymer", "Self Complementary", "Doench Score", "MM0", "MM1", "MM2", "MM3", "MM4")
      sgRNA_data <- sgRNA_data[order(-sgRNA_data$`Doench Score`),]
      ## Creates an empty data table for off-target annotation
      all_offtarget_info <- data.frame("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA")
      colnames(all_offtarget_info) <- c("sgRNA sequence", "Chromosome", "Start", "End", "Mismatches", "Direction", "CFD Score", "Off-target sequence", "Gene ID", "Gene Name", "Sequence Type", "Exon Number")
      data_list <- c("sgRNA_data" = sgRNA_data, "all_offtarget_info" = all_offtarget_info)
      data_list
    } else {
      ## Check for off-targets in the genome
      ## Creates Function that converts all sgRNAs into a format readable
      ## by Biostrings
      multiple_DNAString <- function(seqlist){
        DNAString(seqlist)
      }
      Biostrings_sgRNA <- lapply(sgRNA_with_PAM, multiple_DNAString)
      Off_targ_Bio_sgRNA <- Biostrings_sgRNA
      ## Define genome
      usegenome <- genomename
      seqnames <- seqnames(usegenome)
      ## Creates a series of lists to store the incoming mismatch information
      mm0_list <- c()
      mm1_list <- c()
      mm2_list <- c()
      mm3_list <- c()
      mm4_list <- c()
      off_start <- c()
      off_end <- c()
      off_direction <- c()
      off_sgRNAseq <- c()
      off_offseq <- c()
      off_chr <- c()
      off_mismatch <- c()
      PAM_test_list <- c("GG", "AG", "CG", "GA", "GC", "GT", "TG")
      for (seqname in seqnames) {
        print(paste("Checking for Off-Targets in", seqname, sep = " "))
        chrmm0_list <- c()
        chrmm1_list <- c()
        chrmm2_list <- c()
        chrmm3_list <- c()
        chrmm4_list <- c()
        revchrmm0_list <- c()
        revchrmm1_list <- c()
        revchrmm2_list <- c()
        revchrmm3_list <- c()
        revchrmm4_list <- c()
        for (pattern in Off_targ_Bio_sgRNA) {
          usepattern <- replaceLetterAt(pattern, 21, "N")
          subject <- usegenome[[seqname]]
          off_info <- matchPattern(usepattern, subject, max.mismatch = 4, min.mismatch = 0, fixed = FALSE)
          mis_info <- mismatch(usepattern, off_info, fixed = FALSE)
          rev_pattern <- reverseComplement(usepattern)
          rev_off_info <- matchPattern(rev_pattern, subject, max.mismatch = 4, min.mismatch = 0, fixed = FALSE)
          rev_mis_info <- mismatch(rev_pattern, rev_off_info, fixed = FALSE)
          if (length(off_info) > 0) {
            for (f in 1:length(off_info)) {
              if (substr(as.character(off_info[[f]]), 22, 23) %in% PAM_test_list) {
                off_start[[length(off_start)+1]] <- start(off_info)[f]
                off_end[[length(off_end)+1]] <- end(off_info)[f]
                off_direction[[length(off_direction)+1]] <- "+"
                off_chr[[length(off_chr)+1]] <- seqname
                off_mismatch[[length(off_mismatch)+1]] <- length(mis_info[[f]])
                off_sgRNAseq[[length(off_sgRNAseq)+1]] <- as.character(pattern)
                off_offseq[[length(off_offseq)+1]] <- as.character(off_info[[f]])
              }
            }
          }
          if (length(rev_off_info) > 0) {
            for (f in 1:length(rev_off_info)) {
              if (substr(as.character(rev_off_info[[f]]), 22, 23) %in% PAM_test_list) {
                off_start[[length(off_start)+1]] <- start(rev_off_info)[f]
                off_end[[length(off_end)+1]] <- end(rev_off_info)[f]
                off_direction[[length(off_direction)+1]] <- "-"
                off_chr[[length(off_chr)+1]] <- seqname
                off_mismatch[[length(off_mismatch)+1]] <- length(rev_mis_info[[f]])
                off_sgRNAseq[[length(off_sgRNAseq)+1]] <- as.character(pattern)
                off_offseq[[length(off_offseq)+1]] <- as.character(rev_off_info[[f]])
              }
            }
          }
          individMM <- c()
          if (length(mis_info) > 0) {
            for (f in 1:length(mis_info)) {
              individMM[[length(individMM)+1]] <- length(mis_info[[f]])
            }
            chrmm0_list[[length(chrmm0_list)+1]] <- sum(individMM == 0)
            chrmm1_list[[length(chrmm1_list)+1]] <- sum(individMM == 1)
            chrmm2_list[[length(chrmm2_list)+1]] <- sum(individMM == 2)
            chrmm3_list[[length(chrmm3_list)+1]] <- sum(individMM == 3)
            chrmm4_list[[length(chrmm4_list)+1]] <- sum(individMM == 4)
          } else {
            chrmm0_list[[length(chrmm0_list)+1]] <- 0
            chrmm1_list[[length(chrmm1_list)+1]] <- 0
            chrmm2_list[[length(chrmm2_list)+1]] <- 0
            chrmm3_list[[length(chrmm3_list)+1]] <- 0
            chrmm4_list[[length(chrmm4_list)+1]] <- 0
          }
          individMM <- c()
          if (length(rev_mis_info) > 0) {
            for (f in 1:length(rev_mis_info)) {
              individMM[[length(individMM)+1]] <- length(rev_mis_info[[f]])
            }
            revchrmm0_list[[length(revchrmm0_list)+1]] <- sum(individMM == 0)
            revchrmm1_list[[length(revchrmm1_list)+1]] <- sum(individMM == 1)
            revchrmm2_list[[length(revchrmm2_list)+1]] <- sum(individMM == 2)
            revchrmm3_list[[length(revchrmm3_list)+1]] <- sum(individMM == 3)
            revchrmm4_list[[length(revchrmm4_list)+1]] <- sum(individMM == 4)
          } else {
            revchrmm0_list[[length(revchrmm0_list)+1]] <- 0
            revchrmm1_list[[length(revchrmm1_list)+1]] <- 0
            revchrmm2_list[[length(revchrmm2_list)+1]] <- 0
            revchrmm3_list[[length(revchrmm3_list)+1]] <- 0
            revchrmm4_list[[length(revchrmm4_list)+1]] <- 0
          }
        }
        if (is.null(mm0_list)) {
          mm0_list <- chrmm0_list + revchrmm0_list
          mm1_list <- chrmm1_list + revchrmm1_list
          mm2_list <- chrmm2_list + revchrmm2_list
          mm3_list <- chrmm3_list + revchrmm3_list
          mm4_list <- chrmm4_list + revchrmm4_list
        } else {
          mm0_list <- chrmm0_list + mm0_list + revchrmm0_list
          mm1_list <- chrmm1_list + mm1_list + revchrmm1_list
          mm2_list <- chrmm2_list + mm2_list + revchrmm2_list
          mm3_list <- chrmm3_list + mm3_list + revchrmm3_list
          mm4_list <- chrmm4_list + mm4_list + revchrmm4_list
        }
      }
      ## Calculates off-target scores for each off target sequence
      print("annotating off-targets")
      CFD_Model_Scores <- read.csv("CFD_Scoring.csv")
      off_model_PAMs <- c("AG", "CG", "GA", "GC", "GT", "TG")
      CFD_PAM_Scores <- data.frame(off_model_PAMs, c(0.259259, 0.107142, 0.069444, 0.022222, 0.016129, 0.038961))
      CFD_Scores <- c()
      for (x in 1:length(off_offseq)) {
        if (off_direction[x] == "-") {
          temporary_off <- DNAString(off_offseq[x])
          temporary_off <- reverseComplement(temporary_off)
          CFDoffsplit <- str_split(temporary_off, "", simplify = TRUE)
        } else {
          CFDoffsplit <- str_split(off_offseq[x], "", simplify = TRUE)
        }
        CFDsgRNAsplit <- str_split(off_sgRNAseq[x], "", simplify = TRUE)
        individ_scores <- c()
        for (g in 1:20) {
          if (CFDsgRNAsplit[g] != CFDoffsplit[g]) {
            index <- which(CFD_Model_Scores$Position==g & CFD_Model_Scores$sgRNA==CFDsgRNAsplit[g] & CFD_Model_Scores$DNA==CFDoffsplit[g])
            individ_scores[[length(individ_scores)+1]] <- CFD_Model_Scores[index,4]
          }
        }
        specific_PAM <- (paste(CFDoffsplit[22], CFDoffsplit[23], sep = ""))
        if (isTRUE(specific_PAM != "GG")){
          if (specific_PAM %in% off_model_PAMs) {
            PAM_index <- which(off_model_PAMs==specific_PAM)
            individ_scores[[length(individ_scores)+1]] <- CFD_PAM_Scores[PAM_index,2]
          } else {
            individ_scores[[length(individ_scores)+1]] <- 0
          }
        }
        if (length(individ_scores) == 0) {
          CFD_Scores[[length(CFD_Scores)+1]] <- 1
        }
        if (length(individ_scores) == 1) {
          CFD_Scores[[length(CFD_Scores)+1]] <- individ_scores[1]
        }
        if (length(individ_scores) == 2) {
          CFD_Scores[[length(CFD_Scores)+1]] <- prod(individ_scores[1], individ_scores[2])
        }
        if (length(individ_scores) == 3) {
          CFD_Scores[[length(CFD_Scores)+1]] <- prod(individ_scores[1], individ_scores[2], individ_scores[3])
        }
        if (length(individ_scores) == 4) {
          CFD_Scores[[length(CFD_Scores)+1]] <- prod(individ_scores[1], individ_scores[2], individ_scores[3], individ_scores[4])
        }
      }
      CFD_Scores <- round(CFD_Scores, digits = 3)
      ## Decides whether to annotate off_targets
      if (((sum(mm0_list) + sum(mm1_list) + sum(mm2_list) + sum(mm3_list)) == 0) || (annotateoffs == FALSE)) {
        ## Put lists in data frame
        sgRNA_data <- data.frame(sgRNA_seq, sgRNA_PAM, sgRNA_fow_or_rev, sgRNA_start, sgRNA_end, GCinstance, TTTTHomopolymerdetect, Homopolymerdetect, self_comp_list, Doench_Score, mm0_list, mm1_list, mm2_list, mm3_list, mm4_list)
        ## Set the names of each column
        colnames(sgRNA_data) <- c("sgRNA sequence", "PAM sequence", "Direction", "Start", "End", "GC content", "TTTT Homopolymer", "Homopolymer", "Self Complementary", "Doench Score", "MM0", "MM1", "MM2", "MM3", "MM4")
        sgRNA_data <- sgRNA_data[order(-sgRNA_data$`Doench Score`),]
        sgRNA_data
      } else {
        ## Creates a function that annotates the off-targets called above
        annotate_genome <- function(ochr, ostart, oend, odir, gtfname) {
          gtf <- import(gtfname)
          seqlevelsStyle(gtf) <- "UCSC"
          seqer <- unlist(ochr)
          starter <- as.numeric(ostart)
          ender <- as.numeric(unlist(oend))
          strander <- unlist(odir)
          off_ranges <- GRanges(seqer, IRanges(starter, ender), strander)
          olaps <- findOverlaps(off_ranges, gtf)
          geneid <- c()
          geneidlist <- c()
          genename <- c()
          genenamelist <- c()
          sequencetype <- c()
          sequencetypelist <- c()
          exonnumber <- c()
          exonnumberlist <- c()
          mcols(off_ranges)$gene_id <- c()
          for (p in 1:length(off_ranges)) {
            if (p %in% queryHits(olaps)) {
              geneid <- mcols(gtf)$gene_id[subjectHits(olaps[which(p == queryHits(olaps))])]
              geneid <- unique(geneid)
              geneidlist[[length(geneidlist)+1]] <- paste(geneid, collapse = ", ")
              genename <- mcols(gtf)$gene_name[subjectHits(olaps[which(p == queryHits(olaps))])]
              genename <- unique(genename)
              genenamelist[[length(genenamelist)+1]] <- paste(genename, collapse = ", ")
              sequencetype <- mcols(gtf)$type[subjectHits(olaps[which(p == queryHits(olaps))])]
              sequencetype <- unique(sequencetype)
              sequencetypelist[[length(sequencetypelist)+1]] <- paste(sequencetype, collapse = ", ")
              exonnumber <- mcols(gtf)$exon_number[subjectHits(olaps[which(p == queryHits(olaps))])]
              exonnumber <- unique(exonnumber)
              exonnumberlist[[length(exonnumberlist)+1]] <- paste(exonnumber, collapse = ", ")
            } else {
              geneidlist[[length(geneidlist)+1]] <- "NA"
              genenamelist[[length(genenamelist)+1]] <- "NA"
              sequencetypelist[[length(sequencetypelist)+1]] <- "NA"
              exonnumberlist[[length(exonnumberlist)+1]] <- "NA"
            }
          }
          mcols(off_ranges)$gene_id <- geneidlist
          more_off_info <- data.frame(geneidlist, genenamelist, sequencetypelist, exonnumberlist)
          more_off_info
        }
        ## Compiles data frame of all off-target annotations
        more_off_info <- annotate_genome(off_chr, off_start, off_end, off_direction, gtfname)
        ## Complies all extra sgRNA info into a separate data frame
        all_offtarget_info <- data.frame(off_sgRNAseq, off_chr, off_start, off_end, off_mismatch, off_direction, CFD_Scores, off_offseq, more_off_info$geneidlist, more_off_info$genenamelist, more_off_info$sequencetypelist, more_off_info$exonnumberlist)
        colnames(all_offtarget_info) <- c("sgRNA sequence", "Chromosome", "Start", "End", "Mismatches", "Direction", "CFD Scores", "Off-target sequence", "Gene ID", "Gene Name", "Sequence Type", "Exon Number")
        ## Put lists in data frame
        sgRNA_data <- data.frame(sgRNA_seq, sgRNA_PAM, sgRNA_fow_or_rev, sgRNA_start, sgRNA_end, GCinstance, TTTTHomopolymerdetect, Homopolymerdetect, self_comp_list, Doench_Score, mm0_list, mm1_list, mm2_list, mm3_list, mm4_list)
        ## Set the names of each column
        colnames(sgRNA_data) <- c("sgRNA sequence", "PAM sequence", "Direction", "Start", "End", "GC content", "TTTT Homopolymer", "Homopolymer", "Self Complementary", "Doench Score", "MM0", "MM1", "MM2", "MM3", "MM4")
        sgRNA_data <- sgRNA_data[order(-sgRNA_data$`Doench Score`),]
        data_list <- c("sgRNA_data" = sgRNA_data, "all_offtarget_info" = all_offtarget_info)
        data_list
      }
    }
  } else {
    data_list <- data.frame()
  }
}
