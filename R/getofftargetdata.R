getofftargetdata <- function(alldata){
  off_targ_frame <- data.frame(alldata[15:26])
  colnames(off_targ_frame) <- c("sgRNA sequence", "Chromosome", "Start", "End", "Mismatches", "Direction", "CFD Scores",
                                "Off-target sequence", "Gene ID", "Gene Name", "Sequence Type", "Exon Number")
  off_targ_frame
}
