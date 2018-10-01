getsgRNAdata <- function(alldata){
  sgRNA_frame <- data.frame(alldata[1:15])
  colnames(sgRNA_frame) <- c("sgRNA sequence", "PAM sequence", "Direction", "Start", "End", "GC content",
                                      "TTTT Homopolymer", "Homopolymer", "Self Complementary", "Doench Score", "MM0", "MM1", "MM2", "MM3", "MM4")
  sgRNA_frame
}
