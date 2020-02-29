## Code for the User Interface through Shiny

ui <- fluidPage(
  navbarPage("crispRdesignR",
             tabPanel("sgRNA Designer",
                      titlePanel("sgRNA Designer"),
                      sidebarLayout(
                        sidebarPanel(
                          textInput("sequence", "Target Sequence", placeholder = "Paste target DNA sequence here"),
                          checkboxInput("fasta", "Use FASTA or txt file as target sequence", value = FALSE),
                          tags$div(id = "placeholder1"),
                          selectInput("genome_select", "Select Genome",
                                      c("Homo sapiens (UCSC.hg19)" = "BSgenome.Hsapiens.UCSC.hg19",
                                        "Homo sapiens (UCSC.hg38)" = "BSgenome.Hsapiens.UCSC.hg38",
                                        "Saccharomyces cerevisiae (UCSC.sacCer2)" = "BSgenome.Scerevisiae.UCSC.sacCer2",
                                        "Mus musculus (UCSC.mm10)" = "BSgenome.Mmusculus.UCSC.mm10",
                                        "Drosphila melanogaster (UCSC.dm6)" = "BSgenome.Dmelanogaster.UCSC.dm6",
                                        "Pan troglodytes (UCSC.panTro5)" = "BSgenome.Ptroglodytes.UCSC.panTro5",
                                        "Escheria coli (NCBI.20080805)" = "BSgenome.Ecoli.NCBI.20080805",
                                        "Caenorhabditis elegans (UCSC.ce11)" = "BSgenome.Celegans.UCSC.ce11",
                                        "Ratus norvegicus (UCSC.rn6)" = "BSgenome.Rnorvegicus.UCSC.rn6",
                                        "Arabidopsis thaliana (TAIR.04232008)" = "BSgenome.Athaliana.TAIR.04232008")),
                          fileInput("gtf_file", "Choose genome annotation file (.gtf)",
                                    multiple = FALSE),
                          checkboxInput("options_toggle", "Additional Options", value = FALSE),
                          tags$div(id = "placeholder5"),
                          actionButton("run", "Find sgRNA", icon("paper-plane"))
                        ),
                        mainPanel(
                          tags$div(id = "placeholder3"),
                          DT::dataTableOutput("sgRNA_data"),
                          tags$div(id = "placeholder4"),
                          DT::dataTableOutput("offtarget_data"),
                          titlePanel("About"),
                          column(12, HTML("The Cas9 Guide Finder designs guide RNA sequences (sgRNA) for Cas9 DNA editing.
                                          To begin, enter a sequence into the sequence box, select a genome to search for
                                          Off-Targets, provide a genome annotation file (.gtf) specific to your genome, and click find sgRNA. <br/><br/> Note about Off-target calling in large genomes: When using a large genome like
                                          Homo sapiens, we reccomend using sequences under 250 base pairs. The time it can take
                                          to search these genomes can be multiple hours if too many sgRNA are generated."))
                          )
                          )
                          )
                          )
                        )

server <- function(input, output) {
  ## Increases the maximum file size that can be uploaded to Shiny to accomadate .gtf files
  options(shiny.maxRequestSize=150*1024^2)

  ## Creates a list of reactive values that allows the program to
  ## update only when the action button is pressed
  maindf <- reactiveValues(data = NULL)
  downloadmaindf <- reactiveValues(data = NULL)
  offtargetdf <- reactiveValues(data = NULL)
  downloadofftargetdf <- reactiveValues(data = NULL)

  ## Creates default values for the arguments in the find sgRNA function
  callofftargets <- "yes_off"
  annotateofftargets <- "yes_annotate"
  givenPAM <- "NGG"

  ## Creates a variable for the gene annotation file
  gtf_datapath <<- 0
  gene_annotation_file <<- 0

  ## Runs the sgRNA_design function when the action button is pressed
  observeEvent(input$run, {
    callofftargets <- input$'toggle_off_targets'
    annotateofftargets <- input$'toggle_off_annotation'
    if (input$options_toggle == TRUE) {
      if (nchar(paste(input$'customPAM')) <= 6) {
        givenPAM <- paste(input$'customPAM')
      } else {
        showModal(modalDialog(
          title = "Error",
          "Custom PAM must be 6 base pairs or less"
        ))
      }
    }
    ## Sets default parameteres for the sgRNA_design function
    if (is.null(callofftargets)){
      callofftargets <- "yes_off"
    }
    if (is.null(annotateofftargets)){
      annotateofftargets <- "yes_annotate"
    }
    if (input$'fasta' == TRUE) {
      if (isTRUE(grep("*.fasta", input$'fastafile'$datapath) == 1)) {
        sequence <- rtracklayer::import(input$'fastafile'$datapath, format = "fasta")
        sequence <- as.character(sequence)
      } else if (isTRUE(grep("*.txt", input$'fastafile'$datapath) == 1)) {
        sequence <- read.table(input$'fastafile'$datapath)
        sequence <- paste(sequence[1:nrow(sequence), 1], collapse = "")
      } else {
        showModal(modalDialog(
          title = "Error",
          "Unsuported Target Sequence File Type"
        ))
      }
    } else {
      sequence <- paste(input$'sequence', collapse = "")
      sequence <- stringr::str_replace_all(sequence, stringr::fixed(" "), "")
    }
    # Check to see if input is valid
    if (isTRUE(try(class(Biostrings::DNAString(sequence)) == "DNAString"))) {
      # Create a Progress object
      designprogress <- shiny::Progress$new()
      designprogress$set(message = "Preparing gene annotation file", value = 0, detail = "This may take a while")
      # Close the progress when this reactive exits (even if there's an error)
      on.exit(designprogress$close())
      if (callofftargets == "no_off" | annotateofftargets == "no_annotate") {
        annotating <- FALSE
      } else {
        annotating <- TRUE
      }
      if (annotating == FALSE) {
        gene_annotation_file <<- "placeholder"
      }
      if ((annotating == TRUE) & (is.null(input$'gtf_file'$datapath) == TRUE)) {
        showModal(modalDialog(
          title = "Error",
          "Please provide a genome annotation file (.gtf.gz)"
        ))
      } else {
        if ((annotating != FALSE) & (gtf_datapath == 0)) {
          gtf_datapath <<- input$'gtf_file'$datapath
          gene_annotation_file <<- rtracklayer::import.gff(input$'gtf_file'$datapath)
        }
        if (annotating != FALSE) {
          if (gtf_datapath != input$'gtf_file'$datapath) {
            gtf_datapath <<- input$'gtf_file'$datapath
            gene_annotation_file <<- rtracklater::import.gff(input$'gtf_file'$datapath)
          }
        }
        ## Inititates sgRNA_design_function
        all_data <- sgRNA_design_function(userseq = sequence, genomename = input$'genome_select', gtf = gene_annotation_file, userPAM = givenPAM, designprogress,
                                           calloffs = callofftargets, annotateoffs = annotateofftargets)
        if ((length(all_data) == 0) == FALSE) {
          ## Starts creating the sgRNA table
          int_sgRNA_data <- data.frame(all_data[1:15])
          ## Adds color to indicate unfavorable GC content
          GCinstance <- unlist(int_sgRNA_data[6])
          GCindex <- which(GCinstance >=.8 | GCinstance <=.3)
          GCinstance_color <- as.character(GCinstance)
          for (G in 1:length(GCinstance)){
            if (G %in% GCindex){
              GCinstance_color[G] <- paste('<span style="color:red">', GCinstance_color[G], '<span>', sep = "")
            }
          }
          ## Adds color to indicate unfavorable homopolymers
          Homopolymerdetect <- unlist(int_sgRNA_data[7])
          Homopolymerdetect_color <- as.character(Homopolymerdetect)
          for (H in 1:length(Homopolymerdetect)){
            if (Homopolymerdetect[H] == "TRUE"){
              Homopolymerdetect_color[H] <- paste('<span style="color:red">', Homopolymerdetect[H], '<span>', sep = "")
            }
          }
          proc_sgRNA_data <- data.frame(int_sgRNA_data[1:5], GCinstance_color, Homopolymerdetect_color, int_sgRNA_data[8:15])
          colnames(int_sgRNA_data) <- c("sgRNA sequence", "PAM sequence", "Direction", "Start", "End", "GC content",
                                        "Homopolymer", "Self Complementary", "Efficiency Score", "MM0", "MM1", "MM2", "MM3", "MM4", "Notes")
          colnames(proc_sgRNA_data) <- c("sgRNA sequence", "PAM sequence", "Direction", "Start", "End", "GC content",
                                         "Homopolymer", "Self Complementary", "Efficiency Score", "MM0", "MM1", "MM2", "MM3", "MM4", "Notes")
          if (input$run == 1) {
            insertUI(
              selector = "#placeholder3",
              where = "afterEnd",
              ui = tags$div(id = 'sgRNAdftext',
                            titlePanel("sgRNA Table"),
                            downloadButton("Download_sgRNA", "Download sgRNA")
              )
            )
          }
          ## Outputs the Table
          maindf$sgRNA_data <- proc_sgRNA_data
          downloadmaindf$sgRNA_data <- int_sgRNA_data
          ## Starts creating the off-target annotation table
          int_offtarget_data <- data.frame(all_data[16:27])
          ## Adds code to color mismatches red within the off target sequences
          off_offseq <- as.character(unlist(int_offtarget_data[8]))
          off_sgRNAseq <- as.character(unlist(int_offtarget_data[1]))
          for (x in 1:length(off_offseq)) {
            justsgRNA <- off_sgRNAseq[x]
            justoff <- off_offseq[x]
            splitjustsgRNA <- stringr::str_split(justsgRNA, "", simplify = TRUE)
            splitoffsgRNA <- stringr::str_split(justoff, "", simplify = TRUE)
            mismatches <- which(splitjustsgRNA != splitoffsgRNA)
            splitlistoffsgRNA <- as.list(splitoffsgRNA)
            if (length(mismatches) != 0){
              for (g in length(mismatches):1) {
                splitlistoffsgRNA <- append(splitlistoffsgRNA, '</span>', after = mismatches[g])
                splitlistoffsgRNA <- append(splitlistoffsgRNA, '<span style="color:red">', after = mismatches[g]-1)
              }
              off_offseq[x] <- paste(splitlistoffsgRNA, sep="", collapse = "")
            }
          }
          proc_offtarget_data <- data.frame(int_offtarget_data[1:7], off_offseq, int_offtarget_data[9:12])
          colnames(int_offtarget_data) <- c("sgRNA sequence", "Chromosome", "Start", "End", "Mismatches", "Direction", "CFD Scores",
                                            "Off-target sequence", "Gene ID", "Gene Name", "Sequence Type", "Exon Number")
          colnames(proc_offtarget_data) <- c("sgRNA sequence", "Chromosome", "Start", "End", "Mismatches", "Direction", "CFD Scores",
                                            "Off-target sequence", "Gene ID", "Gene Name", "Sequence Type", "Exon Number")
          if (input$run == 1) {
            insertUI(
              selector = "#placeholder4",
              where = "afterEnd",
              ui = tags$div(id = 'sgRNAofftext',
                            titlePanel("Off-target Information"),
                            column(12, "Note: this program may report sequences in the target region as potential off-target sequences"),
                            downloadButton("Download_off", "Download Off-Targets")
              )
            )
          }
          offtargetdf$data <- proc_offtarget_data
          downloadofftargetdf$data <- int_offtarget_data
        } else {
          showModal(modalDialog(
            title = "Error",
            "No sgRNA were generated from sequence"
          ))
        }
      }
    } else {
      showModal(modalDialog(
        title = "Error",
        "Sequence may contain unsupported characters"
      ))
    }
  })

  output$Download_sgRNA <- downloadHandler(
    filename = function(){"sgRNA.csv"},
    content = function(file) {
      write.csv(downloadmaindf$sgRNA_data, file, row.names = TRUE)
    }
  )

  output$Download_off <- downloadHandler(
    filename = function(){"Offtarget.csv"},
    content = function(file) {
      write.csv(downloadofftargetdf$data, file, row.names = TRUE)
    }
  )

  ## Reactively outputs an sgRNA table when the function is complete
  output$sgRNA_data <- DT::renderDataTable({maindf$sgRNA_data}, escape = FALSE)
  output$offtarget_data <- DT::renderDataTable({offtargetdf$data}, escape = FALSE)

  ## Add fasta file input to the UI
  observeEvent(input$fasta, {
    if (input$fasta == TRUE) {
      insertUI(
        selector = "#placeholder1",
        where = "afterEnd",
        ui = tags$div(id = 'fastainput',
                      fileInput("fastafile", "Choose fasta file",
                                multiple = FALSE)
        )
      )
    } else {
      removeUI(
        selector = 'div#fastainput',
        multiple = FALSE
      )
    }
  })

  ## Add Additional Options input to the UI
  observeEvent(input$options_toggle, {
    if (input$options_toggle == TRUE) {
      insertUI(
        selector = "#placeholder5",
        where = "afterEnd",
        ui = tags$div(id = 'optionsmenu',
                      column(12, HTML("Warning: Doench score not accurate for custom PAMS")),
                      textInput("customPAM", "Custom PAM (Max 6bp)", value = "NGG"),
                      selectInput("toggle_off_targets", "Call Off-Targets?",
                                  c("Yes" = "yes_off",
                                    "No" = "no_off"),
                                  selected = "yes_off"),
                      selectInput("toggle_off_annotation", "Annotate Off-Targets?",
                                  c("Yes" = "yes_annotate",
                                    "No" = "no_annotate"),
                                  selected = "yes_annotate")
        )
      )
    } else {
      removeUI(
        selector = 'div#optionsmenu',
        multiple = TRUE
      )
    }
  })

}

shinyApp(ui=ui, server=server)
