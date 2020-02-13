

# Program 1 - makeOrg() - directory connections -----------------------------------
#  Goal is to create the output directories I need
#  If the directories are already there, do nothing
makeOrg <- function(home){
  # Create the output directories
  dirs <- list.dirs(home, recursive = T)
  baf <- file.path(home,"BAF")
  out <- file.path(home,"Processed files")
  src <- file.path(home,"Output Tables")
  reports <- file.path(home, "Summary Reports")
  bafFiles <- file.path(baf,"Files")
  bafFigures <- file.path(baf,"Figures")

  if (!baf %in% dirs) dir.create(baf)
  if (!out %in% dirs) dir.create(out)
  if (!reports %in% dirs) dir.create(reports)
  if (!bafFiles %in% dirs) dir.create(bafFiles)
  if (!bafFigures %in% dirs) dir.create(bafFigures)

  rm(dirs,home)
}


# Program 2 - GSProcess() Process the GS output files ---------------------------------
# Three types of output files - each batch has all three
#  1.  Sample summary file - this is small, usually 24-96 samples per file
#  2.  FullData table - This is quite large, this is where i get the BAF stuff and some genotyping scores per sample
#  3.  SNP table -  stats are for each SNP across the entire batch
GSProcess <- function(home){

  # Step 1 - generate a list of files I must process ------------------------

  fileSrc <- file.path(home, "Output tables") # raw output files
  outputFolder <- file.path(home,"Processed Files") # already processed

  # Pull the files
  step1a <- lapply(c(fileSrc, outputFolder), function(x){

    # 1.  Generate a list of files
    foo <- data.frame(Files = list.files(x), stringsAsFactors=F)

    # Find Batch number
    foo$Batch <- with(foo, ifelse(
      toupper(substr(Files,1,8)) == "FULLDATA", sub("FullDataTable_", "", foo$Files, ignore.case = T), ifelse(
        toupper(substr(Files,1,7)) == "SAMPLES", sub("Samples_","",Files, ignore.case = T),
        sub("SNPTable_","",Files, ignore.case = T))))%>%
      sub(".txt","",.)

    # 2.  Get file paths
    foo$Path <- lapply(foo$Files, function(y){
      list.files(file.path(x), pattern = y, full.names=T)
    }) %>% as.character()

    # 3.  Get Time last edited
    foo$Time <- file.mtime(file.path(foo$Path))

    # 4.  Split into a list based on the file type
    foo$Type <- with(foo, ifelse(
      toupper(substr(Files,1,8)) == "FULLDATA",1, ifelse(
        toupper(substr(Files,1,7)) == "SAMPLES", 2, 3))) %>%
      factor(1:3, c("FullData","Samples","SNP"))

    # 5.  Split the files into sections
    files <- split(foo,foo$Type)
    return(files)
  })

  names(step1a) <- c("Original","Processed")

  # Figure out what needs to be processed by comparing time stamps
  step1b <- lapply(names(step1a$Original), function(x){
    one <- step1a$Original[[x]]
    two <- step1a$Processed[[x]]

    one$MergeVar <- sub(".txt","",one$Files)
    two$MergeVar <- sub(".RDS","",two$Files)

    joined <- dplyr::full_join(
                        dplyr::select(one, MergeVar, Files, Path, Type, origTime = Time, Batch),
                        dplyr::select(two, MergeVar, newTime = Time),
                        "MergeVar")
    joined$Convert <- with(joined, ifelse(
      is.na(newTime), "Process", ifelse(
        newTime > origTime, "Complete", "Process")))
    table(joined$Convert)

    convert <- dplyr::filter(joined, Convert == "Process") %>%
      dplyr::select(Files, Path, Type, Batch)
    return(convert)
  }) %>% do.call("rbind",.)

  # Split up by type
  step1b <- split(step1b, step1b$Type)
  sampleFiles <- step1b$Samples
  SNPFiles <- step1b$SNP
  FullData <- step1b$FullData


  # Step 2 - read in sample  files ------------------------------------------
  lapply(sampleFiles$Path, function(x){
    batch <- sampleFiles$Batch[sampleFiles$Path == x]
    df <- read.table(x, header=T, sep="\t", stringsAsFactors=F) %>%
      dplyr::select(Sample_ID = Sample.ID,
                    CR = Call.Rate,
                    p10_GC = p10.GC,
                    BarcodeID = Array.Info.Sentrix.ID,
                    Position = Array.Info.Sentrix.Position)
    file.name <- paste0("Samples_",batch,".RDS")
    saveRDS(df, file = file.path(home, "Processed Files",file.name))
  })



  # Step 3 - SNP Summary files ----------------------------------------------
  lapply(SNPFiles$Path, function(x){
    batch <- SNPFiles$Batch[SNPFiles$Path == x]
    df <- data.frame(data.table::fread(x)) %>%
      dplyr::select(SNP = Name,
                    Chromosome = Chr,
                    Position,
                    HWE_p = ChiTest100,
                    MAF = Minor.Freq,
                    Alleles = SNP,
                    Strand = Plus.Minus.Strand)

    # I need the genescore off the FullDataTable
    path2 <- dplyr::filter(FullData, Batch == batch)$Path
    df2 <- data.frame(data.table::fread(path2)) %>%
      dplyr::select(SNP = Name,
                    GenTrain = GenTrain.Score)

    # Merge annotation and the GenTrain score
    df3 <- dplyr::left_join(df,df2,"SNP") %>%
      dplyr::select(SNP, Chromosome, Position, HWE_p, MAF, GenTrain,
                    Alleles, Strand)

    file.name <- paste0("SNPTable_",batch,".RDS")
    saveRDS(df3, file = file.path(home, "Processed Files", file.name))
  })



  # Step 4 - full data table ------------------------------------------------

  lapply(FullData$Path, function(x){
    batch <- FullData$Batch[FullData$Path == x]
    df <- data.frame(data.table::fread(x))
    one <- dplyr::select(df,
                         SNP = Name,
                         Chromosome = Chr,
                         Position,
                         GenTrain = GenTrain.Score,
                         Frac_A = Frac.A,
                         Frac_C = Frac.C,
                         Frac_G = Frac.G,
                         Frac_T = Frac.T)
    one$Batch <- batch

    # Now I need to pull out each sample and create a BAF output file

    # 1.  Generate a list of variable names I want to keep based on sample IDs
    samples <- readRDS(file.path(home,"Processed Files",paste0("Samples_",batch,".RDS")))
    ID <- paste0("X",samples$Sample_ID)
    ID <- gsub("-",".",ID)
    variables <- lapply(ID, function(x){
      paste0(x,c(".Score", ".Log.R.Ratio",".B.Allele.Freq"))
    })




    # 2. Create a list object of each of the samples
    samples2 <- lapply(variables, function(x){
      foo <- dplyr::select(df,
                           SNP = Name,
                           Chromosome = Chr,
                           Position,
                           x)
      names(foo) <- c("SNP","Chromosome","Position","GSScore",
                      "LRR","BAF")
      return(foo)
    })
    names(samples2) <- ID

    # 3.  Create an output file for each of these samples
    dir.create(file.path(home,"BAF","Files",batch))


    lapply(names(samples2), function(x){
      samples2[[x]]$ID <- x

      saveRDS(samples2[[x]],
              file = file.path(home, "BAF", "Files",batch,paste0(x,".RDS")))
    })

    # Save the FullData table
    file.name <- paste0("FullDataTable_",batch,".RDS")
    saveRDS(one, file = file.path(home, "Processed Files", file.name))
  })
}



# Program 3 - movePlink() - Organize the PLINK files ------------------------------------
# Process PLINK Files
# Genome Studio outputs all the plink files "PLINK.ped"....
#  I want to rename and reorganize the files to match the batches, then do a quick QC analysis
movePlink <- function(home){

  plink <- file.path(home,"PLINK")


  # Step 1 - each PLINK project has an output file I've named for batch
  #  Let's try to read the batch file and get the directions to the actual data
  files <- list.files(plink,".txt")


  lapply(files, function(x){

    batch <- sub(".txt","",x)

    doc <- read.table(file.path(plink,x), skip = 2, header=F,
                      stringsAsFactors=F, sep="\t")
    path <- doc$V1

    # First remove the PLINK.Bat and PLINK.phenotype file
    file.remove(file.path(path,"PLINK.BAT"))
    file.remove(file.path(path,"PLINK.phenotype"))
    file.remove(file.path(path,"PLINK.script"))

    # List the .ped and .map files in the path and rename tnem
    ped <- list.files(path, ".ped", full.names=F)
    map <- list.files(path, ".map", full.names=F)
    file.rename(from = file.path(path, ped),
                to = file.path(path, paste0(batch,".ped")))
    file.rename(from = file.path(path, map),
                to = file.path(path, paste0(batch,".map")))


    # Copy the renamed plink files into a new location
    dir.create(file.path(plink,batch))
    oldfiles <- list.files(path, full.names=T)
    lapply(oldfiles, function(x) file.copy(from=x, to = file.path(plink, batch), overwrite = T))

    # Remove the old copies
    lapply(oldfiles, file.remove)
    unlink(path, recursive = T, force = T)

    # Copy the original report file to the new location
    file.copy(file.path(plink,paste0(batch,".txt")),
              file.path(plink,batch,paste0(batch,".txt")))
    file.remove(file.path(plink,paste0(batch,".txt")))
  })
}



# Program 4 - PlinkQC() - QC Analysis of PLINK files ----------------------------------
PlinkQC <- function(home){

  plink <- file.path(home,"PLINK")  # location of PLINK directories
  groups <- list.dirs(plink, recursive = F, full.names=F)  # list of projects needing processing

  littleFun <- function(batch){

    # Create and output path for saving PLINK output
    batchdir <- file.path(plink,batch)
    output <- file.path(batchdir,"output")
    dir.create(output)


    # Need a list of files so I can delete everything when I'm done
    keep <- paste0(batch,c(".map",".ped"))


    # Copy files over to output folder
    move <- function(){
      f <- list.files(batchdir)
      f <- f[!f %in% c(keep,"output")]
      lapply(f, function(x){
        file.copy(from = file.path(batchdir,x),
                  to = output,
                  recursive = T)
      })
      lapply(f, function(x) file.remove(file.path(batchdir, x)))
    }


    path2plink <- plinkQC::checkPlink()


    # Convert the ped files to binary
    com <- c("--file", file.path(batchdir,batch),
             "--make-bed", "--out", file.path(batchdir, batch))
    sys::exec_wait(path2plink, com)
    rm(com)


    # Calculate frequencies
    com <- c("--bfile", file.path(batchdir,batch), "--missing", "--out",file.path(batchdir,batch))
    sys::exec_wait(path2plink, com)
    rm(com)

    # Calculate heterozygosity
    com <- c("--bfile", file.path(batchdir,batch), "--het", "--out",file.path(batchdir,batch))
    sys::exec_wait(path2plink, com)
    rm(com)




    # Calculate IBD for all pairs
    com1 <- c("--bfile", file.path(batchdir,batch),
              "--indep-pairwise", 50, 5, .2,
              "--out",
              file.path(batchdir,"LDprune"))
    com2 <- c("--bfile", file.path(batchdir, batch),
              "--extract",
              file.path(batchdir,"LDprune.prune.in"),
              "--genome",
              "--out",
              file.path(batchdir, batch))

    sys::exec_wait(path2plink, com1)
    sys::exec_wait(path2plink, com2)
    rm(com1, com2)


    # Pull together all the Sample QC statistics

    # Call rates
    imiss <- data.table::fread(file.path(batchdir,paste0(batch,".imiss")))
    imiss$Geno_Rate <- 1-imiss$F_MISS
    imiss <- dplyr::select(imiss, IID, Geno_Rate)

    # Heterozygosity
    het <- data.table::fread(file.path(batchdir,paste0(batch,".het")))
    het$HET_Z <- (mean(het$F) - het$F) / sd(het$F) # recodes F to a z score
    het <- dplyr::select(het, IID, HET_Z)

    # IBD pairs
    ibd <- data.table::fread(file.path(batchdir,paste0(batch,".genome")))
    ibd <- dplyr::filter(ibd, PI_HAT > 0.5)
    IID <- unique(ibd$IID1)
    ibd_Pairs <- lapply(IID, function(x){
      df <- dplyr::filter(ibd, IID1 == x)
      HighIBD <- paste(df$IID2, collapse = ", ")
      Pi_Hat <- paste0(df$PI_HAT, collapse = ", ")
      data.frame(IID = x, HighIBD = HighIBD, PI_HAT = Pi_Hat, stringsAsFactors=F)
    }) %>%
      do.call("rbind",.)
    if (is.null(ibd_Pairs)) {
      ibd_Pairs <-
        data.frame(IID = "",
                   HighIBD = "",
                   PI_HAT = "", stringsAsFactors=F)
    }
    rm(ibd, IID)

    samplesQC <- Reduce(function(x,y) dplyr::full_join(x,y,"IID"),
                        list(imiss, het, ibd_Pairs))


    # QC Analysis of the SNPs

    # Allele Frequencies
    com <- c("--bfile", file.path(batchdir,batch), "--freq", "--out", file.path(batchdir,batch))
    sys::exec_wait(path2plink, com)
    rm(com)

    # Hardy-Weinberg
    com <- c("--bfile", file.path(batchdir,batch), "--hardy", "--out", file.path(batchdir,batch))
    sys::exec_wait(path2plink, com)
    rm(com)



    # Process the output files
    maf <- data.table::fread(file.path(batchdir,paste0(batch,".frq"))) %>%
      dplyr::select(SNP, MAF)


    hwe <- data.table::fread(file.path(batchdir,paste0(batch,".hwe")))  %>%
      dplyr::filter(TEST == "ALL(NP)") %>%
      dplyr::select(SNP, HWE_P = P)

    miss <- data.table::fread(file.path(batchdir,paste0(batch,".lmiss"))) %>%
      dplyr::select(SNP, PropMissing = F_MISS)

    SNPQC <- Reduce(function(x,y) dplyr::full_join(x,y,"SNP"),
                    list(maf,hwe,miss))



    # Move everything to the output folder
    move()
    rm(het,hwe,ibd_Pairs,imiss,maf,miss)

    # Flag my problematic samples and SNPs
    samplesQC$Flag <- with(samplesQC, ifelse(
      Geno_Rate < 0.95, 1, ifelse(
        abs(HET_Z) > 3, 2, ifelse(
          !is.na(HighIBD), 3, 4)))) %>%
      factor(levels = 1:4,
             labels = c("Call Rate <0.95",
                        ">3SD Heterozygosity",
                        "High IBD Pair",
                        "Good Sample"))
    samplesQC <- dplyr::filter(samplesQC, IID != "")


    # SNPs
    SNPQC$Flag <- with(SNPQC, ifelse(
      HWE_P < 10e-6, 1, ifelse(
        PropMissing > 0.05, 2, ifelse(
          is.na(MAF), 4, ifelse(
            MAF < 0.01, 3, 5))))) %>%
      factor(levels = 1:5,
             labels = c("HWE < 10e-5",
                        "PropMissing > 0.05",
                        "MAF < 0.01",
                        "Missing MAF",
                        "Good Sample"))
    SNPQC$Flag[is.na(SNPQC$Flag) & is.na(SNPQC$MAF)] <- "Missing MAF"

    SampleSNP_QC <- list(SampleQC=samplesQC, SNPQC=SNPQC)
    save(SampleSNP_QC, file = file.path(plink,paste0(batch,"_QC.rdata")))
  }

  # before I run littlefun() what have i already run?
  finished <- list.files(plink,"QC")
  finished <- sub("_QC","",finished)
  finished <- sub(".rdata","",finished)
  groups <- groups[!groups %in% finished]
  lapply(groups, littleFun)
}



# Program 5 - BAFPlots() - B-allele plots ----------------------------------------------

BAFPlots <- function(home){

  require(ggplot2)

  # Define the file paths
  files <- file.path(home,"BAF","Files")
  plots <- file.path(home,"BAF","Figures")
  processed <- file.path(home,"Processed Files")

  # Figure out the batch organization
  samplesBatch <- data.frame(File = list.files(processed,"Samples", ignore.case=T), stringsAsFactors=F)
  samplesBatch$batch <- sub("Samples_","",samplesBatch$File) %>%
    sub(".RDS","",.)

  # Create the file.paths for each batch
  plotsDir <- list.files(file.path(home,"BAF"), full.names=F, recursive = F) %>%
    sub("_BAF.RDS","",.)
  batch <- samplesBatch$batch[!samplesBatch$batch %in% plotsDir]

  lapply(batch, function(x) suppressWarnings(dir.create(file.path(plots,x))))

  # Get the list of files that need to be processed
  samplesBatch <- dplyr::filter(samplesBatch, batch %in% batch)

  filelist <- split(samplesBatch, samplesBatch$batch)
  filelist <- lapply(filelist, function(x){
    mybatch <- x[["batch"]]
    foo <- readRDS(file.path(processed,x[["File"]]))
    foo$Batch <- mybatch
    foo$Sample_ID <- sub("-",".",foo$Sample_ID)
    return(dplyr::select(foo, Sample_ID, BarcodeID, Batch))
  })
  names(filelist)

  # Subset it only to the files that have not been made
  a <- lapply(list.dirs(plots), function(x){
    b <- data.frame(File = list.files(x), stringsAsFactors=F)
    b$File <- sub(".png","",b$File)
    return(b)
  }) %>% do.call("rbind",.)

  filelist <- lapply(filelist, function(x){
    dplyr::filter(x, !Sample_ID %in% a$File)
  })


  master <- Reduce(function(x,y) rbind(x,y), filelist)

  # Function to make the figures
  makeFigures <- function(sampleID){

    # Get batch
    batch <- dplyr::filter(master, Sample_ID==sampleID)$Batch

    # Load the dataset
    df <- readRDS(file.path(files,batch,paste0("X",sampleID,".RDS")))

    # Flag the sample ID and create a filename
    Sample_ID <- sub("X","",df$ID[1])
    filename <- paste0(Sample_ID,".png")
    title1 <- paste0("B Allele Frequency plots for sample: ",Sample_ID)

    # Drop Chromosome 0 - these are all copy number variants
    df <- dplyr::filter(df, Chromosome != 0)
    df$Chromosome <- factor(df$Chromosome,
                            levels = c(1:22,"X","XY","Y","MT"),
                            labels = c(1:22,"X","XY","Y","MT"))


    # Split up by chromosome so I can recalculate position
    lst <- split(df,df$Chromosome)
    for (i in 2:length(lst)){
      maxBP <- as.numeric(max(lst[[i-1]]$Position))
      lst[[i]]$Position <- as.numeric(lst[[i]]$Position + maxBP)
    }
    df2 <- Reduce(function(x,y) rbind(x,y), lst)


    # Define my messy statistic - only for autosomes
    df2$PropMess <- ifelse(is.na(df2$BAF), NA, ifelse(
      df2$Chromosome %in% c("X","XY","Y","MT"),9, ifelse(
        df2$BAF >0.85 | df2$BAF < .15, 1, ifelse(
          df$BAF > 0.6 & df2$BAF < 0.4, 1, 2))))


    # This is the percentage of SNPs with CNV problems - add it to a subtitle
    CNV <- nrow(df2[df2$PropMess == 2 & !is.na(df2$PropMess),]) / nrow(df2)
    subtitle <- paste0("BAF Mess = ", round(CNV,4))

    # find chromosome position centers for labels
    centers <- sapply(unique(df2$Chromosome), function(x){
      median(df2$Position[df2$Chromosome==x])
    })


    # G1 - create the BAF Plot
    g1 <- ggplot(df2, aes(x = Position, y = BAF)) +
      geom_point(size = 0.01, shape = 16, color = "royalblue") +
      facet_wrap(~ Chromosome, nrow = 1,
                 scales = "free_x",
                 strip.position = "bottom") +
      labs(title = title1, subtitle = subtitle) +
      theme(axis.text.x = element_blank(),
            axis.title = element_blank(),
            plot.title = element_text(size = 10, hjust = 0),
            strip.text = element_text(size = 4),
            panel.spacing = unit(1,"points"))

    # G2 - Creates a histogram to see the distribution of BAF
    g2 <- ggplot(df2, aes(BAF))+ geom_histogram() +
      labs(title = paste0("BAF Histogram")) +
      scale_x_continuous(label = seq(0,1,0.1), breaks = seq(0,1,0.1)) +
      theme(axis.text.x = element_text(size = 10),
            axis.text.y = element_blank(),
            axis.title = element_blank(),
            plot.title = element_text(size = 10, hjust = 0))

    # Arrange them on one page
    g3 <- ggpubr::ggarrange(g1, g2, ncol = 1)

    # Save the figure to the plots folder
    ggplot2::ggsave(file.path(plots,batch,filename), g3, device = "png",
                    height = 6, width = 11, dpi = 250)

    # Return the summary statistic - useful for exclusions
    return(data.frame(Sample_ID = Sample_ID, BafMess = CNV))

  }

  # Run the function for each batch
  BAFStats <- lapply(filelist, function(x){
    foo <- x
    lapply(foo$Sample_ID,function(y){
      makeFigures(y)
    }) %>% do.call("rbind",.)
  })
  # Save summary statistic for each batch
  lapply(names(BAFStats), function(x){
    if (!is.null(BAFStats[[x]])){
    saveRDS(BAFStats[[x]], file=file.path(home,"BAF",paste0(x,"_BAF.RDS")))
    }
      })
}



# Program 6. finalReport() - create summary QC files ----------------------
finalReport <- function(home){

  # Define file paths
  output <- file.path(home,"Summary Reports")
  processed <- file.path(home,"Processed files")
  plink <- file.path(home,"PLINK")
  BAF <- file.path(home,"BAF")


  # Get the batches I need to process
  batches <- list.dirs(plink, recursive=F, full.names=F)

  # Subset out the ones I've already completed
  summaryFiles <- list.files(output) %>%
    sub("_FinalReport.RDS","",.)
  batches <- batches[!batches %in% summaryFiles]

  # Run the reports for all batches
  reports <- lapply(batches, function(x){

    # Filenames I need to load
    samplesdata <- paste0("Samples_",x,".RDS")
    SNPTable <- paste0("SNPTable_",x,".RDS")
    plinkFile <- paste0(x,"_QC.Rdata")
    BAFFile <- paste0(x,"_BAF.RDS")

    # Read in the datasets
    samplesdata <- readRDS(file.path(processed,samplesdata))
    SNPTable <- readRDS(file.path(processed,SNPTable))
    load(file.path(plink,plinkFile)) # called SampleSNP_QC - PLINK Reports
      sampleQC <- SampleSNP_QC$SampleQC
      SNPQC <- SampleSNP_QC$SNPQC
      rm(SampleSNP_QC)
    BAFFile <- readRDS(file.path(BAF,BAFFile))
      BAFFile$Sample_ID <- as.character(BAFFile$Sample_ID)


    # Step 1 - merge the Sample QC files
    sampleQC <- dplyr::full_join(
      dplyr::select(samplesdata, -CR),
      dplyr::select(sampleQC, Sample_ID = IID,
             Geno_Rate,
             Heterozygosity_Z = HET_Z,
             HighIBD, PI_HAT, Flag),
      "Sample_ID") %>%
      dplyr::full_join(.,BAFFile, "Sample_ID")


    # Summary statistics
    averageCR <- mean(sampleQC$Geno_Rate, na.rm=T)
    averagep10GC <- mean(sampleQC$p10_GC, na.rm=T)
    averageBAF <- mean(sampleQC$BafMess, na.rm=T)


    # Exclusion table (add some exclusions to what is already done):
    sampleQC$exclusion <- with(sampleQC, ifelse(
      BafMess > 0.2, 1, ifelse(
        Geno_Rate < 0.95, 2, ifelse(
          abs(Heterozygosity_Z) >= 3, 3, ifelse(
            !is.na(HighIBD), 4, 6)))))
    sampleQC$exclusion[averageCR < 0.95] <- 5
    sampleQC$exclusion <- factor(sampleQC$exclusion,
                                 levels = 1:6,
                                 labels=c("BAF > 0.2", "Call rate <0.95",
                                          ">3SD Heterozygosity",
                                          "High IBD Pair", "Plate fail: Avg CR < 0.95",
                                          "Good sample"))

    # Format an output exclusion table
    exclusion <- table(sampleQC$exclusion, useNA="ifany") %>%
      data.frame(stringsAsFactors=F)
    names(exclusion) <- c("Exclusion", "N")
    exclusion$Exclusion <- as.character(exclusion$Exclusion)
    exclusion <- rbind(c(paste0("Total in batch ",x), sum(exclusion$N)),
                       exclusion)

    sample_exclusion <- data.frame(Group = c("Samples Exclusion",
                                             rep(NA, nrow(exclusion))),
                                   rbind(NA, exclusion),
                                   stringsAsFactors=F)
    rm(exclusion)

    # Step 2 - SNP QC
    SNPQC <- dplyr::full_join(SNPTable,
                       dplyr::select(SNPQC, SNP,
                              PLINK_MAF = MAF,
                              PLINK_HWE = HWE_P,
                              PropMissing,
                              Flag),
                       "SNP")

  # Exclusions
  SNPQC$Exclusion <- with(SNPQC, ifelse(
   GenTrain < averagep10GC, 1, ifelse(
    PLINK_HWE < 1e-04, 2, ifelse(
      PropMissing > 0.05, 3, ifelse(
        MAF < 0.01, 4, 6))))) %>%
  factor(levels = 1:6,
         labels = c("<10th percentile GS", "HWE < 10e-5",
                    "Missing > 0.05", "MAF < 0.01",
                    "Missing metrics","Good sample"))
SNPQC$Exclusion[is.na(SNPQC$Exclusion)] <- "Missing metrics"

    # Format exclusion table
    exclusion <- table(SNPQC$Exclusion, useNA="ifany") %>%
      data.frame(stringsAsFactors=F)
    names(exclusion) <- c("Exclusion", "N")
    exclusion$Exclusion <- as.character(exclusion$Exclusion)
    exclusion <- rbind(c(paste0("Total SNPs in batch ",x), sum(exclusion$N)),
                       exclusion)

    SNP_exclusion <- data.frame(Group = c("SNP Exclusion", rep(NA, nrow(exclusion))),
                                rbind(NA, exclusion),
                                stringsAsFactors=F)

    # Final output report tables
    exclusions <- rbind(sample_exclusion, SNP_exclusion)

    lst <- list(
      SummaryReport = exclusions,
      SampleSummary = sampleQC,
      SNPSummary = SNPQC
    )

    filename <- paste0(x,"_FinalReport.RDS")
    saveRDS(lst, file = file.path(output,filename))
    return(exclusions)  # Just want to return the summary table
  })
  names(reports) <- batches
  return(reports)
}









# Program 7 - QCPipeline() - run the pipeline (wrapper) -------------------
QCPipeline <- function(home){
  # Run the pipeline functions
  makeOrg(home)
  GSProcess(home)
  movePlink(home)
  PlinkQC(home)
  BAFPlots(home)
}

