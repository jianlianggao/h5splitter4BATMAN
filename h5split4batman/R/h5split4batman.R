h5split4batman<-function(BrukerZippedFile, spec_interval)
{
  #written by Dr. Jianliang Gao & Dr. Jie Hao, at Imperial College London, July 2017
  #Splitting NMR big dataset for BATMAN use
  
  warnDef<-options("warn")$warn
  warnRead<-options(warn = -1)
  datapath<-paste(getwd(), "/Bruker", sep="")
  
  unzipfile <- NULL
  ## unzip file
  zipfile <-BrukerZippedFile
  unzip(zipfile, exdir = datapath)
  unzipfile <- substr(zipfile, 1, nchar(zipfile)-4)
  
  
  
  ## read in bruker spectra
  ## find the data files
  ppm <- NULL
  swp <- NULL #for checking the necessary for interpolation. Added by J Gao.
  pfile <-list.files(path = datapath, pattern = "^procs$", all.files = FALSE,full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  rfile <-list.files(path = datapath, pattern = "^1r$", all.files = FALSE,full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  
  ps <- substr(pfile, 1, nchar(pfile)-5)
  rs <- substr(rfile, 1, nchar(rfile)-2)
  
  if (length(setdiff(rs,ps)) == 0 & length(setdiff(ps,rs)) > 0)
  {
    pfile <- paste(rs,"procs",sep = "")
  } 
  
  L<-length(pfile)
  Lr<-length(rfile)
  sa <- NULL
  snam <- NULL
  if (L==0 || Lr==0 || L!=Lr)
  {   
    return (cat("Bruker file does not exist in datapath, or other problems with bruker files...\n"))
  } else {
    #detect the number of rows and cols for matrix to store all spectra
    cols=L+1
    con  <- file(pfile[1], open = "r")
    aLine <- readLines(con, n = -1, warn = FALSE)
    myV <- strsplit(aLine, "=")    
    close(con)
    
    ftsize <- 131072
    for (j in 1:length(myV))
    {
      if (match("##$FTSIZE",myV[[j]][1],nomatch = 0))
      {
        ftsize <- as.numeric(myV[[j]][2])
        break
      }
    }
    
    rows<-ftsize
    
    #pre-allocate memory to sa  ; by J Gao
    sa<-matrix(nrow=rows, ncol=cols)
    
    for (i in 1:L)
    {    
      con  <- file(pfile[i], open = "r")
      aLine <- readLines(con, n = -1, warn = FALSE)
      myV <- strsplit(aLine, "=")    
      close(con)
      
      ftsize <- 70000
      for (j in 1:length(myV))
      {
        if (match("##$OFFSET",myV[[j]][1],nomatch = 0))
        {    offset <- as.numeric(myV[[j]][2])
        }
        if (match("##$SW_p",myV[[j]][1],nomatch = 0))
        {    sw <- as.numeric(myV[[j]][2])
        }
        if (match("##$SF",myV[[j]][1],nomatch = 0))
        {
          sf <- as.numeric(myV[[j]][2])
        }
        if (match("##$SI",myV[[j]][1],nomatch = 0))
        {  
          si <- as.numeric(myV[[j]][2])
        }
        if (match("##$BYTORDP",myV[[j]][1],nomatch = 0))
        {    bytordp <- as.numeric(myV[[j]][2])
        }
        if (match("##$NC_proc",myV[[j]][1],nomatch = 0))
        {
          ncproc <- as.numeric(myV[[j]][2])
        }
        if (match("##$FTSIZE",myV[[j]][1],nomatch = 0))
        {
          ftsize <- as.numeric(myV[[j]][2])
        }
      }
      
      if (bytordp==0){machine_format =  "little"}
      else {machine_format = "big"}
      #read NMR resonance data
      s<-readBin(rfile[i], what="int",n = ftsize, size = 4, signed = T, endian =machine_format)
      s<- ((2^ncproc)* s)
      nspec <- length(s)
      
      tmpppm <- ppm
      tmpswp <- swp  #for checking if necessary to do interpolation. Added by J Gao.
      
      swp <- sw/sf
      dppm <- swp/(nspec-1)
      ppm<-offset
      ppm<-seq(offset,(offset-swp),by=-dppm)
  
      ## interpolation
      if (!is.null(tmpppm))
      {
        #if (length(tmpppm) != length(ppm))
        if ((tmpppm[1] != ppm[1]) || tmpswp != swp || length(tmpppm) != length(ppm))
        {
          sinter <- approx(ppm, s, xout = tmpppm)
          s <- sinter$y
          s[is.na(s)]<-0
          ppm <- tmpppm
          swp <- tmpswp #Added for checking if necessary to do interpolation. by J Gao.
        }
      }
      if (all(is.na(sa[,1]))) {
        sa[,1]<-ppm
      }
      if (i>1) {
        sa[,i]<-s
      }
      
      #sa<- cbind(sa,s)
      ## find corresponding title
      stitle<-paste(substr(rfile[i],1,nchar(rfile[i])-2),"title",sep="")
      if (!file.exists(stitle))
        stitle<-paste(substr(rfile[i],1,nchar(rfile[i])-2),"TITLE",sep="")
      if (file.exists(stitle))
      {
        if (!file.info(stitle)$size == 0)
        {
          con<-file(stitle,open="r")
          ntem <- readLines(con, n = 1, warn = FALSE)
          close(con)
        } else {
          sT <- strsplit(rfile[i], "/")
          sTitle <-sT[[1]]         
          lsT<-length(sTitle)
          if (lsT>4)
            ntem<-paste(sTitle[lsT-4],"_",sTitle[lsT-3],"_",sTitle[lsT-1],sep="")
          else if (lsT>3)
            ntem<-paste(sTitle[lsT-3],"_",sTitle[lsT-1],sep="")
          else if (lsT>=1)
            ntem<-paste(sTitle[lsT-1],sep="")
          else
            ntem<-i
        }
      } else {
        sT <- strsplit(rfile[i], "/")
        sTitle <-sT[[1]]         
        lsT<-length(sTitle)
        if (lsT>4)
          ntem<-paste(sTitle[lsT-4],"_",sTitle[lsT-3],"_",sTitle[lsT-1],sep="")
        else if (lsT>3)
          ntem<-paste(sTitle[lsT-3],"_",sTitle[lsT-1],sep="")
        else if (lsT>=1)
          ntem<-paste(sTitle[lsT-1],sep="")
        else
          ntem<-i
      }
      snam<- cbind(snam, ntem)            
    }
  }
  snam <- cbind("ppm", snam)
  #sa <- cbind(ppm,sa)
  #cat(dim(snam))
  #cat("\n")
  #cat(dim(sa))
  colnames(sa)<- snam
  #warnRead<-options(warn = warnDef)
  #library(h5)
  cat("start saving into .h5 file....\n")
  filename<-"nmrh5.h5"
  nmrfile<-h5file(filename)
  
  nmrfile["nmr/data"]<-sa
  cat("data is saved. \n")
  nmrfile["nmr/meta"]<-snam
  cat("title info is saved. \n")
  h5close(nmrfile)
  cat(".h5 file is saved. \n")
  
  if (length(unzipfile)>0)
  {
    unlink(unzipfile, recursive = TRUE, force = FALSE)
  }
  unlink(datapath, recursive = TRUE, force = FALSE)
  
  
  nmr_fileh5<-h5file(filename)
  total_spec<-dim(nmr_fileh5["nmr/meta"][])[2]-1
  spec_interval<-as.numeric(spec_interval)
  number_files<-floor(total_spec/spec_interval)

  for (ii in seq(1, number_files)){
    split_range<-c(1, seq((ii-1)*spec_interval+2, ii*spec_interval+1))
    nmr_spec<-nmr_fileh5["nmr/data"][,split_range]
    nmr_spec_colname<-nmr_fileh5["nmr/meta"][,split_range]
    colnames(nmr_spec)<-nmr_spec_colname
    new_filenames<-paste("NMRdata",(ii-1)*spec_interval+1,"-",ii*spec_interval, ".txt",sep="")
    write.table(nmr_spec,file=new_filenames,row.names=FALSE,col.names=TRUE,quote=FALSE,sep = "\t")
  }
  h5close(nmr_fileh5)
}
