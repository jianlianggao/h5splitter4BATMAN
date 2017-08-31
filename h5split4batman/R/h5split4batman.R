h5split4batman<-function(filename, spec_interval)
{
  #by Dr. Jianliang Gao at Imperial College London, July 2017
  #for Splitting NMR big dataset in .h5 file format
  require(h5)
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
