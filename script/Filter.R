if(!suppressMessages(require("optparse"))){
  install.packages("optparse")
  require("optparse")
}
option_list <- list(
  make_option(c("-S", "--simple"), default = FALSE,
              action = "store_true", help = "only simple bubble"),
  make_option(c("-o", "--outprefix"), type = "character", default = "filtered",
              action = "store", help = "output prefix"),
  make_option(c("-i", "--inprefix"), type = "character", default = "input",
              action = "store", help = "input prefix"),
  make_option(c("-l", "--low"), type = "integer", default = 0,
              action = "store", help = "lower coverage cutoff value"),  
  make_option(c("-u", "--up"), type = "integer", default = 10000,
              action = "store", help = "upper coverage cutoff value"),  
  make_option(c("-I", "--indel"), default = FALSE,
              action = "store_true", help = "filter indel"),
  make_option(c("-P", "--snp"), default = FALSE,
              action = "store_true", help = "filter snp"),
  make_option(c("-n", "--num"), type = "integer", default = 10000,
              action = "store", help = "VarNum cutoff value"),
  make_option(c("-d", "--distance"), type = "integer", default = -1,
              action = "store", help = "VarDistance cutoff value"),
  make_option(c("-s", "--size"), type = "integer", default = 10000,
              action = "store", help = "VarSize cutoff value"),
  make_option(c("-q", "--frequency"), type = "double", default = 0.05,
              action = "store", help = "frequency range(frequency,1-frequency)")
)
opt = parse_args(OptionParser(option_list = option_list, usage = " Process the coverage file "))
if(opt$frequency>0.5){
  message("frequency should < 0.5 ")
  q()
}
file.prefix <- opt$inprefix
#####read file#####
bifile <- NULL
trifile <- NULL
tetrafile <- NULL
pentafile <- NULL
if(file.exists(paste0(file.prefix,"_bicov.txt"))){
  bifile <- read.table(paste0(file.prefix,"_bicov.txt"),
                       col.names = c("CovA","CovB","isStrict","VarType","VarId",
                                     "VarNum","VarDis"),stringsAsFactors = FALSE,
                       colClasses = c("CovA"="numeric","CovB"="numeric"))
}else{
  message(paste0("This file ( ",paste0(file.prefix,"_bicov.txt")," ) does not exists !"))
  q()
}

if(file.exists(paste0(file.prefix,"_tricov.txt"))){
  trifile <- read.table(paste0(file.prefix,"_tricov.txt"),
                        col.names = c("CovA","CovB","CovC","isStrict","VarType","VarId",
                                      "VarNum","VarDis"))
}else{
  message(paste0("This file ( ",paste0(file.prefix,"_tricov.txt")," ) does not exists !"))
  q()
}

if(file.exists(paste0(file.prefix,"_tetracov.txt"))){
  tetrafile <- read.table(paste0(file.prefix,"_tetracov.txt"),
                          col.names = c("CovA","CovB","CovC","CovD","isStrict","VarType","VarId",
                                        "VarNum","VarDis"))
}else{
  message(paste0("This file ( ",paste0(file.prefix,"_tetracov.txt")," ) does not exists !"))
  q()
}
if(file.exists(paste0(file.prefix,"_pentacov.txt"))){
  pentafile <- read.table(paste0(file.prefix,"_pentacov.txt"),
                          col.names = c("CovA","CovB","CovC","CovD","CovE","isStrict","VarType","VarId",
                                        "VarNum","VarDis"))
}else{
  message(paste0("This file ( ",paste0(file.prefix,"_pentacov.txt")," ) does not exists !"))
  q()
}

#####process file#####
out_bifile <- NULL
out_trifile <- NULL
out_tetrafile <- NULL
out_pentafile <- NULL

if(opt$simple){
  bifile <- bifile[bifile$isStrict==1,]
  trifile <- trifile[trifile$isStrict==1,]
  tetrafile <- tetrafile[tetrafile$isStrict==1,]
  pentafile <- pentafile[pentafile$isStrict==1,]
}
if(opt$indel){
  bifile <- bifile[bifile$VarType==0,]
  trifile <- trifile[trifile$VarType==0,]
  tetrafile <- tetrafile[tetrafile$VarType==0,]  
  pentafile <- pentafile[pentafile$VarType==0,]
  
}
if(opt$snp){
  bifile <- bifile[bifile$VarType>0,]
  trifile <- trifile[trifile$VarType>0,]
  tetrafile <- tetrafile[tetrafile$VarType>0,]  
  pentafile <- pentafile[pentafile$VarType>0,]
  
}
bifile <- bifile[bifile$CovA>opt$low&bifile$CovB>opt$low&bifile$CovA<opt$up&bifile$CovB<opt$up&
                   bifile$VarNum<opt$num&bifile$VarDis>opt$distance&bifile$VarType<opt$size,]
trifile <- trifile[trifile$CovA>opt$low&trifile$CovB>opt$low&trifile$CovC>opt$low&trifile$CovA<opt$up&trifile$CovB<opt$up&trifile$CovC<opt$up&
                     trifile$VarNum<opt$num&trifile$VarDis>opt$distance&trifile$VarType<opt$size,]
tetrafile <- tetrafile[tetrafile$CovA>opt$low&tetrafile$CovB>opt$low&tetrafile$CovC>opt$low&tetrafile$CovD>opt$low&tetrafile$CovA<opt$up&tetrafile$CovB<opt$up&tetrafile$CovC<opt$up&tetrafile$CovD<opt$up&
                         (tetrafile$CovA+tetrafile$CovB+tetrafile$CovC+tetrafile$CovD)<opt$up&
                         tetrafile$VarNum<opt$num&tetrafile$VarDis>opt$distance&tetrafile$VarType<opt$size,]
pentafile <- pentafile[pentafile$CovA>opt$low&pentafile$CovB>opt$low&pentafile$CovC>opt$low&
                         pentafile$CovD>opt$low&pentafile$CovE>opt$low&pentafile$CovA<opt$up&
                         pentafile$CovB<opt$up&pentafile$CovC<opt$up&pentafile$CovD<opt$up&pentafile$CovE<opt$up&
                         (pentafile$CovA+pentafile$CovB+pentafile$CovC+pentafile$CovD)<opt$up&
                         pentafile$VarNum<opt$num&pentafile$VarDis>opt$distance&pentafile$VarType<opt$size,]
write.table(bifile,file=paste0(opt$outprefix,"_bicov.txt"),col.names = FALSE,sep = "\t",quote = FALSE,row.names = FALSE)
write.table(trifile,file=paste0(opt$outprefix,"_tricov.txt"),col.names = FALSE,sep = "\t",quote = FALSE,row.names = FALSE)
write.table(tetrafile,file=paste0(opt$outprefix,"_tetracov.txt"),col.names = FALSE,sep = "\t",quote = FALSE,row.names = FALSE)
write.table(pentafile,file=paste0(opt$outprefix,"_pentacov.txt"),col.names = FALSE,sep = "\t",quote = FALSE,row.names = FALSE)

fre_all <- NULL
bifre <- c()
trifre <- c()
tetrafre <- c()
pentafre <- c()
if(dim(bifile)[1]!=0)
{
  bifre <- apply(bifile,1,FUN=function(x){
    cov.sum <- as.numeric(x[1])+as.numeric(x[2])
    return (c(as.numeric(x[1])/cov.sum,as.numeric(x[2])/cov.sum))
  })
  fre_all <- c(fre_all,c(bifre[1,],bifre[2,]))
}
if(dim(trifile)[1]!=0)
{
  trifre <- apply(trifile,1,FUN=function(x){
    cov.sum <- as.numeric(x[1])+as.numeric(x[2])+as.numeric(x[3])
    return (c(as.numeric(x[1])/cov.sum,as.numeric(x[2])/cov.sum,as.numeric(x[3])/cov.sum))
  })
  fre_all <- c(fre_all,c(trifre[1,],trifre[2,],trifre[3,]))
  
}
if(dim(tetrafile)[1]!=0)
{
  tetrafre <- apply(tetrafile,1,FUN=function(x){
    cov.sum <- as.numeric(x[1])+as.numeric(x[2])+as.numeric(x[3])+as.numeric(x[4])
    return (c(as.numeric(x[1])/cov.sum,as.numeric(x[2])/cov.sum,as.numeric(x[3])/cov.sum,as.numeric(x[4])/cov.sum))
  })
  fre_all <- c(fre_all,c(tetrafre[1,],tetrafre[2,],tetrafre[3,],tetrafre[4,]))
  
}
if(dim(pentafile)[1]!=0)
{
  pentafre <- apply(pentafile,1,FUN=function(x){
    cov.sum <- as.numeric(x[1])+as.numeric(x[2])+as.numeric(x[3])+as.numeric(x[4])+as.numeric(x[5])
    return (c(as.numeric(x[1])/cov.sum,as.numeric(x[2])/cov.sum,as.numeric(x[3])/cov.sum,as.numeric(x[4])/cov.sum,as.numeric(x[5])/cov.sum))
  })
  fre_all <- c(fre_all,c(pentafre[1,],pentafre[2,],pentafre[3,],pentafre[4,],pentafre[5,]))
  
}
write.table(round(fre_all[fre_all>opt$frequency&fre_all<(1-opt$frequency)],7),file=paste0(opt$outprefix,"_allele_frequency.txt"),col.names = FALSE,row.names = FALSE)
