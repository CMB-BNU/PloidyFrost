##### function: read single sample coverage files #####
readcov <- function(file.prefix=x)
{
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
  }
  
  if(file.exists(paste0(file.prefix,"_tricov.txt"))){
    trifile <- read.table(paste0(file.prefix,"_tricov.txt"),
                          col.names = c("CovA","CovB","CovC","isStrict","VarType","VarId",
                                        "VarNum","VarDis"))
  }else{
    message(paste0("This file ( ",paste0(file.prefix,"_tricov.txt")," ) does not exists !"))
  }
  
  if(file.exists(paste0(file.prefix,"_tetracov.txt"))){
    tetrafile <- read.table(paste0(file.prefix,"_tetracov.txt"),
                            col.names = c("CovA","CovB","CovC","CovD","isStrict","VarType","VarId",
                                          "VarNum","VarDis"))
  }else{
    message(paste0("This file ( ",paste0(file.prefix,"_tetracov.txt")," ) does not exists !"))
  }
  if(file.exists(paste0(file.prefix,"_pentacov.txt"))){
    pentafile <- read.table(paste0(file.prefix,"_pentacov.txt"),
                            col.names = c("CovA","CovB","CovC","CovD","CovE","isStrict","VarType","VarId",
                                          "VarNum","VarDis"))
  }else{
    message(paste0("This file ( ",paste0(file.prefix,"_pentacov.txt")," ) does not exists !"))
  }  
  
  
  cov <- NULL
  
  if(dim(bifile)[1]!=0)
  {
    cov <- rbind(cov,data.frame(coverage=bifile$CovA+bifile$CovB,VarNum=bifile$VarNum,VarSize=bifile$VarType))
  }
  if(dim(trifile)[1]!=0)
  {
    cov <- rbind(cov,data.frame(coverage=trifile$CovA+trifile$CovB+trifile$CovC,VarNum=trifile$VarNum,VarSize=trifile$VarType))
  }
  if(dim(tetrafile)[1]!=0)
  {
    cov <- rbind(cov,data.frame(coverage=tetrafile$CovA+tetrafile$CovB+tetrafile$CovC+tetrafile$CovD,VarNum=tetrafile$VarNum,VarSize=tetrafile$VarType))
  }
  if(dim(pentafile)[1]!=0)
  {
    cov <- rbind(cov,data.frame(coverage=pentafile$CovA+pentafile$CovB+pentafile$CovC+pentafile$CovD+pentafile$CovE,VarNum=pentafile$VarNum,VarSize=pentafile$VarType))
  }
  
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
    fre_all <- rbind(fre_all,data.frame(fre=c(bifre[1,],bifre[2,]),VarNum=rep(bifile$VarNum,time=2),VarSize=bifile$VarType))
    
  }
  if(dim(trifile)[1]!=0)
  {
    trifre <- apply(trifile,1,FUN=function(x){
      cov.sum <- as.numeric(x[1])+as.numeric(x[2])+as.numeric(x[3])
      return (c(as.numeric(x[1])/cov.sum,as.numeric(x[2])/cov.sum,as.numeric(x[3])/cov.sum))
    })
    fre_all <- rbind(fre_all,data.frame(fre=c(trifre[1,],trifre[2,],trifre[3,]),VarNum=rep(trifile$VarNum,time=3),VarSize=trifile$VarType))
    
  }
  if(dim(tetrafile)[1]!=0)
  {
    tetrafre <- apply(tetrafile,1,FUN=function(x){
      cov.sum <- as.numeric(x[1])+as.numeric(x[2])+as.numeric(x[3])+as.numeric(x[4])
      return (c(as.numeric(x[1])/cov.sum,as.numeric(x[2])/cov.sum,as.numeric(x[3])/cov.sum,as.numeric(x[4])/cov.sum))
    })
    fre_all <- rbind(fre_all,data.frame(fre=c(tetrafre[1,],tetrafre[2,],tetrafre[3,],tetrafre[4,]),VarNum=rep(tetrafile$VarNum,time=4),VarSize=tetrafile$VarType))
    
  }
  if(dim(pentafile)[1]!=0)
  {
    pentafile <- apply(pentafile,1,FUN=function(x){
      cov.sum <- as.numeric(x[1])+as.numeric(x[2])+as.numeric(x[3])+as.numeric(x[4])+as.numeric(x[5])
      return (c(as.numeric(x[1])/cov.sum,as.numeric(x[2])/cov.sum,as.numeric(x[3])/cov.sum,as.numeric(x[4])/cov.sum,as.numeric(x[5])/cov.sum))
    })
    fre_all <- rbind(fre_all,data.frame(fre=c(pentafile[1,],pentafile[2,],pentafile[3,],pentafile[4,],pentafile[5,]),VarNum=rep(pentafile$VarNum,time=5),VarSize=pentafile$VarType))
    
  }
  return (list(cov,fre_all))
}

##### function: read multiple samples coverage files #####

colour.readcov <- function(file.prefix=x)
{
  bifile <- NULL
  trifile <- NULL
  tetrafile <- NULL
  pentafile <- NULL
  if(file.exists(paste0(file.prefix,"_bicov.txt"))){
    bifile <- read.table(paste0(file.prefix,"_bicov.txt"),
                         col.names = c("CovA","CovB","color","isStrict","VarType","VarId",
                                       "VarNum","coe","VarDis"),stringsAsFactors = FALSE,
                         colClasses = c("CovA"="numeric","CovB"="numeric"))
  }else{
    message(paste0("This file ( ",paste0(file.prefix,"_bicov.txt")," ) does not exists !"))
  }
  
  if(file.exists(paste0(file.prefix,"_tricov.txt"))){
    trifile <- read.table(paste0(file.prefix,"_tricov.txt"),
                          col.names = c("CovA","CovB","CovC","color","isStrict","VarType","VarId",
                                        "VarNum","coe","VarDis"))
  }else{
    message(paste0("This file ( ",paste0(file.prefix,"_tricov.txt")," ) does not exists !"))
  }
  
  if(file.exists(paste0(file.prefix,"_tetracov.txt"))){
    tetrafile <- read.table(paste0(file.prefix,"_tetracov.txt"),
                            col.names = c("CovA","CovB","CovC","CovD","color","isStrict","VarType","VarId",
                                          "VarNum","coe","VarDis"))
  }else{
    message(paste0("This file ( ",paste0(file.prefix,"_tetracov.txt")," ) does not exists !"))
  }
  
  if(file.exists(paste0(file.prefix,"_pentacov.txt"))){
    pentafile <- read.table(paste0(file.prefix,"_pentacov.txt"),
                            col.names = c("CovA","CovB","CovC","CovD","CovE","color","isStrict","VarType","VarId",
                                          "VarNum","coe","VarDis"))
  }else{
    message(paste0("This file ( ",paste0(file.prefix,"_pentacov.txt")," ) does not exists !"))
  }
  
  
  
  cov <- NULL
  
  if(dim(bifile)[1]!=0)
  {
    cov <- rbind(cov,data.frame(coverage=bifile$CovA+bifile$CovB,VarNum=bifile$VarNum,coe=bifile$coe,color=bifile$color,VarSize=bifile$VarType))
  }
  if(dim(trifile)[1]!=0)
  {
    cov <- rbind(cov,data.frame(coverage=trifile$CovA+trifile$CovB+trifile$CovC,VarNum=trifile$VarNum,coe=trifile$coe,color=trifile$color,VarSize=trifile$VarType))
  }
  if(dim(tetrafile)[1]!=0)
  {
    cov <- rbind(cov,data.frame(coverage=tetrafile$CovA+tetrafile$CovB+tetrafile$CovC+tetrafile$CovD,VarNum=tetrafile$VarNum,coe=tetrafile$coe,color=tetrafile$color,VarSize=tetrafile$VarType))
  }
  if(dim(pentafile)[1]!=0)
  {
    cov <- rbind(cov,data.frame(coverage=pentafile$CovA+pentafile$CovB+pentafile$CovC+pentafile$CovD+pentafile$covE,VarNum=pentafile$VarNum,coe=pentafile$coe,color=pentafile$color,VarSize=pentafile$VarType))
  }  
  
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
    fre_all <- rbind(fre_all,data.frame(fre=c(bifre[1,],bifre[2,]),VarNum=rep(bifile$VarNum,time=2),coe=rep(bifile$coe,time=2),color=rep(bifile$color,time=2),VarSize=bifile$VarType))
    
  }
  if(dim(trifile)[1]!=0)
  {
    trifre <- apply(trifile,1,FUN=function(x){
      cov.sum <- as.numeric(x[1])+as.numeric(x[2])+as.numeric(x[3])
      return (c(as.numeric(x[1])/cov.sum,as.numeric(x[2])/cov.sum,as.numeric(x[3])/cov.sum))
    })
    fre_all <- rbind(fre_all,data.frame(fre=c(trifre[1,],trifre[2,],trifre[3,]),VarNum=rep(trifile$VarNum,time=3),coe=rep(trifile$coe,time=3),color=rep(trifile$color,time=3),VarSize=trifile$VarType))
    
  }
  if(dim(tetrafile)[1]!=0)
  {
    tetrafre <- apply(tetrafile,1,FUN=function(x){
      cov.sum <- as.numeric(x[1])+as.numeric(x[2])+as.numeric(x[3])+as.numeric(x[4])
      return (c(as.numeric(x[1])/cov.sum,as.numeric(x[2])/cov.sum,as.numeric(x[3])/cov.sum,as.numeric(x[4])/cov.sum))
    })
    fre_all <- rbind(fre_all,data.frame(fre=c(tetrafre[1,],tetrafre[2,],tetrafre[3,],tetrafre[4,]),VarNum=rep(tetrafile$VarNum,time=4),coe=rep(tetrafile$coe,time=4),color=rep(tetrafile$color,time=4),VarSize=tetrafile$VarType))
    
  }
  if(dim(pentafile)[1]!=0)
  {
    pentafre <- apply(pentafile,1,FUN=function(x){
      cov.sum <- as.numeric(x[1])+as.numeric(x[2])+as.numeric(x[3])+as.numeric(x[4])+as.numeric(x[5])
      return (c(as.numeric(x[1])/cov.sum,as.numeric(x[2])/cov.sum,as.numeric(x[3])/cov.sum,as.numeric(x[4])/cov.sum,as.numeric(x[5])/cov.sum))
    })
    fre_all <- rbind(fre_all,data.frame(fre=c(pentafre[1,],pentafre[2,],pentafre[3,],pentafre[4,],pentafre[5,]),VarNum=rep(pentafile$VarNum,time=5),coe=rep(pentafile$coe,time=5),color=rep(pentafile$color,time=5),VarSize=pentafile$VarType))
    
  }
  return (list(cov,fre_all))
  
}


#####Analysis of diploid Cyclocarya paliurus SNJ17#####
SNJ17dir <- "SNJ17/PloidyFrost_output/SNJ17"
snj17 <- readcov(SNJ17dir) 
cov <- 39.3 #monoploid coverage c
p <- 2 #ploidy level p

frequency <- snj17[[2]]
frequency$type <- "all"
frequency_num5 <- frequency[frequency$VarNum<=5&frequency$VarSize<=10,]
frequency_num5$type <- "VarNum<=5&VarSize<=10"
frequency_num1 <- frequency[frequency$VarNum==1&frequency$VarSize<=10,]
frequency_num1$type <- "VarNum=1&VarSize<=10"
frequency <- rbind(frequency,frequency_num5)
frequency <- rbind(frequency,frequency_num1)
frequency$type <- factor(frequency$type,levels=c("all","VarNum<=5&VarSize<=10","VarNum=1&VarSize<=10"))

#unfiltered results
coverage_all <- snj17[[1]]
coverage_all$type <- "all"
#after filtering with number of sites in a superbubble larger than 5 or indel size larger than 10
coverage_num5 <- coverage_all[coverage_all$VarNum<=5&coverage_all$VarSize<=10,]
coverage_num5$type <- "VarNum<=5&VarSize<=10"
#after filtering with number of sites in a superbubble larger than 1 or indel size larger than 10 
coverage_num1 <- coverage_all[coverage_all$VarNum==1&coverage_all$VarSize<=10,]
coverage_num1$type <- "VarNum=1&VarSize<=10"
coverage <- coverage_all
coverage <- rbind(coverage,coverage_num5)
coverage <- rbind(coverage,coverage_num1)
coverage$type <- factor(coverage$type,levels=c("all","VarNum<=5&VarSize<=10","VarNum=1&VarSize<=10"))

site.dt <- NULL
site.dt <- data.frame(num.unfiltered=c(nrow(coverage_all)), #total number of sites
                      num.unfiltered.outrange=c(nrow(coverage_all[coverage_all$coverage<(p-1)*cov|coverage_all$coverage>(p+1)*cov,])),
                      num.varnum5=c(nrow(coverage_num5)), # number of remaining sites (VarNum<=5&VarSize<=10)
                      remain.propotion.varnum5=c(nrow(coverage_num5)/nrow(coverage_all)), # propotion of remaining sites in total number (VarNum<=5&VarSize<=10)
                      num.varnum5.filter.outrange=c(nrow(coverage_all[coverage_all$coverage<(p-1)*cov|coverage_all$coverage>(p+1)*cov,])-
                                                      nrow(coverage_num5[coverage_num5$coverage<(p-1)*cov|coverage_num5$coverage>(p+1)*cov,])),
                      filter.propotion.out_range.varnum5=c(1-nrow(coverage_num5[coverage_num5$coverage<(p-1)*cov|coverage_num5$coverage>(p+1)*cov,])/
                                                             nrow(coverage_all[coverage_all$coverage<(p-1)*cov|coverage_all$coverage>(p+1)*cov,])  ),
                      num.varnum1=c(nrow(coverage_num1)), # number of remaining sites (VarNum=1&VarSize<=10)
                      remain.propotion.varnum1=c(nrow(coverage_num1)/nrow(coverage_all)),# propotion of remaining sites in total number (VarNum=1&VarSize<=10)
                      num.varnum1.filter.outrange=c(nrow(coverage_all[coverage_all$coverage<(p-1)*cov|coverage_all$coverage>(p+1)*cov,])-
                                                      nrow(coverage_num1[coverage_num1$coverage<(p-1)*cov|coverage_num1$coverage>(p+1)*cov,])),
                      
                      filter.propotion.out_range.varnum1=c(1-nrow(coverage_num1[coverage_num1$coverage<(p-1)*cov|coverage_num1$coverage>(p+1)*cov,])/
                                                             nrow(coverage_all[coverage_all$coverage<(p-1)*cov|coverage_all$coverage>(p+1)*cov,])  )) 

rownames(site.dt)[nrow(site.dt)] <- "C.paliurus SNJ17"

v <- c()
if(p>=2){
  for(i in 1:(p-1)){
    v <- c(v,(i/p))
  }
}
#sina plot
p1 <- ggplot(data = frequency,aes(x=type,y=fre,color=type))+
  geom_sina(alpha=0.5,shape='.',maxwidth=0.8)+
  scale_y_continuous(breaks = seq(0.2,0.8,0.2),limits = c(0,1))+
  theme_bw()+scale_x_discrete(guide = guide_axis(position = "top"))+
  geom_hline(mapping = aes(yintercept=y),data = data.frame(y=v),linetype=4)+
  labs(y = 'allele frequency',
       x = '')+theme(plot.title = element_text(size = 25, hjust = 0.5,face = "bold"),
                     plot.subtitle = element_text(size = 35, hjust = -0.1,face = "bold"),
                     axis.title = element_text(size=45),
                     axis.text.x = element_text(size=35),  
                     axis.text.y = element_text(size=35),
                     axis.title.x = element_blank(),
                     panel.grid.major.y = element_blank(),
                     panel.grid.minor.y = element_blank(),
                     panel.grid.major.x = element_blank(),
                     panel.grid.minor.x = element_blank(),
                     strip.text.x = element_text(size = 25),
                     strip.text.y = element_text(size = 25),
                     legend.position='none')

p1 <- ggplot(data = frequency,aes(x=fre,y=..density..,fill=type))+
  geom_density(mapping = aes(x=fre))+
  scale_x_continuous(breaks = seq(0.2,0.8,0.2))+
  theme_bw()+facet_grid(~type)+
  geom_vline(mapping = aes(xintercept=x),data = data.frame(x=v),linetype=4)+
  labs(x = 'allele frequency',
       y = 'density')+theme(plot.title = element_text(size = 30, hjust = 0.5,face = "bold"),
                            plot.subtitle = element_text(size = 40, hjust = -0.1,face = "bold"),
                            axis.title = element_text(size=50),
                            axis.text.x = element_text(size=40),  
                            axis.text.y = element_text(size=40),
                            panel.grid.major.y = element_blank(),
                            panel.grid.minor.y = element_blank(),
                            panel.grid.major.x = element_blank(),
                            panel.grid.minor.x = element_blank(),
                            strip.text.x = element_text(size = 30),
                            strip.text.y = element_text(size = 30),
                            legend.position='none')

p2 <- ggplot(data = coverage,aes(x=coverage,y=..count..,fill=type))+
  geom_density(mapping = aes(x=coverage))+
  scale_x_continuous(limits = c(0,quantile(coverage$coverage,0.99)))+
  geom_vline(mapping = aes(xintercept=x),data = data.frame(x=c(cov*(p-1),cov*(p+1))),linetype=4)+
  theme_bw()+facet_grid(~type)+
  labs(x = 'k-mer coverage',
       y = 'count')+theme(plot.title = element_text(size = 30, hjust = 0.5,face = "bold"),
                          plot.subtitle = element_text(size = 40, hjust = -0.1,face = "bold"),
                          axis.title = element_text(size=50),
                          axis.text.x = element_text(size=40),  
                          axis.text.y = element_text(size=40),
                          panel.grid.major.y = element_blank(),
                          panel.grid.minor.y = element_blank(),
                          panel.grid.major.x = element_blank(),
                          panel.grid.minor.x = element_blank(),
                          strip.text.x = element_text(size = 30),
                          strip.text.y = element_text(size = 30),
                          legend.position='none')


all <- c(1.00005,0.356266,0.300164,0.215398,0.177853,0.150212,0.130034,0.114602,0.102419)
s11n6 <- c(0.990932,0.356502,0.298067,0.21587,0.177891,0.150271,0.130087,0.114653,0.102468)
s11n2 <- c(1.09362,0.362962,0.300708,0.2179,0.17947,0.151635,0.131299,0.115754,0.103485)
filter.method <- rep(c("all","VarNum<=5&VarSize<=10","VarNum=1&VarSize<=10"),each=9)
ploidy <- rep(2:10,times=3)
dt <- data.frame(avg.loglikelihood=c(all,s11n6,s11n2),filter=filter.method,ploidy=ploidy)

p3 <- ggplot(dt,aes(x=ploidy,y=avg.loglikelihood,color=filter))+geom_line()+
  geom_hline(aes(color=color,yintercept=y),data.frame(y=c(all[p-1],s11n6[p-1],s11n2[p-1]),color=c("all","VarNum<=5&VarSize<=10","VarNum=1&VarSize<=10")),linetype=3)+
  theme_classic()+
  theme(plot.title = element_text(size = 30, hjust = 0.5,face = "bold"),
        plot.subtitle = element_text(size = 40, hjust = -0.1,face = "bold"),
        axis.title = element_text(size=50),
        axis.text.x = element_text(size=40),  
        axis.text.y = element_text(size=40),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.text.x = element_text(size = 40),
        strip.text.y = element_text(size = 40),
        legend.position = c(0.72,0.8),
        legend.title =element_text(size=40),
        legend.text =element_text(size=35),
        legend.background= element_rect(fill = "transparent",colour = NA))+
  geom_point(size=2)+labs(y="average\nlog-likelihood",color="")+
  scale_x_continuous(breaks = c(2:10))

#####Analysis of tetraploid Cyclocarya paliurus LSX118#####
LSX118dir <- "LSX118/PloidyFrost_output/LSX118"
lsx118 <- readcov(LSX118dir) 
cov <- 18 #monoploid coverage c
p <- 4 #ploidy level p
v <- c()
if(p>=2){
  for(i in 1:(p-1)){
    v <- c(v,(i/p))
  }
}


frequency <- lsx118[[2]]
frequency$type <- "all"
frequency_num5 <- frequency[frequency$VarNum<=5&frequency$VarSize<=10,]
frequency_num5$type <- "VarNum<=5&VarSize<=10"
frequency_num1 <- frequency[frequency$VarNum==1&frequency$VarSize<=10,]
frequency_num1$type <- "VarNum=1&VarSize<=10"
frequency <- rbind(frequency,frequency_num5)
frequency <- rbind(frequency,frequency_num1)
frequency$type <- factor(frequency$type,levels=c("all","VarNum<=5&VarSize<=10","VarNum=1&VarSize<=10"))


#unfiltered results
coverage_all <- lsx118[[1]]
coverage_all$type <- "all"
#after filtering with number of sites in a superbubble larger than 5 or indel size larger than 10
coverage_num5 <- coverage_all[coverage_all$VarNum<=5&coverage_all$VarSize<=10,]
coverage_num5$type <- "VarNum<=5&VarSize<=10"
#after filtering with number of sites in a superbubble larger than 1 or indel size larger than 10 
coverage_num1 <- coverage_all[coverage_all$VarNum==1&coverage_all$VarSize<=10,]
coverage_num1$type <- "VarNum=1&VarSize<=10"
coverage <- coverage_all
coverage <- rbind(coverage,coverage_num5)
coverage <- rbind(coverage,coverage_num1)
coverage$type <- factor(coverage$type,levels=c("all","VarNum<=5&VarSize<=10","VarNum=1&VarSize<=10"))


site.dt <- rbind(site.dt,
                 data.frame(num.unfiltered=c(nrow(coverage_all)), #total number of sites
                            num.unfiltered.outrange=c(nrow(coverage_all[coverage_all$coverage<(p-1)*cov|coverage_all$coverage>(p+1)*cov,])),
                            num.varnum5=c(nrow(coverage_num5)), # number of remaining sites (VarNum<=5&VarSize<=10)
                            remain.propotion.varnum5=c(nrow(coverage_num5)/nrow(coverage_all)), # propotion of remaining sites in total number (VarNum<=5&VarSize<=10)
                            num.varnum5.filter.outrange=c(nrow(coverage_all[coverage_all$coverage<(p-1)*cov|coverage_all$coverage>(p+1)*cov,])-
                                                            nrow(coverage_num5[coverage_num5$coverage<(p-1)*cov|coverage_num5$coverage>(p+1)*cov,])),
                            filter.propotion.out_range.varnum5=c(1-nrow(coverage_num5[coverage_num5$coverage<(p-1)*cov|coverage_num5$coverage>(p+1)*cov,])/
                                                                   nrow(coverage_all[coverage_all$coverage<(p-1)*cov|coverage_all$coverage>(p+1)*cov,])  ),
                            num.varnum1=c(nrow(coverage_num1)), # number of remaining sites (VarNum=1&VarSize<=10)
                            remain.propotion.varnum1=c(nrow(coverage_num1)/nrow(coverage_all)),# propotion of remaining sites in total number (VarNum=1&VarSize<=10)
                            num.varnum1.filter.outrange=c(nrow(coverage_all[coverage_all$coverage<(p-1)*cov|coverage_all$coverage>(p+1)*cov,])-
                                                            nrow(coverage_num1[coverage_num1$coverage<(p-1)*cov|coverage_num1$coverage>(p+1)*cov,])),
                            filter.propotion.out_range.varnum1=c(1-nrow(coverage_num1[coverage_num1$coverage<(p-1)*cov|coverage_num1$coverage>(p+1)*cov,])/
                                                                   nrow(coverage_all[coverage_all$coverage<(p-1)*cov|coverage_all$coverage>(p+1)*cov,])  )) 
                 
)
rownames(site.dt)[nrow(site.dt)] <- "C.paliurus LSX118"


#sina plot
p1 <- ggplot(data = frequency,aes(x=type,y=fre,color=type))+
  geom_sina(alpha=0.5,shape='.',maxwidth=0.8)+
  scale_y_continuous(breaks = seq(0.2,0.8,0.2),limits = c(0,1))+
  theme_bw()+scale_x_discrete(guide = guide_axis(position = "top"))+
  geom_hline(mapping = aes(yintercept=y),data = data.frame(y=v),linetype=4)+
  labs(y = 'allele frequency',
       x = '')+theme(plot.title = element_text(size = 25, hjust = 0.5,face = "bold"),
                     plot.subtitle = element_text(size = 35, hjust = -0.1,face = "bold"),
                     axis.title = element_text(size=45),
                     axis.text.x = element_text(size=35),  
                     axis.text.y = element_text(size=35),
                     axis.title.x = element_blank(),
                     panel.grid.major.y = element_blank(),
                     panel.grid.minor.y = element_blank(),
                     panel.grid.major.x = element_blank(),
                     panel.grid.minor.x = element_blank(),
                     strip.text.x = element_text(size = 25),
                     strip.text.y = element_text(size = 25),
                     legend.position='none')

p1 <- ggplot(data = frequency,aes(x=fre,y=..density..,fill=type))+
  geom_density(mapping = aes(x=fre))+
  scale_x_continuous(breaks = seq(0.2,0.8,0.2))+
  theme_bw()+facet_grid(~type)+
  geom_vline(mapping = aes(xintercept=x),data = data.frame(x=v),linetype=4)+
  labs(x = 'allele frequency',
       y = 'density')+theme(plot.title = element_text(size = 30, hjust = 0.5,face = "bold"),
                            plot.subtitle = element_text(size = 40, hjust = -0.1,face = "bold"),
                            axis.title = element_text(size=50),
                            axis.text.x = element_text(size=40),  
                            axis.text.y = element_text(size=40),
                            panel.grid.major.y = element_blank(),
                            panel.grid.minor.y = element_blank(),
                            panel.grid.major.x = element_blank(),
                            panel.grid.minor.x = element_blank(),
                            strip.text.x = element_text(size = 30),
                            strip.text.y = element_text(size = 30),
                            legend.position='none')

p2 <- ggplot(data = coverage,aes(x=coverage,y=..count..,fill=type))+
  geom_density(mapping = aes(x=coverage))+
  scale_x_continuous(limits = c(0,quantile(coverage$coverage,0.99)))+
  geom_vline(mapping = aes(xintercept=x),data = data.frame(x=c(cov*(p-1),cov*(p+1))),linetype=4)+
  theme_bw()+facet_grid(~type)+
  labs(x = 'k-mer coverage',
       y = 'count')+theme(plot.title = element_text(size = 30, hjust = 0.5,face = "bold"),
                          plot.subtitle = element_text(size = 40, hjust = -0.1,face = "bold"),
                          axis.title = element_text(size=50),
                          axis.text.x = element_text(size=40),  
                          axis.text.y = element_text(size=40),
                          panel.grid.major.y = element_blank(),
                          panel.grid.minor.y = element_blank(),
                          panel.grid.major.x = element_blank(),
                          panel.grid.minor.x = element_blank(),
                          strip.text.x = element_text(size = 30),
                          strip.text.y = element_text(size = 30),
                          legend.position='none')

all <- c(0.227424,0.293669,0.363846,0.288452,0.283436,0.291667,0.315141,0.240532,0.237677)
s11n6 <- c(0.210516,0.295988,0.365156,0.279061,0.281015,0.315482,0.316616,0.237769,0.234488)
s11n2 <- c(0.219527,0.308364,0.380471,0.284347,0.290132,0.307597,0.330356,0.246293,0.243295)
filter.method <- rep(c("all","VarNum<=5&VarSize<=10","VarNum=1&VarSize<=10"),each=9)
p <- 4
ploidy <- rep(2:10,times=3)
dt <- data.frame(avg.loglikelihood=c(all,s11n6,s11n2),filter=filter.method,ploidy=ploidy)

p3 <- ggplot(dt,aes(x=ploidy,y=avg.loglikelihood,color=filter))+geom_line()+
  geom_hline(aes(color=color,yintercept=y),data.frame(y=c(all[p-1],s11n6[p-1],s11n2[p-1]),color=c("all","VarNum<=5&VarSize<=10","VarNum=1&VarSize<=10")),linetype=3)+
  theme_classic()+
  theme(plot.title = element_text(size = 30, hjust = 0.5,face = "bold"),
        plot.subtitle = element_text(size = 40, hjust = -0.1,face = "bold"),
        axis.title = element_text(size=50),
        axis.text.x = element_text(size=40),  
        axis.text.y = element_text(size=40),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.text.x = element_text(size = 40),
        strip.text.y = element_text(size = 40),
        legend.position = c(0.72,0.8),
        legend.title =element_text(size=40),
        legend.text =element_text(size=35),
        legend.background= element_rect(fill = "transparent",colour = NA))+
  geom_point(size=2)+labs(y="average\nlog-likelihood",color="")+
  scale_x_continuous(breaks = c(2:10))

#####Analysis of octoploid Fragaria ¡Á ananassa#####


Fanadir <- "Fananassa/PloidyFrost_output/Fana"
Fananassa <- readcov(Fanadir) 
cov <- 52 #monoploid coverage c
p <- 8 #ploidy level p
v <- c()
if(p>=2){
  for(i in 1:(p-1)){
    v <- c(v,(i/p))
  }
}
frequency <- Fananassa[[2]]
frequency$type <- "all"
frequency_num5 <- frequency[frequency$VarNum<=5&frequency$VarSize<=10,]
frequency_num5$type <- "VarNum<=5&VarSize<=10"
frequency_num1 <- frequency[frequency$VarNum==1&frequency$VarSize<=10,]
frequency_num1$type <- "VarNum=1&VarSize<=10"
frequency <- rbind(frequency,frequency_num5)
frequency <- rbind(frequency,frequency_num1)

frequency$type <- factor(frequency$type,levels=c("all","VarNum<=5&VarSize<=10","VarNum=1&VarSize<=10"))



#unfiltered results
coverage_all <- Fananassa[[1]]
coverage_all$type <- "all"
#after filtering with number of sites in a superbubble larger than 5 or indel size larger than 10
coverage_num5 <- coverage_all[coverage_all$VarNum<=5&coverage_all$VarSize<=10,]
coverage_num5$type <- "VarNum<=5&VarSize<=10"
#after filtering with number of sites in a superbubble larger than 1 or indel size larger than 10 
coverage_num1 <- coverage_all[coverage_all$VarNum==1&coverage_all$VarSize<=10,]
coverage_num1$type <- "VarNum=1&VarSize<=10"
coverage <- coverage_all
coverage <- rbind(coverage,coverage_num5)
coverage <- rbind(coverage,coverage_num1)
coverage$type <- factor(coverage$type,levels=c("all","VarNum<=5&VarSize<=10","VarNum=1&VarSize<=10"))

site.dt <- rbind(site.dt,
                 data.frame(num.unfiltered=c(nrow(coverage_all)), #total number of sites
                            num.unfiltered.outrange=c(nrow(coverage_all[coverage_all$coverage<(p-1)*cov|coverage_all$coverage>(p+1)*cov,])),
                            num.varnum5=c(nrow(coverage_num5)), # number of remaining sites (VarNum<=5&VarSize<=10)
                            remain.propotion.varnum5=c(nrow(coverage_num5)/nrow(coverage_all)), # propotion of remaining sites in total number (VarNum<=5&VarSize<=10)
                            num.varnum5.filter.outrange=c(nrow(coverage_all[coverage_all$coverage<(p-1)*cov|coverage_all$coverage>(p+1)*cov,])-
                                                            nrow(coverage_num5[coverage_num5$coverage<(p-1)*cov|coverage_num5$coverage>(p+1)*cov,])),
                            filter.propotion.out_range.varnum5=c(1-nrow(coverage_num5[coverage_num5$coverage<(p-1)*cov|coverage_num5$coverage>(p+1)*cov,])/
                                                                   nrow(coverage_all[coverage_all$coverage<(p-1)*cov|coverage_all$coverage>(p+1)*cov,])  ),
                            num.varnum1=c(nrow(coverage_num1)), # number of remaining sites (VarNum=1&VarSize<=10)
                            remain.propotion.varnum1=c(nrow(coverage_num1)/nrow(coverage_all)),# propotion of remaining sites in total number (VarNum=1&VarSize<=10)
                            num.varnum1.filter.outrange=c(nrow(coverage_all[coverage_all$coverage<(p-1)*cov|coverage_all$coverage>(p+1)*cov,])-
                                                            nrow(coverage_num1[coverage_num1$coverage<(p-1)*cov|coverage_num1$coverage>(p+1)*cov,])),
                            filter.propotion.out_range.varnum1=c(1-nrow(coverage_num1[coverage_num1$coverage<(p-1)*cov|coverage_num1$coverage>(p+1)*cov,])/
                                                                   nrow(coverage_all[coverage_all$coverage<(p-1)*cov|coverage_all$coverage>(p+1)*cov,])  )) 
                 
)
rownames(site.dt)[nrow(site.dt)] <- "F.ananassa"



#sina plot
p1 <- ggplot(data = frequency,aes(x=type,y=fre,color=type))+
  geom_sina(alpha=0.5,shape='.',maxwidth=0.8)+
  scale_y_continuous(breaks = seq(0.2,0.8,0.2),limits = c(0,1))+
  theme_bw()+scale_x_discrete(guide = guide_axis(position = "top"))+
  geom_hline(mapping = aes(yintercept=y),data = data.frame(y=v),linetype=4)+
  labs(y = 'allele frequency',
       x = '')+theme(plot.title = element_text(size = 25, hjust = 0.5,face = "bold"),
                     plot.subtitle = element_text(size = 35, hjust = -0.1,face = "bold"),
                     axis.title = element_text(size=45),
                     axis.text.x = element_text(size=35),  
                     axis.text.y = element_text(size=35),
                     axis.title.x = element_blank(),
                     panel.grid.major.y = element_blank(),
                     panel.grid.minor.y = element_blank(),
                     panel.grid.major.x = element_blank(),
                     panel.grid.minor.x = element_blank(),
                     strip.text.x = element_text(size = 25),
                     strip.text.y = element_text(size = 25),
                     legend.position='none')

p1 <- ggplot(data = frequency,aes(x=fre,y=..density..,fill=type))+
  geom_density(mapping = aes(x=fre))+
  scale_x_continuous(breaks = seq(0.2,0.8,0.2))+
  theme_bw()+facet_grid(~type)+
  geom_vline(mapping = aes(xintercept=x),data = data.frame(x=v),linetype=4)+
  labs(x = 'allele frequency',
       y = 'density')+theme(plot.title = element_text(size = 30, hjust = 0.5,face = "bold"),
                            plot.subtitle = element_text(size = 40, hjust = -0.1,face = "bold"),
                            axis.title = element_text(size=50),
                            axis.text.x = element_text(size=30),  
                            axis.text.y = element_text(size=30),
                            panel.grid.major.y = element_blank(),
                            panel.grid.minor.y = element_blank(),
                            panel.grid.major.x = element_blank(),
                            panel.grid.minor.x = element_blank(),
                            strip.text.x = element_text(size = 30),
                            strip.text.y = element_text(size = 30),
                            legend.position='none')

p2 <- ggplot(data = coverage,aes(x=coverage,y=..count..,fill=type))+
  geom_density(mapping = aes(x=coverage))+
  scale_x_continuous(limits = c(0,quantile(coverage$coverage,0.99)))+
  geom_vline(mapping = aes(xintercept=x),data = data.frame(x=c(cov*(p-1),cov*(p+1))),linetype=4)+
  theme_bw()+facet_grid(~type)+
  labs(x = 'k-mer coverage',
       y = 'count')+theme(plot.title = element_text(size = 30, hjust = 0.5,face = "bold"),
                          plot.subtitle = element_text(size = 40, hjust = -0.1,face = "bold"),
                          axis.title = element_text(size=50),
                          axis.text.x = element_text(size=30),  
                          axis.text.y = element_text(size=30),
                          panel.grid.major.y = element_blank(),
                          panel.grid.minor.y = element_blank(),
                          panel.grid.major.x = element_blank(),
                          panel.grid.minor.x = element_blank(),
                          strip.text.x = element_text(size = 30),
                          strip.text.y = element_text(size = 30),
                          legend.position='none')

all <- c(0.16978,0.181245,0.248488,0.218651,0.245583,0.225221,0.274842,0.212538,0.205768)
s11n6 <- c(0.143511,0.165359,0.236592,0.205225,0.22885,0.211037,0.267179,0.197518,0.189292)
s11n2 <- c(0.181215,0.173444,0.272752,0.226054,0.268693,0.235357,0.305577,0.238818,0.229267)
filter.method <- rep(c("all","VarNum<=5&VarSize<=10","VarNum=1&VarSize<=10"),each=9)
p <- 8
ploidy <- rep(2:10,times=3)
dt <- data.frame(avg.loglikelihood=c(all,s11n6,s11n2),filter=filter.method,ploidy=ploidy)

p3 <- ggplot(dt,aes(x=ploidy,y=avg.loglikelihood,color=filter))+geom_line()+
  geom_hline(aes(color=color,yintercept=y),data.frame(y=c(all[p-1],s11n6[p-1],s11n2[p-1]),color=c("all","VarNum<=5&VarSize<=10","VarNum=1&VarSize<=10")),linetype=3)+
  theme_classic()+
  theme(plot.title = element_text(size = 30, hjust = 0.5,face = "bold"),
        plot.subtitle = element_text(size = 40, hjust = -0.1,face = "bold"),
        axis.title = element_text(size=50),
        axis.text.x = element_text(size=40),  
        axis.text.y = element_text(size=40),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.text.x = element_text(size = 40),
        strip.text.y = element_text(size = 40),
        legend.position = c(0.7,0.3),
        legend.title =element_text(size=40),
        legend.text =element_text(size=35),
        legend.background= element_rect(fill = "transparent",colour = NA))+
  geom_point(size=2)+labs(y="average\nlog-likelihood",color="")+
  scale_x_continuous(breaks = c(2:10))




#####Analysis of three diploid Cyclocarya paliurus samples SNJ17,SNJ18,SNJ110#####

snjdir <- "multiple_samples_Cpaliurus/PloidyFrost_output/snj"
SNJ <- colour.readcov(snjdir)
cov1 <- 39.3 #monoploid coverage c SNJ17 
cov2 <- 40 #monoploid coverage c SNJ18 
cov3 <- 45.8 #monoploid coverage c SNJ110 

p <- 2 #ploidy level p


v <- c()
if(p>=2){
  for(i in 1:(p-1)){
    v <- c(v,(i/p))
  }
}

frequency <- SNJ[[2]]
frequency$type <- "all"
frequency_num5 <- frequency[frequency$VarNum<=5&frequency$VarSize<=10,]
frequency_num5$type <- "VarNum<=5&VarSize<=10"
frequency_dayu <- frequency[frequency$coe>=0.25,]
frequency_dayu$type <- "Cramer's V >= 0.25"
frequency_xiaoyu <- frequency[frequency$coe<0.25,]
frequency_xiaoyu$type <- "Cramer's V < 0.25"
frequency <- rbind(frequency,frequency_num5)
frequency <- rbind(frequency,frequency_dayu)
frequency <- rbind(frequency,frequency_xiaoyu)
frequency$color[frequency$color==0] <- "SNJ17"
frequency$color[frequency$color==1] <- "SNJ18"
frequency$color[frequency$color==2] <- "SNJ110"
frequency$color <- factor(frequency$color,levels=name)
frequency$type <- factor(frequency$type,levels=c("all","VarNum<=5&VarSize<=10","Cramer's V >= 0.25","Cramer's V < 0.25"))

#unfiltered results
coverage_all <- SNJ[[1]]
coverage_all$type <- "all"
#after filtering with number of sites in a superbubble larger than 5 or indel size larger than 10
coverage_num5 <- coverage_all[coverage_all$VarNum<=5&coverage_all$VarSize<=10,]
coverage_num5$type <- "VarNum<=5&VarSize<=10"
#after filtering with Cramer's V value greater than 0.25
coverage_coe_larger <- coverage_all[coverage_all$coe>=0.25,]
coverage_coe_larger$type <- "Cramer's V >= 0.25"
coverage_coe_smaller <- coverage_all[coverage_all$coe<0.25,]
coverage_coe_smaller$type <- "Cramer's V < 0.25"

coverage <- coverage_all
coverage <- rbind(coverage,coverage_num5)
coverage <- rbind(coverage,coverage_coe_larger)
coverage <- rbind(coverage,coverage_coe_smaller)
coverage$color[coverage$color==0] <- "SNJ17"
coverage$color[coverage$color==1] <- "SNJ18"
coverage$color[coverage$color==2] <- "SNJ110"
coverage$color <- factor(coverage$color,levels=name)
coverage$type <- factor(coverage$type,levels=c("all","VarNum<=5&VarSize<=10","Cramer's V >= 0.25","Cramer's V < 0.25"))


site.dt.multi.samples <- NULL
coverages <- c(cov1,cov2,cov3)

for (i in 0:2){
  cov <- coverages[i+1]
  site.dt.multi.samples <- 
    rbind(site.dt.multi.samples,
          data.frame(num.unfiltered=c(nrow(coverage_all[coverage_all$color==i,])), #total number of sites
                     num.unfiltered.outrange=c(nrow(coverage_all[coverage_all$color==i&(coverage_all$coverage<(p-1)*cov|coverage_all$coverage>(p+1)*cov),])),
                     num.varnum5=c(nrow(coverage_num5[coverage_num5$color==i,])), # number of remaining sites (VarNum<=5&VarSize<=10)
                     remain.propotion.varnum5=c(nrow(coverage_num5[coverage_num5$color==i,])/nrow(coverage_all[coverage_all$color==i,])), # propotion of remaining sites in total number (VarNum<=5&VarSize<=10)
                     num.varnum5.filter.outrange=c(nrow(coverage_all[coverage_all$color==i&(coverage_all$coverage<(p-1)*cov|coverage_all$coverage>(p+1)*cov),])-
                                                     nrow(coverage_num5[coverage_num5$color==i&(coverage_num5$coverage<(p-1)*cov|coverage_num5$coverage>(p+1)*cov),])),
                     filter.propotion.out_range.varnum5=c(1-nrow(coverage_num5[coverage_num5$color==i&(coverage_num5$coverage<(p-1)*cov|coverage_num5$coverage>(p+1)*cov),])/
                                                            nrow(coverage_all[coverage_all$color==i&(coverage_all$coverage<(p-1)*cov|coverage_all$coverage>(p+1)*cov),])),
                     num.coe.larger=c(nrow(coverage_coe_larger[coverage_coe_larger$color==i,])), # number of remaining sites (Cramer's V>0.25)
                     remain.propotion.coe=c(nrow(coverage_coe_larger[coverage_coe_larger$color==i,])/nrow(coverage_all[coverage_all$color==i,])), # propotion of remaining sites in total number (Cramer's V>0.25)
                     num.coe.filter.outrange=c(nrow(coverage_all[coverage_all$color==i&(coverage_all$coverage<(p-1)*cov|coverage_all$coverage>(p+1)*cov),])-
                                                 nrow(coverage_coe_larger[coverage_coe_larger$color==i&(coverage_coe_larger$coverage<(p-1)*cov|coverage_coe_larger$coverage>(p+1)*cov),])),
                     filter.propotion.out_range.coe=c(1-nrow(coverage_coe_larger[coverage_coe_larger$color==i&(coverage_coe_larger$coverage<(p-1)*cov|coverage_coe_larger$coverage>(p+1)*cov),])/
                                                        nrow(coverage_all[coverage_all$color==i&(coverage_all$coverage<(p-1)*cov|coverage_all$coverage>(p+1)*cov),]))
                     
          )
    )
  
}
rownames(site.dt.multi.samples) <- c("SNJ17","SNJ18","SNJ110")


#sina plot
p1 <- ggplot(data = frequency,aes(x=type,y=fre,color=type))+
  geom_sina(alpha=0.5,shape='.',maxwidth=0.8)+
  scale_y_continuous(breaks = seq(0.2,0.8,0.2),limits = c(0,1))+  
  theme_bw()+scale_x_discrete(guide = guide_axis(position = "top"))+facet_grid(color~.)+
  geom_hline(mapping = aes(yintercept=y),data = data.frame(y=v),linetype=4)+
  labs(y = 'allele frequency',
       x = '')+theme(plot.title = element_text(size = 30, hjust = 0.5,face = "bold"),
                     plot.subtitle = element_text(size = 40, hjust = -0.1,face = "bold"),
                     axis.title = element_text(size=50),
                     axis.text.x = element_text(size=40),  
                     axis.text.y = element_text(size=28),
                     panel.grid.major.y = element_blank(),
                     panel.grid.minor.y = element_blank(),
                     panel.grid.major.x = element_blank(),
                     panel.grid.minor.x = element_blank(),
                     strip.text.x = element_text(size = 40),
                     strip.text.y = element_text(size = 40),
                     legend.position='none')

p1 <- ggplot(data = frequency,aes(x=fre,y=..density..,fill=type))+
  geom_density(mapping = aes(x=fre))+
  scale_x_continuous(breaks = seq(0.2,0.8,0.2))+
  theme_bw()+facet_grid(color~type)+
  geom_vline(mapping = aes(xintercept=x),data = data.frame(x=v),linetype=4)+
  labs(x = 'allele frequency',
       y = 'density')+theme(plot.title = element_text(size = 30, hjust = 0.5,face = "bold"),
                            plot.subtitle = element_text(size = 40, hjust = -0.1,face = "bold"),
                            axis.title = element_text(size=50),
                            axis.text.x = element_text(size=34),  
                            axis.text.y = element_text(size=28),
                            panel.grid.major.y = element_blank(),
                            panel.grid.minor.y = element_blank(),
                            panel.grid.major.x = element_blank(),
                            panel.grid.minor.x = element_blank(),
                            strip.text.x = element_text(size = 28),
                            strip.text.y = element_text(size = 28),
                            legend.position='none')


p2 <- ggplot(data = coverage,aes(x=coverage,y=..count..,fill=type))+
  geom_density(mapping = aes(x=coverage))+
  scale_x_continuous(limits = c(0,quantile(coverage$coverage,0.99)))+
  geom_vline(mapping = aes(xintercept=x),data = data.frame(x=c(cov1*(p-1),cov1*(p+1),cov2*(p-1),cov2*(p+1),cov3*(p-1),cov3*(p+1)),color=rep(name,each=2)),linetype=4)+
  theme_bw()+facet_grid(color~type)+
  labs(x = 'k-mer coverage',
       y = 'count')+theme(plot.title = element_text(size = 30, hjust = 0.5,face = "bold"),
                          plot.subtitle = element_text(size = 40, hjust = -0.1,face = "bold"),
                          axis.title = element_text(size=50),
                          axis.text.x = element_text(size=34),  
                          axis.text.y = element_text(size=28),
                          panel.grid.major.y = element_blank(),
                          panel.grid.minor.y = element_blank(),
                          panel.grid.major.x = element_blank(),
                          panel.grid.minor.x = element_blank(),
                          strip.text.x = element_text(size = 28),
                          strip.text.y = element_text(size = 28),
                          legend.position='none')

sample <- c("SNJ17","SNJ18","SNJ110")
s17 <- c(1.0574,0.360985,0.30208,0.217353,0.179301,0.151465,0.13115,0.115621,0.103365,
         1.04441,0.360758,0.300283,0.217412,0.179129,0.151334,0.131033,0.115514,0.103266,
         1.04407,0.362366,0.301053,0.218638,0.179934,0.152036,0.13166,0.116089,0.103802)

s18 <- c(1.07353,0.361553,0.304204,0.217726,0.179726,0.151819,0.131466,0.115911,0.103635,
         1.06002,0.36129,0.302181,0.217742,0.179501,0.151646,0.131311,0.115769,0.103504,
         1.05896,0.362522,0.302633,0.218633,0.180071,0.152146,0.131758,0.116178,0.103885)

s110 <- c(1.12308,0.361355,0.308793,0.216898,0.179665,0.151725,0.131383,0.115832,0.10356,
         1.10662,0.360675,0.307136,0.216847,0.1794,0.151519,0.131199,0.115665,0.103405,
         1.11387,0.362524,0.308679,0.218257,0.180379,0.152369,0.131958,0.11636,0.104053)

filter.method <- rep(rep(c("all","VarNum<=5&VarSize<=10","Cramer's V >= 0.25"),each=9),times=3)
p <- 2
ploidy <- rep(2:10,times=9)
dt <- data.frame(avg.loglikelihood=c(s17,s18,s110),filter=filter.method,ploidy=ploidy,sample=rep(sample,each=27))
dt$sampe <- factor(sample,levels=name,labels=name)
p3 <- ggplot(dt,aes(x=ploidy,y=avg.loglikelihood,color=filter))+geom_line()+
  geom_hline(aes(yintercept=y,color=filter),data.frame(y=c(s17[p-1],s17[p-1+9],s17[p-1+18],
                                                           s18[p-1],s18[p-1+9],s18[p-1+18],
                                                           s110[p-1],s110[p-1+9],s110[p-1+18]),
                                                       filter=rep(c("all","VarNum<=5&VarSize<=10","Cramer's V >= 0.25"),times=3),
                                                       sample=factor(rep(sample,each=3),levels=name,labels=name)),
             linetype=3)+
  facet_grid(factor(sample,levels=name,labels=name)~.)+
  theme_classic()+theme(plot.title = element_text(size = 30, hjust = 0.5,face = "bold"),
                        plot.subtitle = element_text(size = 40, hjust = -0.1,face = "bold"),
                        axis.title = element_text(size=50),
                        axis.text.x = element_text(size=34),  
                        axis.text.y = element_text(size=28),
                        panel.grid.major.y = element_blank(),
                        panel.grid.minor.y = element_blank(),
                        panel.grid.major.x = element_blank(),
                        panel.grid.minor.x = element_blank(),
                        strip.text.x = element_text(size = 40),
                        strip.text.y = element_text(size = 40),
                        legend.position = c(0.7,0.9),
                        legend.title =element_text(size=40),
                        legend.text =element_text(size=35),
                        legend.background= element_rect(fill = "transparent",colour = NA))+
  geom_point(size=2)+labs(y="average\nlog-likelihood",color="")+
  scale_x_continuous(breaks = c(2:10))


