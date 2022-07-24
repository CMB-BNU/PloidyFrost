if(!suppressMessages(require("optparse"))){
  install.packages("optparse")
  require("optparse")
}
if(!suppressMessages(require("ggplot2"))){
  install.packages("ggplot2")
  require("ggplot2")
}
option_list <- list(
  make_option(c("-f", "--file"), type = "character", default = FALSE,
              action = "store", help = "allelic frequency file"),
  make_option(c("-t", "--title"), type = "character", default = "title",
              action = "store", help = "histogram title"),
  make_option(c("-o", "--outprefix"), type = "character", default = "allele_frequency",
              action = "store", help = "histogram output prefix"),
  make_option(c("-p", "--ploidy"), type = "integer", default =0,
              action = "store", help = "ploidy level")

)
opt = parse_args(OptionParser(option_list = option_list, usage = "Rscript Drawfreq.R"))
if(file.exists(opt$f)){
  data <- read.table(opt$f,col.names = c("ratio"))
}else{
  cat(paste0("This file:",opt$f," is not exists!\n"))
  q()
}

v <- c();
if(opt$p>=2){
  for(i in 1:(opt$p-1)){
    v <- c(v,(i/opt$p));
  }
}

ggplot(data = data,aes(x=ratio))+
  geom_density(mapping = aes(x=ratio,y=..density..),fill="#6EBFEC",color="black")+
  scale_x_continuous(breaks = seq(0,1,0.1))+
  theme_bw()+
  labs(title=opt$t, 
       x = 'frequency',
       y = 'density')+
  theme(plot.title = element_text(family = "Times",size = 32, hjust = 0.5,face = "bold.italic"),
        plot.subtitle = element_text(size = 24, hjust = 0.5),
        axis.text = element_text(family = "Times",size=24),  
        axis.title = element_text(family = "Times",size=28),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),)+
  geom_vline(xintercept = v,linetype =2)
suppressMessages(ggsave(paste0(opt$o,"_allele_frequency.png")))
while (!is.null(dev.list()))  
  dev.off()
