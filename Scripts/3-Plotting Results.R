#create manhattan plots for analyses
library(qqman)
library(ggplot2)
library(dplyr)

#read in the lmm output 
LMMouput.file<-fread(paste0(SNPset,"_Tables/Assoc_files/",SNPset,"_",datafilenames,".assoc.txt"))



# Prepare data for plotting (section needs to work on sending to appropriate location )
# ------------------------------------------------------------------------------
# add linkage group column (chr)
chrm<-as.numeric(gsub("Ha412HOChr|:.+","",LMMouput.file$rs))
LMMouput.file$chrmNumber<-chrm

#threshold set based on goa threshold using the effective nunmber of tests identified from Temme et al. 2020 
#effective number of tests = 20562
#adding suggestive threshold for top 1% of SNPs 
suggestedcutoff <-as.data.frame(quantile(LMMouput.file$p_wald, as.numeric(as.character(0.01))))[1, 1]

#plot Manhattan


png(paste0(SNPset,"_Plots/",SNPset,"_",datafilenames,".LMM.Manhattan.png"), width=11.7,height=8.3,units="in",res=1200)
manhattan(LMMouput.file, chr = "chrmNumber", bp="ps", snp = "rs", p="p_wald",suggestiveline =-log10(suggestedcutoff) ,genomewideline = -log10(0.05/20562), ylim=c(0,6),col=c("green", "purple"))
dev.off()
#plot QQPlot 
png(paste0(SNPset,"_Plots/",SNPset,"_",datafilenames,".LMM.QQ.png"), width=11.7,height=8.3,units="in",res=1200)
qq(LMMouput.file$p_wald)
dev.off()


#read in BSLMMoutput file 

bslmmOutput.file<-fread(paste0(paste0(SNPset,"_Tables/",paste0(SNPset,"_",datafilenames,".Sparseparameters.txt"))),header=T,sep="\t", data.table=F)



# add linkage group column (chr)
chrm<-as.numeric(gsub("Ha412HOChr|:.+","",bslmmOutput.file$rs))
bslmmOutput.file$chrmNumber<-chrm

# leverage custom manhattan plots for BSLMM outputs


don <- bslmmOutput.file %>% 
  
  # Compute chromosome size
  group_by(chrmNumber) %>% 
  summarise(chr_len=max(ps)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(bslmmOutput.file, ., by=c("chrmNumber"="chrmNumber")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chrmNumber, ps) %>%
  mutate( BPcum=ps+tot)


axisdf<- don %>% group_by(chrmNumber) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

#plot as PNG

png(file=paste0(SNPset,"_Plots/",SNPset,"_",datafilenames,".BSLMM.Manhattan.png"), width=11.7,height=8.3,units="in",res=1200)

ggplot(don, aes(x=BPcum, y=gamma)) +
  
  # Show all points
  #change point size based on effectsize
  geom_point( aes(color=as.factor(chrmNumber)), size=((20*don$eff)+(1/abs(log(don$eff))))) +
  scale_color_manual(values = rep(c("green", "purple"), 22 )) +
  
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chrmNumber, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1) ) +     # remove space between plot area and x axis
  
  
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

dev.off()






#read in BSLMMoutput file for pvalues  

PVALUE.bslmmOutput.file<-params<-fread(paste0(paste0(SNPset,"_",datafilenames,"_Pvalue_Tables/",paste0(SNPset,"_",datafilenames,"_Pvalue.Sparseparameters.txt"))),header=T,sep="\t", data.table=F)



# add linkage group column (chr)
chrm<-as.numeric(gsub("Ha412HOChr|:.+","",PVALUE.bslmmOutput.file$rs))
PVALUE.bslmmOutput.file$chrmNumber<-chrm

# leverage custom manhattan plots for BSLMM outputs


don <- PVALUE.bslmmOutput.file %>% 
  
  # Compute chromosome size
  group_by(chrmNumber) %>% 
  summarise(chr_len=max(ps)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(PVALUE.bslmmOutput.file, ., by=c("chrmNumber"="chrmNumber")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chrmNumber, ps) %>%
  mutate( BPcum=ps+tot)


axisdf<- don %>% group_by(chrmNumber) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

#plot as PNG
png(file=paste0(SNPset,"_",datafilenames,"_Pvalue_Plots/",paste0(SNPset,"_",datafilenames,"_Pvalue.bslmm.manhattan.png")), width=11.7,height=8.3,units="in",res=1200)

ggplot(don, aes(x=BPcum, y=gamma)) +
  
  # Show all points
  #change point size based on effectsize
  geom_point( aes(color=as.factor(chrmNumber)), size=((20*don$eff)+(1/abs(log(don$eff))))) +
  scale_color_manual(values = rep(c("green", "purple"), 22 )) +
  
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chrmNumber, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1) ) +     # remove space between plot area and x axis
  
  
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

dev.off()



#read in BSLMMoutput file for Effectsize   

BETA.bslmmOutput.file<-params<-fread(paste0(paste0(SNPset,"_",datafilenames,"_Beta_Tables/",paste0(SNPset,"_",datafilenames,"_Beta.Sparseparameters.txt"))),header=T,sep="\t", data.table=F)



# add linkage group column (chr)
chrm<-as.numeric(gsub("Ha412HOChr|:.+","",BETA.bslmmOutput.file$rs))
BETA.bslmmOutput.file$chrmNumber<-chrm

# leverage custom manhattan plots for BSLMM outputs


don <- BETA.bslmmOutput.file %>% 
  
  # Compute chromosome size
  group_by(chrmNumber) %>% 
  summarise(chr_len=max(ps)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(BETA.bslmmOutput.file, ., by=c("chrmNumber"="chrmNumber")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chrmNumber, ps) %>%
  mutate( BPcum=ps+tot)


axisdf<- don %>% group_by(chrmNumber) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

#plot as PNG
paste0(paste0(SNPset,"_",datafilenames,"_Beta_Plots/",paste0(SNPset,"_",datafilenames,".BETA.bslmm.manhattan.png")))

png(file=paste0(paste0(SNPset,"_",datafilenames,"_Beta_Plots/",paste0(SNPset,"_",datafilenames,".BETA.bslmm.manhattan.png"))), width=11.7,height=8.3,units="in",res=1000)

ggplot(don, aes(x=BPcum, y=gamma)) +
  
  # Show all points
  #change point size based on effectsize
  geom_point( aes(color=as.factor(chrmNumber)), size=((20*don$eff)+(1/abs(log(don$eff))))) +
  scale_color_manual(values = rep(c("green", "purple"), 22 )) +
  
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chrmNumber, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1) ) +     # remove space between plot area and x axis
  
  
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

dev.off()




#read in BSLMMoutput file for Combinatory Effectsize& P-value SNPs 

Combo.bslmmOutput.file<-params<-fread(paste0(paste0(SNPset,"_",datafilenames,"_Combo_Tables/",paste0(SNPset,"_",datafilenames,"_Combo.Sparseparameters.txt"))),header=T,sep="\t", data.table=F)



# add linkage group column (chr)
chrm<-as.numeric(gsub("Ha412HOChr|:.+","",Combo.bslmmOutput.file$rs))
Combo.bslmmOutput.file$chrmNumber<-chrm

# leverage custom manhattan plots for BSLMM outputs


don <- Combo.bslmmOutput.file %>% 
  
  # Compute chromosome size
  group_by(chrmNumber) %>% 
  summarise(chr_len=max(ps)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(Combo.bslmmOutput.file, ., by=c("chrmNumber"="chrmNumber")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chrmNumber, ps) %>%
  mutate( BPcum=ps+tot)


axisdf<- don %>% group_by(chrmNumber) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

#plot as PNG
paste0(paste0(SNPset,"_",datafilenames,"_Combo_Plots/",paste0(SNPset,"_",datafilenames,".Combo.bslmm.manhattan.png")))

png(file=paste0(paste0(SNPset,"_",datafilenames,"_Combo_Plots/",paste0(SNPset,"_",datafilenames,".Combo.bslmm.manhattan.png"))), width=11.7,height=8.3,units="in",res=1000)

ggplot(don, aes(x=BPcum, y=gamma)) +
  
  # Show all points
  #change point size based on effectsize
  geom_point( aes(color=as.factor(chrmNumber)), size=((20*don$eff)+(1/abs(log(don$eff))))) +
  scale_color_manual(values = rep(c("green", "purple"), 22 )) +
  
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chrmNumber, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1) ) +     # remove space between plot area and x axis
  
  
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

dev.off()





