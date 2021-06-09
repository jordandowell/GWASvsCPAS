##### Sunflower SAM population GWAS using GEMMA -
#####   a reimagined pipeline borowing some ideas from the version in Masalia et al 2018 (PLOS).
##### Genome: HA412HO
##### SNPs: XRQv1 SNP set built by Greg Owens (UBC) but reordered to HA412HO
##### Algorithm: GEMMA
library(data.table)

#### read in preferences

  SNPset<-NEWSNPS
  pheno.name<-datafiles



## Read in traits and environments to run
#traits<- as.character(read.table("traits_to_run.txt")[,1])

#envs<-as.character(read.table("environments_to_run.txt")[,1])


#create directories as necessary
pheno.data<-fread(paste("data/Phenotypedata/",pheno.name,sep=""))

dir.create(paste0(SNPset,"_Plots/"))
dir.create(paste0(SNPset,"_Tables/"))
dir.create((paste0(SNPset,"_Tables/Assoc_files/")))




 #if the phenotype does not exist go to the next one

trait.data<-pheno.data



#rewrite the fam file to include the new phenotypic information 
fam.file<-fread(paste0(datafilepath ,SNPset,".fam"))
fam.file$V6<-NULL
fam.file <- merge(fam.file,trait.data,by.x="V1",by.y="Line",all.x=T)
write.table(file=paste0(datafilepath ,SNPset,".fam"),fam.file,col.names=F, row.names=F, quote =F)



#using GEMMA to run a bayesian sparse linear mixed model on the full SNP set 
system(paste("Software/gemma -bslmm 1 -bfile ",paste0(datafilepath,SNPset)," -w 250000"," -rpace 100 -wpace 1000"," -k ",paste0(datafilepath,SNPset,".cXX.txt")," -outdir", paste0(SNPset,"_Tables/"), "-o " ,paste0(SNPset,"_",datafilenames)))



# Examine Hyperparameters for bayesian sparse linear mixed model of KEGG SNPS
# Load hyperparameter file for p-value
# ==============================================================================
hyp.params<-read.table(paste0(paste0(SNPset,"_Tables/",paste0(SNPset,"_",datafilenames,".hyp.txt"))),header = T)

# ==============================================================================
# Get mean, median, and 95% ETPI of hyperparameters
# ==============================================================================
# pve -> proportion of phenotypic variance explained by the genotypes
pve<-c("PVE", mean(hyp.params$pve),quantile(hyp.params$pve, probs=c(0.5,0.025,0.975)))

# pge -> proportion of genetic variance explained by major effect loci
pge<-c("PGE",mean(hyp.params$pge),quantile(hyp.params$pge, probs=c(0.5,0.025,0.975)))

# pi -> proportion of variants with non-zero effects
pi<-c("pi",mean(hyp.params$pi),quantile(hyp.params$pi, probs=c(0.5,0.025,0.975)))

# n.gamma -> number of variants with major effect
n.gamma<-c("n.gamma",mean(hyp.params$n_gamma),quantile(hyp.params$n_gamma, probs=c(0.5,0.025,0.975)))
# rho -> approximation to proportion of genetic variance explained by variants with major effects,
#rho close to zero trait is highly polygenic
rho<-c("rho",mean(hyp.params$rho),quantile(hyp.params$rho, probs=c(0.5,0.025,0.975)))


# ==============================================================================

# get table of hyperparameters
# ==============================================================================
hyp.params.table<-as.data.frame(rbind(pve,pge,pi,n.gamma,rho),row.names=F)
colnames(hyp.params.table)<-c("hyperparam", "mean","median","2.5%", "97.5%")
# show table
View(hyp.params.table)
write.table(hyp.params.table,paste0(paste0(SNPset,"_Tables/",paste0(SNPset,"_",datafilenames,".Hyperparametertable.txt"))))
paste0(paste0(SNPset,"_",datafilenames,"_Pvalue_Tables/",paste0(SNPset,"_",datafilenames,"_Pvalue.Hyperparametertable.txt")))
#
# plot traces and distributions of hyperparameters
# ==============================================================================
pdf(file=paste0(paste0(SNPset,"_Tables/",paste0(SNPset,"_",datafilenames,".Hyperparameters.pdf"))), width=8.3,height=11.7)
layout(matrix(c(1,1,2,3,4,4,5,6), 4, 2, byrow = TRUE))

# PVE
# ------------------------------------------------------------------------------
plot(hyp.params$pve, type="l", ylab="PVE", main="PVE - trace")
hist(hyp.params$pve, main="PVE - posterior distribution", xlab="PVE")
plot(density(hyp.params$pve), main="PVE - posterior distribution", xlab="PVE")
# ------------------------------------------------------------------------------

# PGE
# ------------------------------------------------------------------------------
plot(hyp.params$pge, type="l", ylab="PGE", main="PGE - trace")
hist(hyp.params$pge, main="PGE - posterior distribution", xlab="PGE")
plot(density(hyp.params$pge), main="PGE - posterior distribution", xlab="PGE")
# ------------------------------------------------------------------------------

# pi 
# ------------------------------------------------------------------------------
plot(hyp.params$pi, type="l", ylab="pi", main="pi")
hist(hyp.params$pi, main="pi", xlab="pi")
plot(density(hyp.params$pi), main="pi", xlab="pi")
# ------------------------------------------------------------------------------

# N_gamma #number of large effect snps
# ------------------------------------------------------------------------------
plot(hyp.params$n_gamma, type="l", ylab="n_gamma", main="n_gamma - trace")
hist(hyp.params$n_gamma, main="n_gamma - posterior distribution", xlab="n_gamma")
plot(density(hyp.params$n_gamma), main="n_gamma - posterior distribution", xlab="n_gamma")
# ------------------------------------------------------------------------------

# rho
# ------------------------------------------------------------------------------
plot(hyp.params$rho, type="l", ylab="rho", main="rho - trace")
hist(hyp.params$rho, main="rho - posterior distribution", xlab="rho")
plot(density(hyp.params$rho), main="rho - posterior distribution", xlab="rho")
# ------------------------------------------------------------------------------
#h
# ------------------------------------------------------------------------------
plot(hyp.params$h, type="l", ylab="h", main="h - trace")
hist(hyp.params$h, main="h - posterior distribution", xlab="h")
plot(density(hyp.params$h), main="h - posterior distribution", xlab="h")





dev.off()
# ==============================================================================


# Load parameters output
# ==============================================================================
params<-fread(paste0(paste0(SNPset,"_Tables/",paste0(SNPset,"_",datafilenames,".param.txt"))),header=T,sep="\t", data.table=F)
# ==============================================================================
# Get variants with sparse effect size on phenotypes 
# ==============================================================================
# add sparse effect size (= beta * gamma) to data frame
params["eff"]<-abs(params$beta*params$gamma)

# get variants with effect size > 0
params.effects<-params[params$eff>0,]

# show number of variants with measurable effect
nrow(params.effects)

# sort by descending effect size
params.effects.sort<-params.effects[order(-params.effects$eff),]

# show top 10 variants with highest effect
head(params.effects.sort, 10) 

#write file with sparse effect calculated 

write.table(params,file=paste0(SNPset,"_Tables/",SNPset,"_",datafilenames,".Sparseparameters.txt"), quote=F, row.names=F, sep="\t")
# variants with the highest sparse effects
# ------------------------------------------------------------------------------
# top 1% variants (above 99% quantile)
top1<-params.effects.sort[params.effects.sort$eff>quantile(params.effects.sort$eff,0.99),]
# top 0.1% variants (above 99.9% quantile)
top01<-params.effects.sort[params.effects.sort$eff>quantile(params.effects.sort$eff,0.999),]
# top 0.01% variants (above 99.99% quantile)
top001<-params.effects.sort[params.effects.sort$eff>quantile(params.effects.sort$eff,0.9999),]
# ------------------------------------------------------------------------------

# write tables
write.table(top1, file=paste0(paste0(SNPset,"_Tables/",paste0(SNPset,"_",datafilenames,".top1eff.dsv"))), quote=F, row.names=F, sep="\t")
write.table(top01, file=paste0(paste0(SNPset,"_Tables/",paste0(SNPset,"_",datafilenames,".top0.1eff.dsv"))), quote=F, row.names=F, sep="\t")
write.table(top001, file=paste0(paste0(SNPset,"_Tables/",paste0(SNPset,"_",datafilenames,".top0.01eff.dsv"))), quote=F, row.names=F, sep="\t")


# ==============================================================================
# Get variants with high Posterior Inclusion Probability (PIP) == gamma
# ==============================================================================
# PIP is the frequency a variant is estimated to have a sparse effect in the MCMC

# sort variants by descending PIP
params.pipsort<-params[order(-params$gamma),]

# Show top 10 variants with highest PIP
head(params.pipsort,10)

# sets of variants above a certain threshold
# variants with effect in 1% MCMC samples or more
pip01<-params.pipsort[params.pipsort$gamma>=0.01,]
# variants with effect in 10% MCMC samples or more
pip10<-params.pipsort[params.pipsort$gamma>=0.10,]
# variants with effect in 25% MCMC samples or more
pip25<-params.pipsort[params.pipsort$gamma>=0.25,]
# variants with effect in 50% MCMC samples or more
pip50<-params.pipsort[params.pipsort$gamma>=0.50,]

# write tables
write.table(pip01, file=paste0(paste0(SNPset,"_Tables/",paste0(SNPset,"_",datafilenames,".pip01.dsv"))), quote=F, row.names=F, sep="\t")
write.table(pip10, file=paste0(paste0(SNPset,"_Tables/",paste0(SNPset,"_",datafilenames,".pip10.dsv"))), quote=F, row.names=F, sep="\t")
write.table(pip25, file=paste0(paste0(SNPset,"_Tables/",paste0(SNPset,"_",datafilenames,".pip25.dsv"))), quote=F, row.names=F, sep="\t")
write.table(pip50, file=paste0(paste0(SNPset,"_Tables/",paste0(SNPset,"_",datafilenames,".pip50.dsv"))), quote=F, row.names=F, sep="\t")
# ------------------------------------------------------------------------------
# ==============================================================================





#using GEMMA to run a  linear mixed model 
system(paste("Software/gemma -bfile ",paste0(datafilepath,SNPset)," -k ",paste0(datafilepath,SNPset,".cXX.txt"), "-c ",paste0(datafilepath,SNPset,".PCA_EV"), "-lmm 1 -outdir ",paste0(SNPset,"_Tables/Assoc_files/"), "-o ",paste0(SNPset,"_",datafilenames)))





#read in outputs to make 3 datasets top 1% of SNPs based on p-value ,effect size, and a combination of both newly formed datasets


#read in files tped, map, & lmm output
#we are making an assumption that the SNPlist and SNPvariants are in the same order
SNPMap<-read.table(paste0(datafilepath,SNPset,".map"))
SNPTped<- read.table(paste0(datafilepath,SNPset,".tped"))



#read in the lmm output 
LMMouput.file<-fread(paste0(SNPset,"_Tables/Assoc_files/",SNPset,"_",datafilenames,".assoc.txt"))



#get dimensions

#get 1% number of SNPS
OnePercent<- round(dim(LMMouput.file)[1]*0.01,1)

#sort data frame by pvalue 
PvalueLMM<- LMMouput.file[order(LMMouput.file$p_wald),]


#subset top 1% store just names 

PvalueMarkernames<-PvalueLMM$rs[1:OnePercent]



#subset map tped file based on stored names rename and write to folder 

Pvalue_SNPMap<- SNPMap[SNPMap$V2%in%PvalueMarkernames,]
Pvalue_SNPTped<- SNPTped[SNPTped$V2%in%PvalueMarkernames,]

#write top pvalue SNPs to a suitible folder 
write.table(Pvalue_SNPMap, paste0("data/PvalueSNPS","/",SNPset,"_",datafilenames,"_Pvalue.map"), col.names = F,row.names = F,quote = F, sep = "\t")
write.table(Pvalue_SNPTped, paste0("data/PvalueSNPS","/",SNPset,"_",datafilenames,"_Pvalue.tped"), col.names = F,row.names = F,quote = F)



#sort data frame by beta(effectsize) 

BetaLMM<- LMMouput.file[order(abs(LMMouput.file$beta)),]
#subset top 1% store just names 
BetaMarkernamers<-BetaLMM$rs[1:OnePercent]

#subset tped file based on stored names; rename and write to folder 

Beta_SNPMap<- SNPMap[SNPMap$V2%in%BetaMarkernamers,]
Beta_SNPTped<- SNPTped[SNPTped$V2%in%BetaMarkernamers,]

#write top effectsize SNPs to a suitable folder 
write.table(Beta_SNPMap, paste0("data/EffectsizeSNPS","/",SNPset,"_",datafilenames,"_Beta.map"), col.names = F,row.names = F,quote = F, sep = "\t")
write.table(Beta_SNPTped, paste0("data/EffectsizeSNPS","/",SNPset,"_",datafilenames,"_Beta.tped"), col.names = F,row.names = F,quote = F)



#subset marker lists the SNPs are in the top 1% of betas or pvalues

pval_Beta_SNPMap<- SNPMap[SNPMap$V2%in%BetaMarkernamers|SNPMap$V2%in%PvalueMarkernames,]


pval_Beta_SNPTped<- SNPTped[SNPTped$V2%in%BetaMarkernamers|SNPTped$V2%in%PvalueMarkernames,]


#write top pvalue SNPs to a suitible folder 
write.table(pval_Beta_SNPMap, paste0("data/ComboSNPS","/",SNPset,"_",datafilenames,"_PvalueBeta.map"), col.names = F,row.names = F,quote = F, sep = "\t")
write.table(pval_Beta_SNPTped, paste0("data/ComboSNPS","/",SNPset,"_",datafilenames,"_PvalueBeta.tped"), col.names = F,row.names = F,quote = F)

###add a tfam file to each folder make sure there names are equal afterwards



#use plink to make bed files for each folder 
Pvalue_filePath<-"data/PvalueSNPS/"
PvalueSNPSET<-paste0(SNPset,"_",datafilenames,"_Pvalue")


system(paste("Software/plink --tfile",paste0(Pvalue_filePath, PvalueSNPSET),"--covar",paste0(datafilepath, NEWSNPS,"_COVARIATES.txt") , "--make-bed", "--allow-extra-chr", "--out ",paste0(Pvalue_filePath, PvalueSNPSET)))

#use plink to make bed files for each folder 
Beta_filePath<-"data/EffectsizeSNPS/"
BetaSNPSET<-paste0(SNPset,"_",datafilenames,"_Beta")


system(paste("Software/plink --tfile",paste0(Beta_filePath, BetaSNPSET),"--covar",paste0(datafilepath, NEWSNPS,"_COVARIATES.txt") , "--make-bed", "--allow-extra-chr", "--out ",paste0(Beta_filePath, BetaSNPSET)))

#use plink to make bed files for each folder 
Combo_filePath<-"data/ComboSNPS/"
ComboSNPSET<-paste0(SNPset,"_",datafilenames,"_PvalueBeta")


system(paste("Software/plink --tfile",paste0(Combo_filePath, ComboSNPSET),"--covar",paste0(datafilepath, NEWSNPS,"_COVARIATES.txt") , "--make-bed", "--allow-extra-chr", "--out ",paste0(Combo_filePath, ComboSNPSET)))




#data is now ready to be used for bslmm
#create a directory to store bslmm outputs

dir.create(paste0(SNPset,"_",datafilenames,"_Pvalue_Plots/"))
dir.create(paste0(SNPset,"_",datafilenames,"_Pvalue_Tables/"))

#USE GEMMA to run a bayesian sparse linear mixed model 
system(paste("Software/gemma -bslmm 1 -bfile ",paste0(Pvalue_filePath, PvalueSNPSET)," -w 250000"," -rpace 100 -wpace 1000"," -k ",paste0(datafilepath,SNPset,".cXX.txt")," -outdir", paste0(SNPset,"_",datafilenames,"_Pvalue_Tables/"), "-o " ,paste0(SNPset,"_",datafilenames,"_Pvalue")))




# Load hyperparameter file for p-value
# ==============================================================================
hyp.params<-read.table(paste0(paste0(SNPset,"_",datafilenames,"_Pvalue_Tables/",paste0(SNPset,"_",datafilenames,"_Pvalue.hyp.txt"))),header = T)

# ==============================================================================
# Get mean, median, and 95% ETPI of hyperparameters
# ==============================================================================
# pve -> proportion of phenotypic variance explained by the genotypes
pve<-c("PVE", mean(hyp.params$pve),quantile(hyp.params$pve, probs=c(0.5,0.025,0.975)))

# pge -> proportion of genetic variance explained by major effect loci
pge<-c("PGE",mean(hyp.params$pge),quantile(hyp.params$pge, probs=c(0.5,0.025,0.975)))

# pi -> proportion of variants with non-zero effects
pi<-c("pi",mean(hyp.params$pi),quantile(hyp.params$pi, probs=c(0.5,0.025,0.975)))

# n.gamma -> number of variants with major effect
n.gamma<-c("n.gamma",mean(hyp.params$n_gamma),quantile(hyp.params$n_gamma, probs=c(0.5,0.025,0.975)))
# rho -> approximation to proportion of genetic variance explained by variants with major effects,
#rho close to zero trait is highly polygenic
rho<-c("rho",mean(hyp.params$rho),quantile(hyp.params$rho, probs=c(0.5,0.025,0.975)))


# ==============================================================================

# get table of hyperparameters
# ==============================================================================
hyp.params.table<-as.data.frame(rbind(pve,pge,pi,n.gamma,rho),row.names=F)
colnames(hyp.params.table)<-c("hyperparam", "mean","median","2.5%", "97.5%")
# show table
View(hyp.params.table)
write.table(hyp.params.table,paste0(paste0(SNPset,"_",datafilenames,"_Pvalue_Tables/",paste0(SNPset,"_",datafilenames,"_Pvalue.Hyperparametertable.txt"))))
paste0(paste0(SNPset,"_",datafilenames,"_Pvalue_Tables/",paste0(SNPset,"_",datafilenames,"_Pvalue.Hyperparametertable.txt")))
#
# plot traces and distributions of hyperparameters
# ==============================================================================
pdf(file=paste0(paste0(SNPset,"_",datafilenames,"_Pvalue_Tables/",paste0(SNPset,"_",datafilenames,"_Pvalue.Hyperparameters.pdf"))), width=8.3,height=11.7)
layout(matrix(c(1,1,2,3,4,4,5,6), 4, 2, byrow = TRUE))

# PVE
# ------------------------------------------------------------------------------
plot(hyp.params$pve, type="l", ylab="PVE", main="PVE - trace")
hist(hyp.params$pve, main="PVE - posterior distribution", xlab="PVE")
plot(density(hyp.params$pve), main="PVE - posterior distribution", xlab="PVE")
# ------------------------------------------------------------------------------

# PGE
# ------------------------------------------------------------------------------
plot(hyp.params$pge, type="l", ylab="PGE", main="PGE - trace")
hist(hyp.params$pge, main="PGE - posterior distribution", xlab="PGE")
plot(density(hyp.params$pge), main="PGE - posterior distribution", xlab="PGE")
# ------------------------------------------------------------------------------

# pi 
# ------------------------------------------------------------------------------
plot(hyp.params$pi, type="l", ylab="pi", main="pi")
hist(hyp.params$pi, main="pi", xlab="pi")
plot(density(hyp.params$pi), main="pi", xlab="pi")
# ------------------------------------------------------------------------------

# N_gamma #number of large effect snps
# ------------------------------------------------------------------------------
plot(hyp.params$n_gamma, type="l", ylab="n_gamma", main="n_gamma - trace")
hist(hyp.params$n_gamma, main="n_gamma - posterior distribution", xlab="n_gamma")
plot(density(hyp.params$n_gamma), main="n_gamma - posterior distribution", xlab="n_gamma")
# ------------------------------------------------------------------------------

# rho
# ------------------------------------------------------------------------------
plot(hyp.params$rho, type="l", ylab="rho", main="rho - trace")
hist(hyp.params$rho, main="rho - posterior distribution", xlab="rho")
plot(density(hyp.params$rho), main="rho - posterior distribution", xlab="rho")
# ------------------------------------------------------------------------------
#h
# ------------------------------------------------------------------------------
plot(hyp.params$h, type="l", ylab="h", main="h - trace")
hist(hyp.params$h, main="h - posterior distribution", xlab="h")
plot(density(hyp.params$h), main="h - posterior distribution", xlab="h")





dev.off()
# ==============================================================================


# Load parameters output
# ==============================================================================
params<-fread(paste0(paste0(SNPset,"_",datafilenames,"_Pvalue_Tables/",paste0(SNPset,"_",datafilenames,"_Pvalue.param.txt"))),header=T,sep="\t", data.table=F)
# ==============================================================================
# Get variants with sparse effect size on phenotypes 
# ==============================================================================
# add sparse effect size (= beta * gamma) to data frame
params["eff"]<-abs(params$beta*params$gamma)
#write params file 
write.table(params, file=paste0(paste0(SNPset,"_",datafilenames,"_Pvalue_Tables/",paste0(SNPset,"_",datafilenames,"_Pvalue.Sparseparameters.txt"))), quote=F, row.names=F, sep="\t")


# get variants with effect size > 0
params.effects<-params[params$eff>0,]

# show number of variants with measurable effect
nrow(params.effects)

# sort by descending effect size
params.effects.sort<-params.effects[order(-params.effects$eff),]

# show top 10 variants with highest effect
head(params.effects.sort, 10) 

# variants with the highest sparse effects
# ------------------------------------------------------------------------------
# top 1% variants (above 99% quantile)
top1<-params.effects.sort[params.effects.sort$eff>quantile(params.effects.sort$eff,0.99),]
# top 0.1% variants (above 99.9% quantile)
top01<-params.effects.sort[params.effects.sort$eff>quantile(params.effects.sort$eff,0.999),]
# top 0.01% variants (above 99.99% quantile)
top001<-params.effects.sort[params.effects.sort$eff>quantile(params.effects.sort$eff,0.9999),]
# ------------------------------------------------------------------------------

# write tables
write.table(top1, file=paste0(paste0(SNPset,"_",datafilenames,"_Pvalue_Tables/",paste0(SNPset,"_",datafilenames,"_Pvalue.top1eff.dsv"))), quote=F, row.names=F, sep="\t")
write.table(top01, file=paste0(paste0(SNPset,"_",datafilenames,"_Pvalue_Tables/",paste0(SNPset,"_",datafilenames,"_Pvalue.top0.1eff.dsv"))), quote=F, row.names=F, sep="\t")
write.table(top001, file=paste0(paste0(SNPset,"_",datafilenames,"_Pvalue_Tables/",paste0(SNPset,"_",datafilenames,"_Pvalue.top0.01eff.dsv"))), quote=F, row.names=F, sep="\t")


# ==============================================================================
# Get variants with high Posterior Inclusion Probability (PIP) == gamma
# ==============================================================================
# PIP is the frequency a variant is estimated to have a sparse effect in the MCMC

# sort variants by descending PIP
params.pipsort<-params[order(-params$gamma),]

# Show top 10 variants with highest PIP
head(params.pipsort,10)

# sets of variants above a certain threshold
# variants with effect in 1% MCMC samples or more
pip01<-params.pipsort[params.pipsort$gamma>=0.01,]
# variants with effect in 10% MCMC samples or more
pip10<-params.pipsort[params.pipsort$gamma>=0.10,]
# variants with effect in 25% MCMC samples or more
pip25<-params.pipsort[params.pipsort$gamma>=0.25,]
# variants with effect in 50% MCMC samples or more
pip50<-params.pipsort[params.pipsort$gamma>=0.50,]

# write tables
write.table(pip01, file=paste0(paste0(SNPset,"_",datafilenames,"_Pvalue_Tables/",paste0(SNPset,"_",datafilenames,"_Pvalue.pip01.dsv"))), quote=F, row.names=F, sep="\t")
write.table(pip10, file=paste0(paste0(SNPset,"_",datafilenames,"_Pvalue_Tables/",paste0(SNPset,"_",datafilenames,"_Pvalue.pip10.dsv"))), quote=F, row.names=F, sep="\t")
write.table(pip25, file=paste0(paste0(SNPset,"_",datafilenames,"_Pvalue_Tables/",paste0(SNPset,"_",datafilenames,"_Pvalue.pip25.dsv"))), quote=F, row.names=F, sep="\t")
write.table(pip50, file=paste0(paste0(SNPset,"_",datafilenames,"_Pvalue_Tables/",paste0(SNPset,"_",datafilenames,"_Pvalue.pip50.dsv"))), quote=F, row.names=F, sep="\t")
# ------------------------------------------------------------------------------
# ==============================================================================



#data is now ready to be used for bslmm
#create a directory to store bslmm outputs

dir.create(paste0(SNPset,"_",datafilenames,"_Beta_Plots/"))
dir.create(paste0(SNPset,"_",datafilenames,"_Beta_Tables/"))

#use plink to make bed files for each folder 
Beta_filePath<-"data/EffectsizeSNPS/"
BetaSNPSET<-paste0(SNPset,"_",datafilenames,"_Beta")


#USE GEMMA to run a bayesian sparse linear mixed model 
system(paste("Software/gemma -bslmm 1 -bfile ",paste0(Beta_filePath, BetaSNPSET)," -w 250000"," -rpace 100 -wpace 1000"," -k ",paste0(datafilepath,SNPset,".cXX.txt")," -outdir", paste0(SNPset,"_",datafilenames,"_Beta_Tables/"), "-o " ,paste0(SNPset,"_",datafilenames,"_Beta")))





# Load hyperparameter file for effectsize
# ==============================================================================
hyp.params<-read.table(paste0(paste0(SNPset,"_",datafilenames,"_Beta_Tables/",paste0(SNPset,"_",datafilenames,"_Beta.hyp.txt"))),header = T)

# ==============================================================================
# Get mean, median, and 95% ETPI of hyperparameters
# ==============================================================================
# pve -> proportion of phenotypic variance explained by the genotypes
pve<-c("PVE", mean(hyp.params$pve),quantile(hyp.params$pve, probs=c(0.5,0.025,0.975)))

# pge -> proportion of genetic variance explained by major effect loci
pge<-c("PGE",mean(hyp.params$pge),quantile(hyp.params$pge, probs=c(0.5,0.025,0.975)))

# pi -> proportion of variants with non-zero effects
pi<-c("pi",mean(hyp.params$pi),quantile(hyp.params$pi, probs=c(0.5,0.025,0.975)))

# n.gamma -> number of variants with major effect
n.gamma<-c("n.gamma",mean(hyp.params$n_gamma),quantile(hyp.params$n_gamma, probs=c(0.5,0.025,0.975)))
# rho -> approximation to proportion of genetic variance explained by variants with major effects,
#rho close to zero trait is highly polygenic
rho<-c("rho",mean(hyp.params$rho),quantile(hyp.params$rho, probs=c(0.5,0.025,0.975)))


# ==============================================================================

# get table of hyperparameters
# ==============================================================================
hyp.params.table<-as.data.frame(rbind(pve,pge,pi,n.gamma,rho),row.names=F)
colnames(hyp.params.table)<-c("hyperparam", "mean","median","2.5%", "97.5%")
# show table
View(hyp.params.table)
write.table(hyp.params.table,paste0(paste0(SNPset,"_",datafilenames,"_Beta_Tables/",paste0(SNPset,"_",datafilenames,"_Beta.Hyperparametertable.txt"))))

#
# plot traces and distributions of hyperparameters
# ==============================================================================
pdf(file=paste0(paste0(SNPset,"_",datafilenames,"_Beta_Tables/",paste0(SNPset,"_",datafilenames,"_Beta.Hyperparameters.pdf"))), width=8.3,height=11.7)
layout(matrix(c(1,1,2,3,4,4,5,6), 4, 2, byrow = TRUE))

# PVE
# ------------------------------------------------------------------------------
plot(hyp.params$pve, type="l", ylab="PVE", main="PVE - trace")
hist(hyp.params$pve, main="PVE - posterior distribution", xlab="PVE")
plot(density(hyp.params$pve), main="PVE - posterior distribution", xlab="PVE")
# ------------------------------------------------------------------------------

# PGE
# ------------------------------------------------------------------------------
plot(hyp.params$pge, type="l", ylab="PGE", main="PGE - trace")
hist(hyp.params$pge, main="PGE - posterior distribution", xlab="PGE")
plot(density(hyp.params$pge), main="PGE - posterior distribution", xlab="PGE")
# ------------------------------------------------------------------------------

# pi 
# ------------------------------------------------------------------------------
plot(hyp.params$pi, type="l", ylab="pi", main="pi")
hist(hyp.params$pi, main="pi", xlab="pi")
plot(density(hyp.params$pi), main="pi", xlab="pi")
# ------------------------------------------------------------------------------

# N_gamma #number of large effect snps
# ------------------------------------------------------------------------------
plot(hyp.params$n_gamma, type="l", ylab="n_gamma", main="n_gamma - trace")
hist(hyp.params$n_gamma, main="n_gamma - posterior distribution", xlab="n_gamma")
plot(density(hyp.params$n_gamma), main="n_gamma - posterior distribution", xlab="n_gamma")
# ------------------------------------------------------------------------------

# rho
# ------------------------------------------------------------------------------
plot(hyp.params$rho, type="l", ylab="rho", main="rho - trace")
hist(hyp.params$rho, main="rho - posterior distribution", xlab="rho")
plot(density(hyp.params$rho), main="rho - posterior distribution", xlab="rho")
# ------------------------------------------------------------------------------
#h
# ------------------------------------------------------------------------------
plot(hyp.params$h, type="l", ylab="h", main="h - trace")
hist(hyp.params$h, main="h - posterior distribution", xlab="h")
plot(density(hyp.params$h), main="h - posterior distribution", xlab="h")





dev.off()
# ==============================================================================

# Load parameters output
# ==============================================================================
params<-fread(paste0(paste0(SNPset,"_",datafilenames,"_Beta_Tables/",paste0(SNPset,"_",datafilenames,"_Beta.param.txt"))),header=T,sep="\t", data.table=F)
# ==============================================================================
# Get variants with sparse effect size on phenotypes 
# ==============================================================================
# add sparse effect size (= beta * gamma) to data frame
params["eff"]<-abs(params$beta*params$gamma)
#safe file for plotting
write.table(params,paste0(paste0(SNPset,"_",datafilenames,"_Beta_Tables/",paste0(SNPset,"_",datafilenames,"_Beta.Sparseparameters.txt"))), quote=F, row.names=F, sep="\t")

# get variants with effect size > 0
params.effects<-params[params$eff>0,]

# show number of variants with measurable effect
nrow(params.effects)

# sort by descending effect size
params.effects.sort<-params.effects[order(-params.effects$eff),]

# show top 10 variants with highest effect
head(params.effects.sort, 10) 

# variants with the highest sparse effects
# ------------------------------------------------------------------------------
# top 1% variants (above 99% quantile)
top1<-params.effects.sort[params.effects.sort$eff>quantile(params.effects.sort$eff,0.99),]
# top 0.1% variants (above 99.9% quantile)
top01<-params.effects.sort[params.effects.sort$eff>quantile(params.effects.sort$eff,0.999),]
# top 0.01% variants (above 99.99% quantile)
top001<-params.effects.sort[params.effects.sort$eff>quantile(params.effects.sort$eff,0.9999),]
# ------------------------------------------------------------------------------

# write tables
write.table(top1, file=paste0(paste0(SNPset,"_",datafilenames,"_Beta_Tables/",paste0(SNPset,"_",datafilenames,"_Beta.top1eff.dsv"))), quote=F, row.names=F, sep="\t")
write.table(top01, file=paste0(paste0(SNPset,"_",datafilenames,"_Beta_Tables/",paste0(SNPset,"_",datafilenames,"_Beta.top0.1eff.dsv"))), quote=F, row.names=F, sep="\t")
write.table(top001, file=paste0(paste0(SNPset,"_",datafilenames,"_Beta_Tables/",paste0(SNPset,"_",datafilenames,"_Beta.top0.01eff.dsv"))), quote=F, row.names=F, sep="\t")


# ==============================================================================
# Get variants with high Posterior Inclusion Probability (PIP) == gamma
# ==============================================================================
# PIP is the frequency a variant is estimated to have a sparse effect in the MCMC

# sort variants by descending PIP
params.pipsort<-params[order(-params$gamma),]

# Show top 10 variants with highest PIP
head(params.pipsort,10)

# sets of variants above a certain threshold
# variants with effect in 1% MCMC samples or more
pip01<-params.pipsort[params.pipsort$gamma>=0.01,]
# variants with effect in 10% MCMC samples or more
pip10<-params.pipsort[params.pipsort$gamma>=0.10,]
# variants with effect in 25% MCMC samples or more
pip25<-params.pipsort[params.pipsort$gamma>=0.25,]
# variants with effect in 50% MCMC samples or more
pip50<-params.pipsort[params.pipsort$gamma>=0.50,]

# write tables
write.table(pip01, file=paste0(paste0(SNPset,"_",datafilenames,"_Beta_Tables/",paste0(SNPset,"_",datafilenames,"_Beta.pip01.dsv"))), quote=F, row.names=F, sep="\t")
write.table(pip10, file=paste0(paste0(SNPset,"_",datafilenames,"_Beta_Tables/",paste0(SNPset,"_",datafilenames,"_Beta.pip10.dsv"))), quote=F, row.names=F, sep="\t")
write.table(pip25, file=paste0(paste0(SNPset,"_",datafilenames,"_Beta_Tables/",paste0(SNPset,"_",datafilenames,"_Beta.pip25.dsv"))), quote=F, row.names=F, sep="\t")
write.table(pip50, file=paste0(paste0(SNPset,"_",datafilenames,"_Beta_Tables/",paste0(SNPset,"_",datafilenames,"_Beta.pip50.dsv"))), quote=F, row.names=F, sep="\t")
# ------------------------------------------------------------------------------
# ==============================================================================


#Combination data is now ready to be used for bslmm
#create a directory to store bslmm outputs

dir.create(paste0(SNPset,"_",datafilenames,"_Combo_Plots/"))
dir.create(paste0(SNPset,"_",datafilenames,"_Combo_Tables/"))

#use plink to make bed files for each folder 
Combo_filePath<-"data/ComboSNPS/"
ComboSNPSET<-paste0(SNPset,"_",datafilenames,"_PvalueBeta")


#USE GEMMA to run a bayesian sparse linear mixed model 
system(paste("Software/gemma -bslmm 1 -bfile ",paste0(Combo_filePath, ComboSNPSET)," -w 250000"," -rpace 100 -wpace 1000"," -k ",paste0(datafilepath,SNPset,".cXX.txt")," -outdir", paste0(SNPset,"_",datafilenames,"_Combo_Tables/"), "-o " ,paste0(SNPset,"_",datafilenames,"_Combo")))





# Load hyperparameter file for effectsize
# ==============================================================================
hyp.params<-read.table(paste0(paste0(SNPset,"_",datafilenames,"_Combo_Tables/",paste0(SNPset,"_",datafilenames,"_Combo.hyp.txt"))),header = T)

# ==============================================================================
# Get mean, median, and 95% ETPI of hyperparameters
# ==============================================================================
# pve -> proportion of phenotypic variance explained by the genotypes
pve<-c("PVE", mean(hyp.params$pve),quantile(hyp.params$pve, probs=c(0.5,0.025,0.975)))

# pge -> proportion of genetic variance explained by major effect loci
pge<-c("PGE",mean(hyp.params$pge),quantile(hyp.params$pge, probs=c(0.5,0.025,0.975)))

# pi -> proportion of variants with non-zero effects
pi<-c("pi",mean(hyp.params$pi),quantile(hyp.params$pi, probs=c(0.5,0.025,0.975)))

# n.gamma -> number of variants with major effect
n.gamma<-c("n.gamma",mean(hyp.params$n_gamma),quantile(hyp.params$n_gamma, probs=c(0.5,0.025,0.975)))
# rho -> approximation to proportion of genetic variance explained by variants with major effects,
#rho close to zero trait is highly polygenic
rho<-c("rho",mean(hyp.params$rho),quantile(hyp.params$rho, probs=c(0.5,0.025,0.975)))


# ==============================================================================

# get table of hyperparameters
# ==============================================================================
hyp.params.table<-as.data.frame(rbind(pve,pge,pi,n.gamma,rho),row.names=F)
colnames(hyp.params.table)<-c("hyperparam", "mean","median","2.5%", "97.5%")
# show table
View(hyp.params.table)
write.table(hyp.params.table,paste0(paste0(SNPset,"_",datafilenames,"_Combo_Tables/",paste0(SNPset,"_",datafilenames,"_Combo.Hyperparametertable.txt"))))

#
# plot traces and distributions of hyperparameters
# ==============================================================================
pdf(file=paste0(paste0(SNPset,"_",datafilenames,"_Combo_Tables/",paste0(SNPset,"_",datafilenames,"_Combo.Hyperparameters.pdf"))), width=8.3,height=11.7)
layout(matrix(c(1,1,2,3,4,4,5,6), 4, 2, byrow = TRUE))

# PVE
# ------------------------------------------------------------------------------
plot(hyp.params$pve, type="l", ylab="PVE", main="PVE - trace")
hist(hyp.params$pve, main="PVE - posterior distribution", xlab="PVE")
plot(density(hyp.params$pve), main="PVE - posterior distribution", xlab="PVE")
# ------------------------------------------------------------------------------

# PGE
# ------------------------------------------------------------------------------
plot(hyp.params$pge, type="l", ylab="PGE", main="PGE - trace")
hist(hyp.params$pge, main="PGE - posterior distribution", xlab="PGE")
plot(density(hyp.params$pge), main="PGE - posterior distribution", xlab="PGE")
# ------------------------------------------------------------------------------

# pi 
# ------------------------------------------------------------------------------
plot(hyp.params$pi, type="l", ylab="pi", main="pi")
hist(hyp.params$pi, main="pi", xlab="pi")
plot(density(hyp.params$pi), main="pi", xlab="pi")
# ------------------------------------------------------------------------------

# N_gamma #number of large effect snps
# ------------------------------------------------------------------------------
plot(hyp.params$n_gamma, type="l", ylab="n_gamma", main="n_gamma - trace")
hist(hyp.params$n_gamma, main="n_gamma - posterior distribution", xlab="n_gamma")
plot(density(hyp.params$n_gamma), main="n_gamma - posterior distribution", xlab="n_gamma")
# ------------------------------------------------------------------------------

# rho
# ------------------------------------------------------------------------------
plot(hyp.params$rho, type="l", ylab="rho", main="rho - trace")
hist(hyp.params$rho, main="rho - posterior distribution", xlab="rho")
plot(density(hyp.params$rho), main="rho - posterior distribution", xlab="rho")
# ------------------------------------------------------------------------------
#h
# ------------------------------------------------------------------------------
plot(hyp.params$h, type="l", ylab="h", main="h - trace")
hist(hyp.params$h, main="h - posterior distribution", xlab="h")
plot(density(hyp.params$h), main="h - posterior distribution", xlab="h")





dev.off()
# ==============================================================================


# Load parameters output
# ==============================================================================
params<-fread(paste0(paste0(SNPset,"_",datafilenames,"_Combo_Tables/",paste0(SNPset,"_",datafilenames,"_Combo.param.txt"))),header=T,sep="\t", data.table=F)
# ==============================================================================
# Get variants with sparse effect size on phenotypes 
# ==============================================================================
# add sparse effect size (= beta * gamma) to data frame
params["eff"]<-abs(params$beta*params$gamma)
#safe file for plotting
write.table(params,paste0(paste0(SNPset,"_",datafilenames,"_Combo_Tables/",paste0(SNPset,"_",datafilenames,"_Combo.Sparseparameters.txt"))), quote=F, row.names=F, sep="\t")

# get variants with effect size > 0
params.effects<-params[params$eff>0,]

# show number of variants with measurable effect
nrow(params.effects)

# sort by descending effect size
params.effects.sort<-params.effects[order(-params.effects$eff),]

# show top 10 variants with highest effect
head(params.effects.sort, 10) 

# variants with the highest sparse effects
# ------------------------------------------------------------------------------
# top 1% variants (above 99% quantile)
top1<-params.effects.sort[params.effects.sort$eff>quantile(params.effects.sort$eff,0.99),]
# top 0.1% variants (above 99.9% quantile)
top01<-params.effects.sort[params.effects.sort$eff>quantile(params.effects.sort$eff,0.999),]
# top 0.01% variants (above 99.99% quantile)
top001<-params.effects.sort[params.effects.sort$eff>quantile(params.effects.sort$eff,0.9999),]
# ------------------------------------------------------------------------------

# write tables
write.table(top1, file=paste0(paste0(SNPset,"_",datafilenames,"_Combo_Tables/",paste0(SNPset,"_",datafilenames,"_Combo.top1eff.dsv"))), quote=F, row.names=F, sep="\t")
write.table(top01, file=paste0(paste0(SNPset,"_",datafilenames,"_Combo_Tables/",paste0(SNPset,"_",datafilenames,"_Combo.top0.1eff.dsv"))), quote=F, row.names=F, sep="\t")
write.table(top001, file=paste0(paste0(SNPset,"_",datafilenames,"_Combo_Tables/",paste0(SNPset,"_",datafilenames,"_Combo.top0.01eff.dsv"))), quote=F, row.names=F, sep="\t")


# ==============================================================================
# Get variants with high Posterior Inclusion Probability (PIP) == gamma
# ==============================================================================
# PIP is the frequency a variant is estimated to have a sparse effect in the MCMC

# sort variants by descending PIP
params.pipsort<-params[order(-params$gamma),]

# Show top 10 variants with highest PIP
head(params.pipsort,10)

# sets of variants above a certain threshold
# variants with effect in 1% MCMC samples or more
pip01<-params.pipsort[params.pipsort$gamma>=0.01,]
# variants with effect in 10% MCMC samples or more
pip10<-params.pipsort[params.pipsort$gamma>=0.10,]
# variants with effect in 25% MCMC samples or more
pip25<-params.pipsort[params.pipsort$gamma>=0.25,]
# variants with effect in 50% MCMC samples or more
pip50<-params.pipsort[params.pipsort$gamma>=0.50,]

# write tables
write.table(pip01, file=paste0(paste0(SNPset,"_",datafilenames,"_Combo_Tables/",paste0(SNPset,"_",datafilenames,"_Combo.pip01.dsv"))), quote=F, row.names=F, sep="\t")
write.table(pip10, file=paste0(paste0(SNPset,"_",datafilenames,"_Combo_Tables/",paste0(SNPset,"_",datafilenames,"_Combo.pip10.dsv"))), quote=F, row.names=F, sep="\t")
write.table(pip25, file=paste0(paste0(SNPset,"_",datafilenames,"_Combo_Tables/",paste0(SNPset,"_",datafilenames,"_Combo.pip25.dsv"))), quote=F, row.names=F, sep="\t")
write.table(pip50, file=paste0(paste0(SNPset,"_",datafilenames,"_Combo_Tables/",paste0(SNPset,"_",datafilenames,"_Combo.pip50.dsv"))), quote=F, row.names=F, sep="\t")
# ------------------------------------------------------------------------------
# ==============================================================================



#########analyzing KEGG SNPS
#create appropriate directories
dir.create(paste0(KEGGSNPS,"_Plots/"))
dir.create(paste0(KEGGSNPS,"_Tables/"))
dir.create((paste0(KEGGSNPS,"_Tables/Assoc_files/")))




#if the phenotype does not exist go to the next one

trait.data<-pheno.data



#rewrite the fam file to include the new phenotypic information 
fam.file<-fread(paste0(KEGG.datafilepath ,KEGGSNPS,".fam"))
fam.file$V6<-NULL
fam.file <- merge(fam.file,trait.data,by.x="V1",by.y="Line",all.x=T)
write.table(file=paste0(KEGG.datafilepath ,KEGGSNPS,".fam"),fam.file,col.names=F, row.names=F, quote =F)





#using GEMMA to run a bayesian sparse linear mixed model on the KEGG SNP set 
system(paste("Software/gemma -bslmm 1 -bfile ",paste0(KEGG.datafilepath,KEGGSNPS)," -w 250000"," -rpace 100 -wpace 1000"," -k ",paste0(KEGG.datafilepath,KEGGSNPS,".cXX.txt")," -outdir", paste0(KEGGSNPS,"_Tables/"), "-o " ,paste0(KEGGSNPS,"_",datafilenames)))



# Examine Hyperparameters for bayesian sparse linear mixed model of KEGG SNPS
# Load hyperparameter file for p-value
# ==============================================================================
hyp.params<-read.table(paste0(paste0(KEGGSNPS,"_Tables/",paste0(KEGGSNPS,"_",datafilenames,".hyp.txt"))),header = T)

# ==============================================================================
# Get mean, median, and 95% ETPI of hyperparameters
# ==============================================================================
# pve -> proportion of phenotypic variance explained by the genotypes
pve<-c("PVE", mean(hyp.params$pve),quantile(hyp.params$pve, probs=c(0.5,0.025,0.975)))

# pge -> proportion of genetic variance explained by major effect loci
pge<-c("PGE",mean(hyp.params$pge),quantile(hyp.params$pge, probs=c(0.5,0.025,0.975)))

# pi -> proportion of variants with non-zero effects
pi<-c("pi",mean(hyp.params$pi),quantile(hyp.params$pi, probs=c(0.5,0.025,0.975)))

# n.gamma -> number of variants with major effect
n.gamma<-c("n.gamma",mean(hyp.params$n_gamma),quantile(hyp.params$n_gamma, probs=c(0.5,0.025,0.975)))
# rho -> approximation to proportion of genetic variance explained by variants with major effects,
#rho close to zero trait is highly polygenic
rho<-c("rho",mean(hyp.params$rho),quantile(hyp.params$rho, probs=c(0.5,0.025,0.975)))


# ==============================================================================

# get table of hyperparameters
# ==============================================================================
hyp.params.table<-as.data.frame(rbind(pve,pge,pi,n.gamma,rho),row.names=F)
colnames(hyp.params.table)<-c("hyperparam", "mean","median","2.5%", "97.5%")
# show table
View(hyp.params.table)
write.table(hyp.params.table,paste0(paste0(KEGGSNPS,"_Tables/",paste0(KEGGSNPS,"_",datafilenames,".Hyperparametertable.txt"))))
paste0(paste0(KEGGSNPS,"_",datafilenames,"_Pvalue_Tables/",paste0(KEGGSNPS,"_",datafilenames,"_Pvalue.Hyperparametertable.txt")))
#
# plot traces and distributions of hyperparameters
# ==============================================================================
pdf(file=paste0(paste0(KEGGSNPS,"_Tables/",paste0(KEGGSNPS,"_",datafilenames,".Hyperparameters.pdf"))), width=8.3,height=11.7)
layout(matrix(c(1,1,2,3,4,4,5,6), 4, 2, byrow = TRUE))

# PVE
# ------------------------------------------------------------------------------
plot(hyp.params$pve, type="l", ylab="PVE", main="PVE - trace")
hist(hyp.params$pve, main="PVE - posterior distribution", xlab="PVE")
plot(density(hyp.params$pve), main="PVE - posterior distribution", xlab="PVE")
# ------------------------------------------------------------------------------

# PGE
# ------------------------------------------------------------------------------
plot(hyp.params$pge, type="l", ylab="PGE", main="PGE - trace")
hist(hyp.params$pge, main="PGE - posterior distribution", xlab="PGE")
plot(density(hyp.params$pge), main="PGE - posterior distribution", xlab="PGE")
# ------------------------------------------------------------------------------

# pi 
# ------------------------------------------------------------------------------
plot(hyp.params$pi, type="l", ylab="pi", main="pi")
hist(hyp.params$pi, main="pi", xlab="pi")
plot(density(hyp.params$pi), main="pi", xlab="pi")
# ------------------------------------------------------------------------------

# N_gamma #number of large effect snps
# ------------------------------------------------------------------------------
plot(hyp.params$n_gamma, type="l", ylab="n_gamma", main="n_gamma - trace")
hist(hyp.params$n_gamma, main="n_gamma - posterior distribution", xlab="n_gamma")
plot(density(hyp.params$n_gamma), main="n_gamma - posterior distribution", xlab="n_gamma")
# ------------------------------------------------------------------------------

# rho
# ------------------------------------------------------------------------------
plot(hyp.params$rho, type="l", ylab="rho", main="rho - trace")
hist(hyp.params$rho, main="rho - posterior distribution", xlab="rho")
plot(density(hyp.params$rho), main="rho - posterior distribution", xlab="rho")
# ------------------------------------------------------------------------------
#h
# ------------------------------------------------------------------------------
plot(hyp.params$h, type="l", ylab="h", main="h - trace")
hist(hyp.params$h, main="h - posterior distribution", xlab="h")
plot(density(hyp.params$h), main="h - posterior distribution", xlab="h")





dev.off()
# ==============================================================================


# Load parameters output
# ==============================================================================
params<-fread(paste0(paste0(KEGGSNPS,"_Tables/",paste0(KEGGSNPS,"_",datafilenames,".param.txt"))),header=T,sep="\t", data.table=F)
# ==============================================================================
# Get variants with sparse effect size on phenotypes 
# ==============================================================================
# add sparse effect size (= beta * gamma) to data frame
params["eff"]<-abs(params$beta*params$gamma)

# get variants with effect size > 0
params.effects<-params[params$eff>0,]

# show number of variants with measurable effect
nrow(params.effects)

# sort by descending effect size
params.effects.sort<-params.effects[order(-params.effects$eff),]

# show top 10 variants with highest effect
head(params.effects.sort, 10) 

#write file with sparse effect calculated 

write.table(params,file=paste0(KEGGSNPS,"_Tables/",KEGGSNPS,"_",datafilenames,".Sparseparameters.txt"), quote=F, row.names=F, sep="\t")
# variants with the highest sparse effects
# ------------------------------------------------------------------------------
# top 1% variants (above 99% quantile)
top1<-params.effects.sort[params.effects.sort$eff>quantile(params.effects.sort$eff,0.99),]
# top 0.1% variants (above 99.9% quantile)
top01<-params.effects.sort[params.effects.sort$eff>quantile(params.effects.sort$eff,0.999),]
# top 0.01% variants (above 99.99% quantile)
top001<-params.effects.sort[params.effects.sort$eff>quantile(params.effects.sort$eff,0.9999),]
# ------------------------------------------------------------------------------

# write tables
write.table(top1, file=paste0(paste0(KEGGSNPS,"_Tables/",paste0(KEGGSNPS,"_",datafilenames,".top1eff.dsv"))), quote=F, row.names=F, sep="\t")
write.table(top01, file=paste0(paste0(KEGGSNPS,"_Tables/",paste0(KEGGSNPS,"_",datafilenames,".top0.1eff.dsv"))), quote=F, row.names=F, sep="\t")
write.table(top001, file=paste0(paste0(KEGGSNPS,"_Tables/",paste0(KEGGSNPS,"_",datafilenames,".top0.01eff.dsv"))), quote=F, row.names=F, sep="\t")


# ==============================================================================
# Get variants with high Posterior Inclusion Probability (PIP) == gamma
# ==============================================================================
# PIP is the frequency a variant is estimated to have a sparse effect in the MCMC

# sort variants by descending PIP
params.pipsort<-params[order(-params$gamma),]

# Show top 10 variants with highest PIP
head(params.pipsort,10)

# sets of variants above a certain threshold
# variants with effect in 1% MCMC samples or more
pip01<-params.pipsort[params.pipsort$gamma>=0.01,]
# variants with effect in 10% MCMC samples or more
pip10<-params.pipsort[params.pipsort$gamma>=0.10,]
# variants with effect in 25% MCMC samples or more
pip25<-params.pipsort[params.pipsort$gamma>=0.25,]
# variants with effect in 50% MCMC samples or more
pip50<-params.pipsort[params.pipsort$gamma>=0.50,]

# write tables
write.table(pip01, file=paste0(paste0(KEGGSNPS,"_Tables/",paste0(KEGGSNPS,"_",datafilenames,".pip01.dsv"))), quote=F, row.names=F, sep="\t")
write.table(pip10, file=paste0(paste0(KEGGSNPS,"_Tables/",paste0(KEGGSNPS,"_",datafilenames,".pip10.dsv"))), quote=F, row.names=F, sep="\t")
write.table(pip25, file=paste0(paste0(KEGGSNPS,"_Tables/",paste0(KEGGSNPS,"_",datafilenames,".pip25.dsv"))), quote=F, row.names=F, sep="\t")
write.table(pip50, file=paste0(paste0(KEGGSNPS,"_Tables/",paste0(KEGGSNPS,"_",datafilenames,".pip50.dsv"))), quote=F, row.names=F, sep="\t")
# ------------------------------------------------------------------------------
# ==============================================================================





#using GEMMA to run a  linear mixed model 
system(paste("Software/gemma -bfile ",paste0(KEGG.datafilepath,KEGGSNPS)," -k ",paste0(KEGG.datafilepath,KEGGSNPS,".cXX.txt"), "-c ",paste0(KEGG.datafilepath,KEGGSNPS,".PCA_EV"), "-lmm 1 -outdir ",paste0(KEGGSNPS,"_Tables/Assoc_files/"), "-o ",paste0(KEGGSNPS,"_",datafilenames)))




#end of script












#read in outputs to make 3 datasets top 1% of SNPs based on p-value ,effect size, and a combination of both newly formed datasets


#read in files tped, map, & lmm output
#we are making an assumption that the SNPlist and SNPvariants are in the same order
SNPMap<-read.table(paste0(KEGG.datafilepath,KEGGSNPS,".map"))
SNPTped<- read.table(paste0(KEGG.datafilepath,KEGGSNPS,".tped"))



#read in the lmm output 
LMMouput.file<-fread(paste0(KEGGSNPS,"_Tables/Assoc_files/",KEGGSNPS,"_",datafilenames,".assoc.txt"))









