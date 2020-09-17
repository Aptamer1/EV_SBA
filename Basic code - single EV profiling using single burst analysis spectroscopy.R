####Direct quality control and multidimensional phenotyping of  single extracellular vesicles using burst analysis spectroscopy.####

###Load libraries (to install packages use the install.packages command)
library(ggplot2)
library(dplyr)
library(magrittr)
library(MASS)

###
###All valus that have to be inserted are indicated by 'INSERT ... VALUE HERE'

###Read in the burst data
df <-read.table("INSERT FILE NAME HERE.csv",sep=",",header=T)

###Make a copy
dfcopy <- df

########
##2C cross-talk correction
########
#Channel DD 
df$Number.of.Photons..DD. <- dfcopy$Number.of.Photons..DD.  - INSERT LS VALUE HERE * dfcopy$Number.of.Photons..AA. 

#Channel AA
df$Number.of.Photons..AA. <- dfcopy$Number.of.Photons..AA.  - INSERT LS VALUE HERE * dfcopy$Number.of.Photons..DD. 

##set all negative values to 1
df[c(14, 16, 18, 20)] <- replace(df[c(14, 16, 18, 20)], df[c(14, 16, 18, 20)] < 1, 1)

#######
##2C thresholding
#######

#Get rid of inexplicably high values
df <- subset (df, df$Number.of.Photons..DD. <  2000)
df <- subset (df, df$Number.of.Photons..AA. <  2000)

###Define the threshold values for each channel -- Any signal above this value will be defined as a burst
thresDD <- mean(df$Number.of.Photons..DD.) + 3*sd(df$Number.of.Photons..DD.)
thresAA <- mean(df$Number.of.Photons..AA.) + 3*sd(df$Number.of.Photons..AA.)

###Apply the threshold
df <- subset (df, df$Number.of.Photons..DD. >  thresDD | 
                df$Number.of.Photons..AA. >  thresAA)
#####
##Example code to define marker positivity -- adjust the markers to the markers that have been used in your experiment
####
df['name'] <- with(df, ifelse(Number.of.Photons..DD.  > thresDD &
                                Number.of.Photons..AA.  > thresAA, 'CD9 ICAM pos',
                              ifelse(Number.of.Photons..DD.  > thresDD & 
                                       Number.of.Photons..AA.  < thresAA , 'CD9 pos',
                                     ifelse(Number.of.Photons..DD.  < thresDD & 
                                              Number.of.Photons..AA.  > thresAA , 'ICAM pos', 'Und'))))
df$name <- as.factor(df$name)

######
####Size calculation and plotting of individual bursts
######
df['size']= ((((df$Duration..ms.*10^-3)*8*293.15*1.3806*10^-23)/(6*3.14*(8.9*10^-4)*(INSERT OMEGA r VALUE HERE *10^-6)^2))*10^9)

tiff("INSERT FILE NAME HERE.tiff", units="in", width=4, height=4, res=300)
ggplot(df, aes(x=size)) + geom_histogram(bins = 50) +
  scale_x_continuous(limits = c(0, 500))+
  labs(
    x = "Size (nm)",
    y = "Counts")+
  theme_bw() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=15))
dev.off()

###
##Time trace plotting
###

tiff("INSERT FILE NAME HERE.tiff", units="in", width=4, height=4, res=300)
ggplot() +
  geom_line(data=df, aes(x=Mean.Macrotime..ms., y=Number.of.Photons..DD.), color = 'red') +
  geom_line(data=df, aes(x=Mean.Macrotime..ms., y=Number.of.Photons..AA.), color = 'blue') +
  labs(
    x = "Time (s)",
    y = "Counted photons")+
  geom_hline(yintercept= thresDD, linetype="dashed", color = "red")+
  geom_hline(yintercept=thresAA, linetype="dashed", color = "blue")+
  theme_bw() +
  theme(text = element_text(size=15))
dev.off()


####
##Plotting of individual bursts in 2 channels
####

tiff("INSERT FILE NAME HERE.tiff", units="in", width=4, height=4, res=300)
ggplot() +
  geom_point(data=df, aes(x=log(Number.of.Photons..DD.), y=log(Number.of.Photons..AA.)), size = 0.5) +
  labs(
    x = "INSERT MARKER NAME HERE (ln counted photons)",
    y = "INSERT MARKER NAME HERE (ln counted photons)")+
  geom_hline(yintercept= log(thresAA), linetype="dashed", color = "red")+
  geom_vline(xintercept=log(thresDD), linetype="dashed", color = "red")+
  theme_bw() +
  theme(axis.text=element_text(size=13),
        axis.title=element_text(size=14))
dev.off()

########
##3C cross-talk correction
########

#Channel BB

df$Number.of.Photons..BB. <- dfcopy$Number.of.Photons..BB.  - INSERT LS VALUE HERE * dfcopy$Number.of.Photons..BG. - INSERT LS VALUE HERE * dfcopy$Number.of.Photons..BR.

#Channel BG

df$Number.of.Photons..BG. <- dfcopy$Number.of.Photons..BG.  - INSERT LS VALUE HERE * dfcopy$Number.of.Photons..BB. - INSERT LS VALUE HERE * dfcopy$Number.of.Photons..BR.

#Channel BR

df$Number.of.Photons..BR. <- dfcopy$Number.of.Photons..BR. - INSERT LS VALUE HERE * dfcopy$Number.of.Photons..BB. - INSERT LS VALUE HERE * dfcopy$Number.of.Photons..BG.

##set all negative values to 1
df[c(41, 42, 43, 59,60,61)] <- replace(df[c(41, 42, 43, 59,60,61)], df[c(41, 42, 43, 59,60,61)] < 1, 1)

#Get rid of inexplicably high values
df <- subset (df, df$Number.of.Photons..BB. <  2000)
df <- subset (df, df$Number.of.Photons..BG. <  2000)
df <- subset (df, df$Number.of.Photons..BR. <  2000)

###Define the threshold values for each channel -- Any signal above this value will be defined as a burst
thresBB <- mean(df$Number.of.Photons..BB.) + 3*sd(df$Number.of.Photons..BB.)
thresBG <- mean(df$Number.of.Photons..BG.) + 3*sd(df$Number.of.Photons..BG.)
thresBR <- mean(df$Number.of.Photons..BR.) + 3*sd(df$Number.of.Photons..BR.)

###Apply the threshold
df <- subset (df, df$Number.of.Photons..BB. >  thresBB | 
                df$Number.of.Photons..BG. >  thresBG |
                df$Number.of.Photons..BR. >  thresBR)

#####
##Example code to define marker positivity -- adjust the markers to the markers that have been used in your experiment
####
df['name'] <- with(df, ifelse(Number.of.Photons..BB.  > thresBB &
                                Number.of.Photons..BG.  > thresBG &
                                Number.of.Photons..BR.  > thresBR, 'CD9 CD63 ICAM pos',
                              ifelse(Number.of.Photons..BB.  > thresBB & 
                                       Number.of.Photons..BG.  > thresBG &
                                       Number.of.Photons..BR.  < thresBR, 'CD9 CD63 pos',
                                     ifelse(Number.of.Photons..BB.  > thresBB & 
                                              Number.of.Photons..BG.  < thresBG &
                                              Number.of.Photons..BR.  > thresBR, 'CD9 ICAM pos',
                                            ifelse(Number.of.Photons..BB.  < thresBB & 
                                                     Number.of.Photons..BG.  > thresBG &
                                                     Number.of.Photons..BR.  > thresBR, 'CD63 ICAM pos',
                                                   ifelse(Number.of.Photons..BB.  > thresBB & 
                                                            Number.of.Photons..BG.  < thresBG &
                                                            Number.of.Photons..BR.  < thresBR, 'CD9 pos',
                                                          ifelse(Number.of.Photons..BB.  < thresBB & 
                                                                   Number.of.Photons..BG.  > thresBG &
                                                                   Number.of.Photons..BR.  < thresBR, 'CD63 pos',
                                                                 ifelse(Number.of.Photons..BB.  < thresBB & 
                                                                          Number.of.Photons..BG.  < thresBG &
                                                                          Number.of.Photons..BR.  > thresBR, 'ICAM pos', 'Und'))))))))
df$name <- as.factor(df$name)

#################
###Size calculation
#########
df['size']= ((((df$Duration..ms.*10^-3)*8*293.15*1.3806*10^-23)/(6*3.14*(8.9*10^-4)*(INSERT OMEGA r VALUE HERE *10^-6)^2))*10^9)

tiff("INSERT FILE NAME HERE.tiff", units="in", width=4, height=4, res=300)
ggplot(df, aes(x=size)) + geom_histogram(bins = 50) +
  scale_x_continuous(limits = c(0, 500))+
  labs(
    x = "Size (nm)",
    y = "Counts")+
  theme_bw() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=15))
dev.off()

#######
##Time trace plotting
######

tiff("INSERT FILE NAME HERE.tiff", units="in", width=4, height=4, res=300)
ggplot() +
  geom_line(data=df, aes(x=Mean.Macrotime..s., y=Number.of.Photons..BB.), color = 'red') +
  geom_line(data=df, aes(x=Mean.Macrotime..s., y=Number.of.Photons..BG.), color = 'blue') +
  geom_line(data=df, aes(x=Mean.Macrotime..s., y=Number.of.Photons..BR.), color = 'gray') +
  labs(
    x = "Time (s)",
    y = "Counted photons")+
  geom_hline(yintercept= thresBB, linetype="dashed", color = "red")+
  geom_hline(xintercept=thresBG, linetype="dashed", color = "blue")+
  geom_hline(xintercept=thresBR, linetype="dashed", color = "gray")+
  theme_bw() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=15))
dev.off()

#####
##Duration - number of photons and marker expression plotting
###
tiff("INSERT FILE NAME HERE.tiff", units="in", width=10, height=4, res=300)
ggplot() +
  geom_density_2d(data=df, aes(x=Duration..ms., y=Number.of.Photons), size = 1) +
  labs(
    x = "Duration (ms)",
    y = "Counted photons")+
  theme_bw() +
  theme(axis.text=element_text(size=13),
        axis.title=element_text(size=14))+
  guides(colour = guide_legend(override.aes = list(size=10)))
dev.off()

####
##t-SNE dimensionality reduction - visualization of multidimensional data
####
##Load libraries (use install.packages to install if requried)
library(readr)
library(Rtsne)

##load color palette
cbp1 <- c( "#E69F00", "#56B4E9", "#009E73",
           "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

set.seed(123) # for reproducibility

##scale the data of interest
dfscaled <- scale(df[,c(37, 59,60,61)], center = TRUE, scale = TRUE) #duration 37, Photons 59,60,61, Count rate 41, 42, 43

##Apply t-SNE
tsne <- Rtsne(dfscaled, dims = 2, perplexity=50, verbose=TRUE, max_iter = 500, learning = 200, check_duplicates = FALSE)

tsne = as.data.frame(tsne$Y)

tiff("INSERT FILE NAME HERE.tiff", units="in", width=4, height=4, res=300)
ggplot(tsne, aes(x=V1, y=V2, group = df$sample)) +  
  geom_point(aes(color = factor(df$sample)), size=0.5) +
  xlab("tSNE 1") + ylab("tSNE 2") +
  theme_light(base_size=15) +
  labs(color = "Samples") + 
  guides(color = guide_legend(override.aes = list(size = 5)))+
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=15))+
  scale_colour_manual(values=cbp1)
dev.off()

###Code to automatically fine tune the perplexity parameter
tsne_plot <- function(perpl=30,iterations=500,learning=200){
  set.seed(123) # for reproducibility
  tsne <- Rtsne(dfscaled, dims = 2, perplexity=perpl, verbose=TRUE, max_iter=iterations, eta=learning, check_duplicates = FALSE)
  
  tsne = as.data.frame(tsne$Y)
  
  p <- ggplot(tsne, aes(x=V1, y=V2, color = df$name)) +
    geom_point(size=0.6) +
    xlab("") + ylab("") +
    ggtitle(paste0("perpl " , perpl , ' iterations ' , iterations , ' learning ' , learning )) +
    theme_light(base_size=20) +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank()) +
    theme(legend.position="none")
  theme_bw()
  print(p)
}

perplexity_values <- c(2,5, 30,50,100)
sapply(perplexity_values,function(i){tsne_plot(perpl=i)})



