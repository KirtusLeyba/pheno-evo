## Some code for plotting the output of the BehaviorSpace run
## of the pheno-evo_19030.nlogo script
## on August 30, 2019 - testing combinations of 
## initial.response.error, initial.switch.rate, env.noise, diff.rate
## (no evolution). There were 16 runs, no replicate runs.
## The code below separates each of the runs into its own dataframe
## and then loops through all dataframes to extract and plot each of several outputs.
## This is very similar to some_code_for_plotting_results_190628.R
## except that many things are now written as functions so as to make looping easy.

setwd("/Volumes/GoogleDrive/My Drive/From Desktop/BoL_docs/Complexity/micro-pop-model/results")
library(data.table)
library(ggplot2)
library(stringr)
library(reshape2)
library(DescTools)
library(ggTimeSeries)
library(RColorBrewer)
library(dplyr)

# read in data - this takes a few minutes
data<-fread('pheno-evo_exp190830.csv', skip=6, header=T, 
            sep=',', check.names=T)
colnames(data)[1]<-'run.number'
colnames(data)[11]<-'step'
colnames(data)[c(13,14,15)]<-c('mean.toxin','patch.toxin','degrade.rate')

summary_data<-distinct(data[,c(1,3,4,6,9)])
write.csv(summary_data, file='run_characteristics.csv', row.names=F)

######################################################################################

### here are some functions that parse the lists of values
### spit out by NetLogo at each timepoint

# parse a string of numeric values (such as degrade rates)
parsenumeric<-function(charstring){
  s1<-str_replace(as.character(charstring),'\\[','')
  s2<-str_replace(as.character(s1),'\\]','')
  s3<-as.numeric(strsplit(s2, split=' ')[[1]])
  s4<-hist(s3, breaks=seq(0, 1, 0.01), plot=F)$density
  return(data.frame(s4))
}

# parse a string of character values that will be treated as factors
# (e.g., barcode names, which look like numbers)
parsefactors<-function(charstring){
  s1<-str_replace(as.character(charstring),'\\[','')
  s2<-str_replace(as.character(s1),'\\]','')
  s3<-as.numeric(strsplit(s2, split=' ')[[1]]) # a list of all the barcodes
  s4<-data.frame(table(s3)) # a table of frequencies of all the barcodes (s3 = barcode name)
  s4$s3<-as.character(s4$s3) # make the barcode names into character strings 
  return(data.frame(s4))
}

# this one is like parsefactors but normalizes the abundance of each barcode name at each timepoint
# and calculates the proportion of the total population made up by that barcode name
parsefactors.norm<-function(charstring){
  s1<-str_replace(as.character(charstring),'\\[','')
  s2<-str_replace(as.character(s1),'\\]','')
  s3<-as.numeric(strsplit(s2, split=' ')[[1]]) # a list of all the barcodes
  s4<-data.frame(table(s3)) # a table of frequencies of all the barcodes (s3 = barcode name)
  s4$s3<-as.character(s4$s3) # make the barcode names into character strings 
  s4$prop<-s4$Freq/sum(s4$Freq)
  return(data.frame(s4))
}

######################################################################################

### divide data into individual runs

# no need to keep all the columns of data:
keep.colnames<-c('initial.response.error','env.noise', 'initial.switch.rate','diff.rate',
                 'step','count.turtles','mean.toxin','patch.toxin','degrade.rate')

# loop through each run number and separate it out,
# trim each dataframe to just the entries where the population isn't 0
for (runnum in unique(data$run.number)){
  newdf<-subset(data, run.number==runnum, select=keep.colnames)
  if (length(which(newdf$count.turtles==0))>0){
    lastdatapoint<-min(which(newdf$count.turtles==0))
    newdf<-newdf[1:lastdatapoint,]
    }
  assign(paste0('data',runnum), newdf) 
  }

# delete the original dataframe
rm(data)

######################################################################################

### process data on degradation rates 

# Degradation rates of each turtle through time
process.degrade.rates<-function(df){
  t3<-parsenumeric(df$degrade.rate[1])
  t3$s4<-as.numeric(t3$s4)
  if (dim(df)[1]>1000){
    for (tick in seq(100,dim(df)[1],100)){ # loop through every 100th timepoint 
      t3<-cbind(t3, parsenumeric(df$degrade.rate[tick]))}
    colnames(t3)<-c(1, seq(100,dim(df)[1],100))
    w<-100
  } 
  else {
    for (tick in seq(2,dim(df)[1],1)){ # take every timepoint 
      t3<-cbind(t3, parsenumeric(df$degrade.rate[tick]))}
    colnames(t3)<-c(1, seq(2,dim(df)[1],1))
    w<-1
  }
  t3$degrade.rate<-seq(0.01, 1, 0.01)
  degraderates<-melt(t3, id.vars='degrade.rate', 
                     variable.name='timepoint', value.name='abundance')
  degraderates$timepoint<-as.numeric(as.character(degraderates$timepoint))
  p2<-ggplot(degraderates, aes(x=timepoint, y=degrade.rate, fill=abundance))
  p2<-p2 + geom_tile(width=w, height=0.01)
  p2<-p2 + scale_fill_gradientn(colours=brewer.pal(11, 'Spectral'), name='% of \npopulation')
  #p2<-p2 + scale_fill_gradient(low='#ffffd9', high='#081d58', name='% of \npopulation')
  p2<-p2 + ylab('Degradation rate') + xlab('Time (ticks)')
  p2
  p2<-p2 + ggtitle(substitute(df))
  return(p2)
  }

# Entropy of the population's toxin degradation rates, a summary statistic
process.entropy<-function(df){
  t3<-parsenumeric(df$degrade.rate[1])
  t3$s4<-as.numeric(t3$s4)
  pheno.ent<-data.frame(timepoint=1, entropy=Entropy(t3$s4))
  if (dim(df)[1]>1000){
    for (tick in seq(100,dim(df)[1],100)){ # loop through every 100th timepoint 
      pheno.ent<-rbind(pheno.ent, data.frame(timepoint=tick, 
                                             entropy=Entropy(parsenumeric(df$degrade.rate[tick])$s4)))}
  } 
  else {
    for (tick in seq(2,dim(df)[1],1)){ # take every timepoint 
      pheno.ent<-rbind(pheno.ent, data.frame(timepoint=tick, 
                                             entropy=Entropy(parsenumeric(df$degrade.rate[tick])$s4)))}
  }
  p9<-ggplot(data=pheno.ent, aes(x=timepoint, y=entropy))
  p9<-p9 + geom_line() + theme_bw() + ggtitle('Shannon entropy of population degrade rates')
  p9<-p9 + xlab('Timepoint (ticks)') + ylab('Entropy')
  p9<-p9 + ggtitle(substitute(df))
  return(p9)
}


######################################################################################

## Extract the value of toxin concentration at each patch over time

# This function is similar to parsenumeric, above,
# but collapses the toxin concentration values so that anything above 1 is recoded as 1
# to make plotting easier.
parsetoxin<-function(charstring){
  s1<-str_replace(as.character(charstring),'\\[','')
  s2<-str_replace(as.character(s1),'\\]','')
  s3<-as.numeric(strsplit(s2, split=' ')[[1]])
  s3[which(s3>1)]<-1 # if toxin is above 1, recode it as 1 
  s4<-hist(s3, breaks=seq(0, 1, 0.01), plot=F)$density
  return(data.frame(s4))
}

# This function uses parsetoxin to carry out all toxin processing for the dataframe.
process.patch.toxin<-function(df){
  t3<-parsetoxin(df$patch.toxin[1])
  t3$s4<-as.numeric(t3$s4)
  if (dim(df)[1]>1000){
    for (tick in seq(100,dim(df)[1],100)){ # loop through every 100th timepoint 
      t3<-cbind(t3, parsetoxin(df$patch.toxin[tick]))}
    colnames(t3)<-c(1, seq(100,dim(df)[1],100))
  } 
  else {
    for (tick in seq(2,dim(df)[1],1)){ # take every timepoint 
      t3<-cbind(t3, parsetoxin(df$patch.toxin[tick]))}
    colnames(t3)<-c(1, seq(2,dim(df)[1],1))
  }
  t3$patch.toxin<-seq(0.01, 1, 0.01)
  toxinlevels<-melt(t3, id.vars='patch.toxin', 
                     variable.name='timepoint', value.name='abundance')
  toxinlevels$timepoint<-as.numeric(as.character(toxinlevels$timepoint))
  if (dim(t3)[2]<100){w<-1} else {w<-100}
  p2<-ggplot(toxinlevels, aes(x=timepoint, y=patch.toxin, fill=abundance))
  p2<-p2 + geom_tile(width=w, height=0.01)
  p2<-p2 + scale_fill_gradient(low='white', high='black', name='% of \npatches')
  p2<-p2 + ylab('Patch toxin level') + xlab('Time (ticks)')
  p2<-p2 + ggtitle(substitute(df))
  return(p2)
}

######################################################################################

### now, the easy stuff

process.count.turtles<-function(df){
  df1<-subset(df, select=c(step, mean.toxin, count.turtles))
  if (dim(df1)[1]>2000){
    df1<-df1[dim(df1)[1]-500:dim(df1)[1],]}
  scaling<-max(df1$count.turtles)/max(df1$mean.toxin)
  p8<-ggplot(data=df1, aes(x=step))
  p8<-p8 + geom_line(aes(y=mean.toxin))
  p8<-p8 + geom_line(aes(y=count.turtles/scaling), color='red')
  p8<-p8 + scale_y_continuous(sec.axis = sec_axis(~.*scaling, name='Number of Turtles'))
  p8<-p8 + theme_bw() + xlab('Timepoint (ticks)') + ylab('Average toxin per patch')
  p8<-p8 + ggtitle(substitute(df))
  p8
}


######################################################################################

### export plots
setwd("/Volumes/GoogleDrive/My Drive/From Desktop/BoL_docs/Complexity/micro-pop-model/results/190902")

# I would have liked to set up some code to loop through each of the dataframes...
# but for now it's easier just to do find-and-replace to run the code below
# for each of the 16 datasets.

jpeg('data16_degraderateplot.jpg', width=7, height=7, res=200, units='in')
process.degrade.rates(data16)
dev.off()

jpeg('data16_entropyplot.jpg', width=5, height=4, res=200, units='in')
process.entropy(data16)
dev.off()

jpeg('data6_toxinplot.jpg', width=7, height=7, res=200, units='in')
process.patch.toxin(data6)
dev.off()

jpeg('data4_populationplot.jpg', width=8, height=3, res=200, units='in')
process.count.turtles(data4)
dev.off()
