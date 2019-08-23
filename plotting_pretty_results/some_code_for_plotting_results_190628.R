# Here is some code I've developed for plotting the results of a single Netlogo run,
# data exported from Behaviorspace.
# Here, the model that I ran was pheno-evo_190701.nlogo, using env3 as initial conditions
# and the data file output is called env3-exp_190629-table.csv

setwd("/Volumes/GoogleDrive/My Drive/From Desktop/BoL_docs/Complexity/micro-pop-model/data")
library(data.table)
library(ggplot2)
library(stringr)
library(reshape2)
library(DescTools)
library(ggTimeSeries)
library(RColorBrewer)

# read in data - this takes a few minutes
data<-fread('env3-exp_190629-table.csv', skip=6, header=T, 
            sep=',', check.names=T)
colnames(data)[19]<-'step'
colnames(data)[21:25]<-c('mean.toxin', 'mu.vals','sig.vals','barcodes','degrade.rate')

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

### process data on lineage barcodes
### this may be done either for relative abundances or absolute abundances.
### currently it's set for relative abundances.
### if you want to do absolute abundances, switch out the commented lines as indicated.
t1<-parsefactors.norm(data$barcodes[1])[,c(1,3)] # make a base dataframe (relative abundance)
# t1<-parsefactors(data$barcodes[1]) # (absolute abundance option)
for (timepoint in seq(2,dim(data)[1],100)){ # loop through every 100th timepoint 
  t1<-merge(t1, parsefactors.norm(data$barcodes[timepoint])[,c(1,3)], by='s3', all=T) # (relative abundance)
  #  t1<-merge(t1, parsefactors(data$barcodes[timepoint]), by='s3', all=T) # (absolute abundance option)
}
colnames(t1)<-c('barcode', 1, seq(2,dim(data)[1],100))
barcodes<-melt(t1, id.vars='barcode', 
                   variable.name='timepoint', value.name='abundance')
barcodes$barcode<-factor(barcodes$barcode)
barcodes$timepoint<-as.numeric(as.character(barcodes$timepoint))
barcodes$abundance[is.na(barcodes$abundance)]<-0
rm(t1)

### and plot! 
# Here, we use a stacked bar plot and turn it sideways
p1<-ggplot(barcodes, aes(y=timepoint, x=barcode))
p1<-p1 + geom_bar(aes(fill=abundance, colour=abundance), stat='identity') + coord_flip()
p1<-p1 + scale_fill_gradient2(low='#ffffcc', mid='#41b6c4', high='#08306b', name='Number of\nturtles')
p1<-p1 + scale_colour_gradient2(low='#ffffcc', mid='#41b6c4', high='#08306b', name='Number of\nturtles')
p1<-p1 + ylab('Time (ticks)')
p1<-p1 + theme(axis.text.x=element_blank()) + theme_bw()
p1

# this version uses stacked bar plots the normal way
p3<-ggplot(barcodes, aes(x=timepoint, y=abundance, fill=barcode, color=barcode))
p3<-p3 + geom_bar(stat='identity', position='fill', width=1)
p3<-p3 + xlab('Time (ticks)')
p3<-p3 + theme_bw() + theme(legend.position='none')
p3

# this is the best version: st(r)eam graph from the ggTimeSeries package
p4<-ggplot(barcodes, aes(x=timepoint, y=abundance,
                         group=barcode, fill=barcode))
p4<-p4 + stat_steamgraph(color='black', size=0.1) 
p4<-p4 + scale_fill_manual(values=rep(c(brewer.pal(12, 'Paired'), 
                                        rev(brewer.pal(12, 'Paired'))),5))
p4<-p4 + theme_bw() + theme(legend.position='none') 
p4<-p4 + xlab('Timepoint (ticks)') + ylab('Barcode relative abundance')
p4

# export plot
pdf('barcode_streamgraph.pdf', width=10, height=7)
p4
dev.off()

png('barcode_streamgraph.png', width=10, height=7, res=200, units='in')
p4
dev.off()
######################################################################################

### process data on degradation rates 
t3<-parsenumeric(data$degrade.rate[1])
t3$s4<-as.numeric(t3$s4)
pheno.ent<-data.frame(timepoint=1, entropy=Entropy(t3$s4))
for (tick in seq(2,dim(data)[1],100)){ # loop through every 1000th timepoint 
  t3<-cbind(t3, parsenumeric(data$degrade.rate[tick]))
  pheno.ent<-rbind(pheno.ent, data.frame(timepoint=tick, 
                                         entropy=Entropy(parsenumeric(data$degrade.rate[tick])$s4)))
}
colnames(t3)<-c(1, seq(2,dim(data)[1],100))
t3$degrade.rate<-seq(0.01, 1, 0.01)
degraderates<-melt(t3, id.vars='degrade.rate', 
                   variable.name='timepoint', value.name='abundance')
degraderates$timepoint<-as.numeric(as.character(degraderates$timepoint))
rm(t3)

# plot!
p2<-ggplot(degraderates, aes(x=timepoint, y=degrade.rate, fill=abundance))
p2<-p2 + geom_tile(width=100, height=0.01)
p2<-p2 + scale_fill_gradient(low='#ffffd9', high='#081d58', name='% of \npopulation')
p2<-p2 + ylab('Degradation rate') + xlab('Time (ticks)')
p2
# export plot
pdf('degraderateplot.pdf', width=7, height=7)
p2
dev.off()

png('degraderateplot.png', width=7, height=7, res=200, units='in')
p2
dev.off()

# plot entropy
p9<-ggplot(data=pheno.ent, aes(x=timepoint, y=entropy))
p9<-p9 + geom_line() + theme_bw() + ggtitle('Shannon entropy of population degrade rates')
p9<-p9 + xlab('Timepoint (ticks)') + ylab('Entropy')
p9

png('entropyplot.png', width=5, height=4, res=200, units='in')
p9
dev.off()
######################################################################################

### process data on genotypes
### this is more complicated because each cell has 5 mu values and 5 sigmas
### here, I'm only processing data from timepoint 25000.
mustring<-data$mu.vals[25000]
s1<-str_replace_all(as.character(mustring),'\\] \\[','_')
s2<-str_replace(as.character(s1),'\\[\\[','')
s3<-str_replace(as.character(s2),'\\]\\]','')
s4<-strsplit(s3, split='_')[[1]]
length(unique(s4)) ## before proceeding, check how many unique genotypes there are!
s5<-data.frame(table(s4)) # there are actually only 24 unique genotypes
rm(s1, s2, s3, s4)

### here is a function to calculate the summed distribution
### from any 5 mu and sigma values
indiv.genotype<-function(ind.muvals, ind.sigvals){
  genotype<-0
  for (d.index in 1:length(ind.muvals)){
    newdist<-dunif(x=seq(0, 1, 0.01),
                   min=ind.muvals[d.index]-0.5*ind.sigvals[d.index],
                   max=ind.muvals[d.index]+0.5*ind.sigvals[d.index])
    genotype<-genotype+(newdist/sum(newdist)) # normalize to equal 1
  }
  return(genotype) # so all the densities actually sum to 5, but it's not a big deal
}


fake.sigs<-rep(0.2, 5) # something screwed up with the export of the sigma values, so I'm replacing them for now
timepoint.genotypes<-data.frame(cellnum=numeric(), degrade.rate=numeric(), 
                                density=numeric(), abundance=numeric(), entropy=numeric())
for (cell in 1:dim(s5)[1]){
  cell.mus<-as.numeric(strsplit(as.character(s5$s4[cell]), split=' ')[[1]])
  cell.genotype<-indiv.genotype(cell.mus, fake.sigs)
  cell.abund<-s5$Freq[cell]
  cell.entropy<-Entropy(cell.genotype)
  cell.df<-data.frame(cellnum=rep(cell, 101), degrade.rate=seq(0, 1, 0.01), 
                      density=cell.genotype, abundance=cell.abund, entropy=cell.entropy)
  timepoint.genotypes<-rbind(timepoint.genotypes, cell.df)
}
timepoint.genotypes$cellnum<-factor(timepoint.genotypes$cellnum)
rm(s5, cell.df)

# plot!
p7<-ggplot(data=timepoint.genotypes, aes(x=degrade.rate, y=density, group=cellnum))
p7<-p7 + geom_col(aes(fill=abundance, color=abundance))
p7<-p7 + facet_grid(cellnum~.) + theme_bw() + ggtitle('Genotypes at timepoint 25000')
p7<-p7 + xlab('Degradation rate') + ylab('Probability density')
p7
# export plot
pdf('genotypeplot.pdf', width=7, height=12)
p7
dev.off()

png('genotypeplot.png', width=7, height=12, res=200, units='in')
p7
dev.off()


######################################################################################

### now, the easy stuff
df1<-subset(data, select=c(step, mean.toxin, count.turtles))
df2<-df1[1:1000,]
df3<-df1[24000:25000,]

p8<-ggplot(data=df3, aes(x=step))
p8<-p8 + geom_line(aes(y=mean.toxin))
p8<-p8 + geom_line(aes(y=count.turtles/26000), color='red')
p8<-p8 + scale_y_continuous(sec.axis = sec_axis(~.*26000, name='Number of Turtles'))
p8<-p8 + theme_bw() + xlab('Timepoint (ticks)') + ylab('Average toxin per patch')
p8

# export plot
pdf('populationplot_start.pdf', width=8, height=3)
p8
dev.off()



######################################################################################
######################################################################################
######################################################################################

# playing around with different ways to plot distributions

# a fake dataframe
df<-data.frame(mu=c(0.8, 0.1, 0.2, 0.4, 0.35), sig=rep(0.2, 5), height=c(0.1, 0.4, 0.3, 0.5, 0.5))

# use runif to generate random draws from the distribution and plot a histogram
hist(c(runif(n=10000, min=df$mu[1]-0.5*df$sig[1], max=df$mu[1]+0.5*df$sig[1]),
       runif(n=10000, min=df$mu[2]-0.5*df$sig[2], max=df$mu[2]+0.5*df$sig[2])),
     xlim=c(0, 1))

# use dunif to get the density distributions, sum, and plot
dist.d<-0
for (d.index in 1:length(df$mu)){
  dist.d<-dist.d+
    dunif(x=seq(0, 1, 0.01),
          min=df$mu[d.index]-0.5*df$sig[d.index],
          max=df$mu[d.index]+0.5*df$sig[d.index])
}
plot(seq(0, 1, 0.01), dist.d, type='l', col='red')




