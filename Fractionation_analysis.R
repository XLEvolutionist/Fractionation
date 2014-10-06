rm(list = ls())
require(data.table)
library(Matrix)
library(ggplot2)
library(edgeR)
library(plyr)
library(fields)
library(reshape)
library(ggthemes)
library(grid)
library(grDevices)
library(mgcv)
library(data.table)
library(GenomicFeatures)


######
#######################
##FUNCTIONS
#######################
######


#########
#this function calculates the number of Tc genes that span the block
#########
total_blk_size<-function(x) {
  limits<-match(c("B","E"), x)
  return(length(limits[1]:limits[2])) 
}

#########
#this function calculates the number of Gr genes that span the block
#return the number of genes retaiend and lost
#########
#x is a vector, usually a column from the synteny table (i.e. a block of genes in synteny)
total_gr_size<-function(x) {
  limits<-match(c("B","E"), x)
  x<-x[limits[1]:limits[2]]
  return(cbind(sum(table(x)[2:4]),table(x)[1]))
}

########
#A function to calculate % retention in Gr
########
#x = the number of genes retained in Gr
#y= the total nuber of genes in the block
per_cent_calc<-function(x,y) {
  return((x/y)*100)
}

####
# A function to find out what chromosome each block is on
####

chrDesignation<-function(x) {
  a<-unique(x)
  print(a[3])
}#chrDesignation

###
#A function to return the start and stop co-ordinates of each block
###

BlockCoords<-function(x) {
  c<-match(c("B","E"), x)
  if (length(c) != 2 ) {
    print ( paste("Error block length = ", length(c)))
  }
  return(match(c("B","E"), x))
}

#a slightly modifed version to use for the when Gorai names are used instead "B", "E", "-"

BlockCoordsGr<-function(x) {
  c<-which(x != "-")
  if (length(c) < 2 ) {
    print ( paste("Error block length = ", length(c)))
  }
  return(c(x[c[1]], x[rev(c)[1]]))
}


######
#A function to find if there is overlap between to blocks
######
#x and y are columns of the synteny table
#len is the minimum length of overlap for the comparison
overlap<-function(x,y,len) {
  limitsX<-match(c("B","E"), x)
  #limitsX<-na.omit(limitsX)
  limitsY<-match(c("B","E"), y)
  #limitsY<-na.omit(limitsY)
  #do the ranges overlap
  length_to_test<-length(intersect(limitsX[1]:limitsX[2],limitsY[1]:limitsY[2]))
  if ( length_to_test >= len ) {
    m<-cbind(limitsX,limitsY)
    start<-max(m[1,])
    end<-min(m[2,])
    x<-x[start:end]
    y<-y[start:end]
    xy<-list("x"=x,"y"=y)
    return(xy)
  }#
  else {
    return(0)
  }#else  
}#overlap

#######
#A function to perform a fisher's exact test
#######
#x, and y are vectors of length 2 that contain thenumber of lost 
#and retained genes for each of the blocks that are being compared
f.test<-function(x) {
  dm<-x
  test<-fisher.test(dm, alternative="two.sided")
  return(cbind(test$statistic,test$p.value))
}#f.test
####
#The same with a Chi Sq test
####
c.test<-function(x) {
  dm<-x
  test<-chisq.test(dm)
  return(cbind(test$statistic,test$p.value))
}#f.test

######
#An error test
######

error.test<-function(x) {
  if ( length(unique(x)) == 1 ) {
    print("Error")
    print(unique(x))
  }#if
}#

###
#Compare blocks, figure out shared pairs and return the number of genes upregulated in each block 
###
gene.overlap<-function(x,y,len,exp) {
  #find the position of retained genes
  x.index<-which(x != "-")
  y.index<-which(y != "-")
  #now find shared retained homeologs?
  inter<-intersect(x.index,y.index)
  if ( length(inter) == 0 ) {
    b<-list("x"=0,"y"=0,"p.value"="NA")
    return(b)
    stop("Zero length, none of the genes overlap")
  }#if
  #get the names of the genes to grab
  x.names<-x[inter]
  y.names<-y[inter]
  
  #grab the gene expression per block
  x.exp<-expAve[x.names,]
  y.exp<-expAve[y.names,]
  #calculate differences
  diff.exp<-x.exp-y.exp
  y.large<-length(diff.exp[diff.exp < 0])
  x.large<-length(diff.exp[diff.exp > 0])
  
  if ( is.null(dim(diff.exp)) ) { 
    print(x.large)
    a<-list("x"=x.large,"y"=y.large,"p.value"=binom.test(x.large,(length(diff.exp)), 0.5)$p.value)
    return(a)
  }#if
  else {
    a<-list("x"=x.large,"y"=y.large,"p.value"=binom.test(x.large,(dim(diff.exp)[1]*dim(diff.exp)[2]), 0.5)$p.value)
    return(a)
    
  }#else
}#gene.overlap

#####
#a function that sorts blocks according to level of fractionation
#####
#x and y and block names (numbers)
sortBlock<-function(x,y) {
  #find the % retention
  x_per<-fracMat[x,4]
  y_per<-fracMat[y,4]
  if ( x_per == y_per ) {
    print("Equally fractionated")
  }#if
  else {
    if ( x_per > y_per ) {
      return(list(c(x,y)))
    }
    else {
      return(list(c(y,x)))
    }
  }
}#sortBlock

####
#Return blocks in order of MF to LF
####

mostFrac<-function(x,z) {
  fracData<-z
  pairs<-x
  #print(b)
  #if ( b )
  output<-fracData[as.character(x),4]
  print(output)
  min2max<-pairs[order(output)]
  #print(order(output))
  output<-c(output)
  return(order(output))
}#mostFrac

###
#return gene expr data by  rownames
###

returnExpr<-function(names,exp) {
  return(exp[names,])
}

###
#bin genes into LF and MF
###

binFrac<-function(LF, ListExp) {
  print(LF[,1])
  MFexp<-ListExp[LF[,1],]
}

####
# A function that calculates mean expression per tissue
####

expPerTissue<-function(x) {
  mean1<-mean(x[1:3])
  mean2<-mean(x[4:6])
  mean3<-mean(x[7:9])
  return(c(mean1,mean2,mean3))
}

#a function the compresses a matrix x in to a single column replacing NAs in a 
#a column with the value in the the previous column, works from the last column to 
#the first. The effect is to combines columns, replacing NAs with data from other
#colums
#x<-cbind(c(NA,"BOB",NA),c("BOB",NA,NA))
#> compress(x)
#[1] 2
#[1] "BOB" "BOB" NA   

compress<-function(x, test) {
  test<-x
  #assign('test',test,envir=.GlobalEnv)
  for ( i in (dim(test)[2]:2) ) {
    print(paste("i is:",i))
    print(paste("i-1 is:",i-1))
    print(test)
    inx<-which(is.na(test[,i-1]))
    test[inx,i-1]<-test[inx,i]
    print(test)
  }#for
  assign('test',test,envir=.GlobalEnv)
  return(test)
  #return(y)
}#compress

####
################################################################################
# Data preparation phase
################################################################################
####

###################
#Load in and re-arrange the data 
###################

#load in the syn data with PAC IDs
PacData<-read.table("/Users/simonrenny-byfield/cotton/Fractionation/Gr_vs_Tc_BlastN/gorai_aligned_sections.txt", 
                    header = F , sep = "\t")
#load in the syntenic data
mapping<-read.table("/Users/simonrenny-byfield/cotton/Fractionation/Gr_vs_Tc_BlastN/Tc_aligned_sections.txt", header = F, sep = "\t")

#load in the Kaks data
Kaks<-read.table("/Users/simonrenny-byfield/cotton/Fractionation/Gr_vs_Tc_BlastN/KaKs_data_Gr.txt", header = F , sep = "\t")
colnames(Kaks)<-c("Ks","Ka","Gr","Tc")
Kaks<-Kaks[!duplicated(Kaks[,3]),]
rownames(Kaks)<-Kaks[,3]

#load in the 24nt sRNA data
sRNA_data<-read.table("/Users/simonrenny-byfield/cotton/Fractionation/Gr_vs_Tc_BlastN/Final_coverage_100bp_10bp_5kb_window_in_D5_genome_on_TEs.txt",
                        sep="\t", header = T)
#load in the 24nt WITHIN genes data,
withinTE<-read.table("/Users/simonrenny-byfield/cotton/Fractionation/Final_coverage_24nt_siRNAs_on_TEs_inserted_in_gene_no_strandness_TE_unmasked.txt"
                     , sep="\t", header = T)
#for genes with two or more TEs sum the 24nt RNAs
withinTE.dt<-data.table(withinTE)
withinTE.dt<-withinTE.dt[,list(RNAs=sum(RNAs)), by='gene']

#load in the TE positional data
TEpos<-read.table("/Users/simonrenny-byfield/cotton/Fractionation/Gr_vs_Tc_BlastN/NearestTE/overlap_allowed/gene2nearestTE.txt",
                    sep = "\t", header = F)
colnames(TEpos)<-c("gene","distance")

#load in the gene positional info for ideograms
posGr<-read.table("/Users/simonrenny-byfield/cotton/Fractionation/Gr_vs_Tc_BlastN/Gr_by_chromsome.txt", header = F , sep = "\t")
posGr<-unique(posGr)
colnames(posGr)<-c("pac.id","chr","range.bp")
rownames(posGr)<-posGr[,1]

#load in TE coverage
TEcoverage<-read.table("/Users/simonrenny-byfield/cotton/Fractionation/Gr_vs_Tc_BlastN/Final_TE_proportions_in_sliding_window_5kb_up_downstream_D5_genes.txt", 
                          header = F, sep ="\t")

#add a - to all windows upstream of the genes
TEcoverage[,4][TEcoverage[,3] == "up"]<- -502 + abs(TEcoverage[,4][TEcoverage[,3] == "up"])

###
#Load in the GC data

#this is only a subset of genes, redo to contain the complete set
gcTable<-read.table("/Users/simonrenny-byfield/cotton/Fractionation/Gr_vs_Tc_BlastN/GC/GC_tableR.txt", 
                    sep = "\t" , header = FALSE, fill = TRUE)
#replace NA with zero
gcTable[is.na(gcTable)]<-0
rownames(gcTable)<-gcTable[,1]
gcTable<-gcTable[,-1]

####
#Clean and re-arrange the data
####

#expression data in RPKM
expression_RPKM<-load(file ="/Users/simonrenny-byfield/cotton/Fractionation/Gr_vs_Tc_BlastN/expressionRPKMordered.RData")
test<-load(file ="/Users/simonrenny-byfield/cotton/Fractionation/Gr_vs_Tc_BlastN/expressionRPKMordered.RData")

####
#Calculate mean expression per tissue

petalAve<-rowMeans(expression_RPKM[,1:3])
leafAve<-rowMeans(expression_RPKM[,4:6])
seedAve<-rowMeans(expression_RPKM[,7:9])

#combine the mean expr data
expAve<-cbind(petalAve,leafAve,seedAve)
colnames(expAve)<-c("petal","leaf","seed")
rownames(expAve)<-rownames(expression_RPKM)

###
#Clean the syntenic block info

#name the rows according the Tc gene name
rownames(mapping)<-mapping[,2]
rownames(PacData)<-mapping[,2]

#remove the column of gene names (now in rownames(mapping))
mapping<-mapping[,-2]
#PacData<-PacData[,-length(PacData[1,])]
PacData<-PacData[,-1]

#declare the cotton chromosomes, used for later to identify gene matches
c_chr<-c(1:13)

#remove all those columns which contain nothing but "-"
mapping2<-mapping[, apply(mapping,2,function(x) { length(unique(x))>1 }) ] 
PacData2<-PacData[, apply(mapping,2,function(x) { length(unique(x))>1 }) ]

#remove all rows with only "-". While these have genes in Tc they have no corresponding 
#gene mapped to Gr in any of the aligned syntenic blocks and so may well be insertions in Tc
mapping3<-mapping2[apply(mapping2,1,function(x) { length(unique(x))>1 }), ] 
PacData3<-PacData2[apply(mapping2,1,function(x) { length(unique(x))>1 }), ]
tc_chr<-mapping3[,1]
mapping3<-mapping3[,-1]
PacData3<-PacData3[,-202]
colnames(mapping3)<-1:length(mapping3[1,])
colnames(PacData3)<-colnames(mapping3)

#find all the unique combinations of the blocks. We will see from all of these which can
#be tested. i.e. those have overlaps.
combination<-combn(c(1:dim(mapping3)[2]),m=2)

compListCoords<-NULL

#output a table of the coverage per Tc gene

Tc_depth<-apply(mapping3,1,function(x) length(x[x != "-"]))
Tc_depth<-cbind("chr"=tc_chr,"depth"=Tc_depth,"gene_index"=1:length(tc_chr))

Tc_depth[,3][Tc_depth[,1]==2]<-Tc_depth[,3][Tc_depth[,1]==2]-cumsum(table(tc_chr))[1]
Tc_depth[,3][Tc_depth[,1]==3]<-Tc_depth[,3][Tc_depth[,1]==3]-cumsum(table(tc_chr))[2]
Tc_depth[,3][Tc_depth[,1]==4]<-Tc_depth[,3][Tc_depth[,1]==4]-cumsum(table(tc_chr))[3]
Tc_depth[,3][Tc_depth[,1]==5]<-Tc_depth[,3][Tc_depth[,1]==5]-cumsum(table(tc_chr))[4]
Tc_depth[,3][Tc_depth[,1]==6]<-Tc_depth[,3][Tc_depth[,1]==6]-cumsum(table(tc_chr))[5]
Tc_depth[,3][Tc_depth[,1]==7]<-Tc_depth[,3][Tc_depth[,1]==7]-cumsum(table(tc_chr))[6]
Tc_depth[,3][Tc_depth[,1]==8]<-Tc_depth[,3][Tc_depth[,1]==8]-cumsum(table(tc_chr))[7]
Tc_depth[,3][Tc_depth[,1]==9]<-Tc_depth[,3][Tc_depth[,1]==9]-cumsum(table(tc_chr))[8]
Tc_depth[,3][Tc_depth[,1]==10]<-Tc_depth[,3][Tc_depth[,1]==10]-cumsum(table(tc_chr))[9]

#write out the table for use in circos
write.table(Tc_depth, file="/Users/simonrenny-byfield/cotton/Fractionation/Gr_vs_Tc_BlastN/circos_depth_data.txt", quote=FALSE, sep = "\t")
save(Tc_depth, file="/Users/simonrenny-byfield/cotton/Fractionation/Gr_vs_Tc_BlastN/circos_depth_data.RData")



########
###############################################################################
#Analysis phase
###############################################################################
########

#foreach block record the number of genes in Gr, the total in Tc and the % retained in Gr

#use apply to calculate ancestral block size for each column ( i.e. each "block")
Tc_genes_blk<-apply(mapping3, 2 , function(x) total_blk_size(x))

#use apply to calculate number of retaiend and lost Gr genes for each column ( i.e. each "block")
Gr_genes_blk<-apply(mapping3, 2 , function(x) total_gr_size(x))

#calculate the % gene retention of Gr genes
Gr_per_cent<-per_cent_calc(Gr_genes_blk[1,],Tc_genes_blk)

#find the beginning and end of each block
block_start_end<-apply(mapping3, 2, function(x) BlockCoords(x))
SB<-matrix(block_start_end,ncol = 2 , byrow = T)

#do the same for the PAC numbers, tranlsate to cotton coordinates
block_start_end<-apply(PacData3, 2, function(x) BlockCoordsGr(x))
PB<-matrix(block_start_end,ncol = 2 , byrow = T)
Gr.start<-posGr[as.character(PB),3]
Gr.end<-match(as.character(PB), rownames(posGr))
Gr.end<-matrix(Gr.end, byrow =F , ncol = 2)

#Gr_start_stop<-matrix(posGr[which(as.vector(PB)),3],ncol = 2, byrow=F)

#grab the Tc chr per block
matches<-apply(mapping3,2,function(x) which(x =="B"))
tc_names<-tc_chr[as.numeric(matches)]
which(tc_chr$Gene == names)
#apply(mapping3,2,function(x) chrDesignation(x))
#make a data.frame of block size, retaiend gene number of % retained
fracMat<-matrix(c( Tc_genes_blk , Gr_genes_blk[1,] , Gr_genes_blk[2,], Gr_per_cent,
                      tc_names,as.numeric(SB),as.numeric(apply(mapping3,2,function(x) as.character(chrDesignation(x)))),
                         Gr.end),
                          nrow = length(Tc_genes_blk), ncol=10, byrow = F)
#add some column names
colnames(fracMat)<-c("Tc","Gr.retained", "Gr.lost" , "percent" , "Tc.Chr", "Tc.start" , "Tc.end" , "Gr.Chr", "Gr.start", "Gr.end")
rownames(fracMat)<-1:length(fracMat[,1])
#order the matrix based on Tc, Gr and then start position
fracMatOrd<-fracMat[order(fracMat[,5], fracMat[,8], fracMat[,6]),]

#table has been rearranged by Tc chr, Gr chrthen by Tc star , then Gr start (saved in Fractionation_comparision_v3.R)

#write out the FracMatCircos matrix
save(fracMatCircos, file = "/Users/simonrenny-byfield/cotton/Fractionation/Gr_vs_Tc_BlastN/FracMatCircos.RData")
#load(file = "/Users/simonrenny-byfield/cotton/Fractionation/Gr_vs_Tc_BlastN/FracMatCircos.RData")
write.table(fracMatCircos, file = "/Users/simonrenny-byfield/cotton/Fractionation/Gr_vs_Tc_BlastN/circos_link_data.txt", 
            quote=FALSE, sep = "\t")

#plot a histogram of the percent retained
#for sections over 150 genes long
ret<-fracMat[,"percent"][fracMat[,"Tc"]>0]
percentHist<-hist( log10(ret) , breaks = 15, col = "grey")

#perform a dip tst for unimodality
dip.test(ret, simulate.p.value = FALSE, B = 100000)

####
#Find all the unique block combinations using combn()
combination<-combn(c(1:dim(mapping3)[2]),m=2)

####
#Perform a Fisher's exact test on each comparable pair of blocks
####

mx<-NULL
for ( i in seq(1:length(combination[1,]))) {
  #assign each "block to a variable
  blk1<-mapping3[,combination[,i][1]]
  blk2<-mapping3[,combination[,i][2]]
  #test if the blocks overlap
  ovr<-overlap(blk1,blk2,100)
  #if the test failed go to the next comparison
  if ( !is.list(ovr) ) {
    #print("fail")
    next
  }#if
  else {
    #do a test to see if each section is not just "-", it should not be, something is badly wrong if error messages result
    error.test(ovr$x)
    error.test(ovr$y)
    #make a table of each block
    blk1totals<-table(ovr$x)
    blk2totals<-table(ovr$y)
    dm<-matrix(c(sum(blk1totals[2:4]),sum(blk1totals),sum(blk2totals[2:4]),sum(blk2totals)),
               nrow = 2,
               dimnames = list(gene = c("retained","expected"),
                               block = c(combination[,i][1],combination[,i][2])))
    #load in some info to mx
    mx<-rbind(mx, cbind( combination[,i][1] , combination[,i][2], sum(blk1totals[2:4]) , 
                            sum(blk2totals[2:4]) , sum(blk2totals[1:4]) , c.test(dm) ) ) 
  }#else
}#for


#rename the columns of matrix mx
colnames(mx)<-c("block#1","block#2","retained in 1","retained in 2","total","X-sqaured","p.value")
rownames(mx)<-NULL
#adjust the p.values using "BH" method
mx[,"p.value"]<-p.adjust(mx[,"p.value"], method = "BH")
mx.df<-as.data.frame(mx)
#print out the data in mx
#write.table(file="/Users/simonrenny-byfield/cotton/Fractionation/Gr_vs_Tc_BlastN/chisq.table.txt", mx , quote=FALSE, sep = "\t")


#find the coverage for each block
coverage<-table(c(mx[,1],mx[,2]))
#if a block matches another block, table() will gove the number 1, but the blcok is actaully in copy number 2 
#(i.e. itself AND the matching block), simialr to 3 ,4 and 5 and so on.
cov<-NULL
for ( i in seq(1:length(fracMatCircos[,1]))) {
    print(i)
    a<-coverage[as.character(i)]
    print(a)
    if ( !is.na(a)) {
      cov<-c(cov,a)
    }
    else {
      cov<-c(cov,0)
    }
}#for
cov<-as.vector(cov)+1
fracMat<-cbind(fracMat[,-c(11:12)],"coverage"=cov)

####
#Make a list of all the blocks by bin 
####

#use pretty to figure some nice breaks
breaks<-pretty(range(ret), min.n = 4)

#split the blocks into groups according to these breaks
bins<-.bincode(fracMat[,"percent"], breaks=breaks)
groups<-split(rownames(fracMat.df), bins)

#add the bin identity to the fracMat matrix
fracMat<-cbind(fracMat,"bin" = factor(bins))
fracMat.df<-as.data.frame(fracMat)

#now get the Gorai names for each group
#rapply(groups, function(x) print(unique(PacData3[,x[1:length(x)]])))

Tc.by.group<-apply(mapping3,2, function(x) rownames(mapping3)[which(x != "-")])
Gr.by.group<-apply(PacData3,2, function(x) x[which(x != "-")])

Kaks.bin<-rapply(Gr.by.group, function(x) print(Kaks[x,1]), how = "replace")
#extract expression by bin
exp.bin<-rapply(Gr.by.group, function(x) expression[x,], how = "replace")
exp.bin.melt<-melt(exp.bin)

kaks.df<-melt(Kaks.bin)
#add  the bin group to each Ks value
ks.bin<-bins[as.numeric(kaks.df[,"L1"])]
kaks.df<-cbind(kaks.df,"ks.bin"=ks.bin)


#####
#Manually assemble blocks into "whole chromosomes", where possible, according to the dotplot
#"/Users/simonrenny-byfield/cotton/Fractionation/Gr_vs_Tc_BlastN/figures/dotplot.tiff"
#####

##Thid is where most of the grunt work is done, inclusing for the bias fractionation work in the latest
##cotton fractionation paper.

#a list of the total N# of genes per Tc chromosome
tc_genes_per_chr<-table(tc_chr)

###
# For Tc chr 2
###
#for Tc chr 2/ Gr chr 5, blks 137,138,139
blk25s<-c("137","138","139")
chr25blks<-fracMatOrd[blk25s,]
chr_2_5_sums<-colSums(chr25blks[,1:3])
#for Tc chr 2/ Gr chr 8, blks 179,184,185
blk28s<-c("179","184","185")
chr28blks<-fracMatOrd[blk28s,]
chr_2_8_sums<-colSums(chr28blks[,1:3])

#[,1] [,2]
#Gr.retained  642  929
#Tc          3199 3612

#test for difference between the two

chisq.test(cbind(chr_2_8_sums[2:1],chr_2_5_sums[2:1]))

#Pearson's Chi-squared test with Yates' continuity correction

#data:  cbind(chr_2_8_sums[2:1], chr_2_5_sums[2:1])
#X-squared = 18.904, df = 1, p-value = 1.375e-05

###
# For Tc chr 6
###

#for Tc chr 6/ Gr chr 6
blk66s<-c("149","150")
chr66blks<-fracMatOrd[blk66s,]
chr_6_6_sums<-colSums(chr66blks[,1:3])

#for Tc chr 6/ Gr chr 9
blk69s<-"190"
chr_6_9_sums<-fracMatOrd[blk69s,1:3]

#            [,1] [,2]
#Gr.retained  580  147
#Tc          2353 1584

#compare 6 and 9
chisq.test(cbind(chr_6_9_sums[2:1],chr_6_6_sums[2:1]))

#  Pearson's Chi-squared test with Yates' continuity correction

#data:  cbind(chr_6_9_sums[2:1], chr_6_6_sums[2:1])
#X-squared = 104.4615, df = 1, p-value < 2.2e-16



#for Tc chr6/ Gr chr 10
blk610<-c("33","34","36")
chr610blks<-fracMatOrd[blk610,]
chr_6_10_sums<-colSums(chr610blks[,1:3])

#compare 10 and 6

#[,1] [,2]
#Gr.retained  227  147
#Tc          1973 1584

chisq.test(cbind(chr_6_10_sums[2:1],chr_6_6_sums[2:1]))

#Pearson's Chi-squared test with Yates' continuity correction


#data:  cbind(chr_6_10_sums[2:1], chr_6_6_sums[2:1])
#X-squared = 3.5429, df = 1, p-value = 0.0598


#for 10 and 9

#[,1] [,2]
#Gr.retained  227  580
#Tc          1973 2353

chisq.test(cbind(chr_6_10_sums[2:1],chr_6_9_sums[2:1]))
#Pearson's Chi-squared test with Yates' continuity correction

#data:  cbind(chr_6_10_sums[2:1], chr_6_9_sums[2:1])
#X-squared = 84.1342, df = 1, p-value < 2.2e-16

####
# For Tc chr 7
####

#for Tc chr 7/ Gr chr 2, blks 88,86,89
blk72s<-c("88","86","89")
chr72blks<-fracMatOrd[blk72s,]
chr_7_2_sums<-colSums(chr72blks[,1:3])

#for Tc chr 7/ Gr chr 13 blks 76, 75
blk713s<-c("76","75")
chr713blks<-fracMatOrd[blk713s,]
chr_7_13_sums<-colSums(chr713blks[,1:3])

#[,1] [,2]
#Gr.retained  420  225
#Tc          1843 1621

chisq.test(cbind(chr_7_2_sums[2:1],chr_7_13_sums[2:1]))

#Pearson's Chi-squared test with Yates' continuity correction
#
#data:  cbind(chr_7_2_sums[2:1], chr_7_13_sums[2:1])
#X-squared = 30.7035, df = 1, p-value = 3.006e-08

#for Tc chr 8/ Gr chr 5, blks 133,132
blk85s<-c("133","132")
chr85blks<-fracMatOrd[blk85s,]
chr_8_5_sums<-colSums(chr85blks[,1:3])

#for Tc chr 8/ Gr chr 9, blks 191
blk89s<-c("191")
chr89blks<-fracMatOrd[blk89s,]
chr_8_9_sums<-chr89blks[1:3]

#[,1] [,2]
#Gr.retained  608  236
#Tc          1784 1799

chisq.test(cbind(chr_8_9_sums[2:1],chr_8_5_sums[2:1]))

#Pearson's Chi-squared test with Yates' continuity correction

#data:  cbind(chr_8_9_sums[2:1], chr_8_5_sums[2:1])
#X-squared = 135.2268, df = 1, p-value < 2.2e-16

#####
#for Tc chr 9
#####

#for Tc chr 9/ Gr chr 4, blks 113,114,128
blk85s<-c("113","114","130")
chr85blks<-fracMatOrd[blk85s,]
chr_9_4_sums<-colSums(chr85blks[,1:3])

#for Tc chr 9/ Gr chr 9, blks 188,189
blk99s<-c("188","189")
chr99blks<-fracMatOrd[blk99s,]
chr_9_9_sums<-colSums(chr99blks[,1:3])

#for Tc chr 9/ Gr chr 13, blks 79,80,81,82
blk913s<-c("79","80","81","82")
chr913blks<-fracMatOrd[blk913s,]
chr_9_13_sums<-colSums(chr913blks[,1:3])


#compare 9 and 4

#[,1] [,2]
#Gr.retained  981  409
#Tc          3440 3176

chisq.test(cbind(chr_9_9_sums[2:1],chr_9_4_sums[2:1]))

#Pearson's Chi-squared test with Yates' continuity correction
#
#data:  cbind(chr_9_9_sums[2:1], chr_9_4_sums[2:1])
#X-squared = 159.6202, df = 1, p-value < 2.2e-16

#[,1] [,2]
#Gr.retained  981  400
#Tc          3440 2791

chisq.test(cbind(chr_9_9_sums[2:1],chr_9_13_sums[2:1]))

#Pearson's Chi-squared test with Yates' continuity correction

#data:  cbind(chr_9_9_sums[2:1], chr_9_13_sums[2:1])
#X-squared = 115.6658, df = 1, p-value < 2.2e-16

#[,1] [,2]
#Gr.retained  409  400
#Tc          3176 2791

chisq.test(cbind(chr_9_4_sums[2:1],chr_9_13_sums[2:1]))

#  Pearson's Chi-squared test with Yates' continuity correction

#data:  cbind(chr_9_4_sums[2:1], chr_9_13_sums[2:1])
#X-squared = 1.9324, df = 1, p-value = 0.1645

####
# For Tc chr 10
####

#for Tc chr 10/ Gr chr 9, blks 44
blk109s<-c("44")
chr109blks<-fracMatOrd[blk109s,]
chr_10_9_sums<-chr109blks[1:3]

#for Tc chr 10/ Gr chr 11, blks 195
blk1011s<-c("195")
chr1011blks<-fracMatOrd[blk1011s,]
chr_10_11_sums<-chr1011blks[1:3]

#[,1] [,2]
#Gr.retained  397  170
#Tc          1858 1792

chisq.test(cbind(chr_10_9_sums[2:1],chr_10_11_sums[2:1]))

#Pearson's Chi-squared test with Yates' continuity correction

#data:  cbind(chr_10_9_sums[2:1], chr_10_11_sums[2:1])
#X-squared = 71.2973, df = 1, p-value < 2.2e-16

#####
#Do the same chisq-test but with the total chromosome length, not just "overlap"
#####
###
#Chr 2
###
#Tc chr 2 Gr 5 8 
chisq.test(cbind(c(929,3641),c(642,3641)))
##
#Chr 6
##
#Tc chr 6 Gr 6 /
chisq.test(cbind(c(147,2637),c(580,2637)))
#Tc chr 6 Gr 6/10
chisq.test(cbind(c(147,2637),c(227,2637)))
#Tc chr 6 Gr 9/10
chisq.test(cbind(c(580,2637),c(227,2637)))

##
#Chr 7
##

#Tc chr 7 2/13
chisq.test(cbind(c(420,1873),c(225,1873)))

###
#Chr 8
###

#Tc chr 8 Gr 5/9
chisq.test(cbind(c(236,2040),c(608,2040)))

###
#Chr 9
###

#Tc chr 9 Gr 4/9
chisq.test(cbind(c(343,3599),c(981,3599)))

#Tc chr 9 Gr 9/13
chisq.test(cbind(c(400,3599),c(981,3599)))

#Tc chr 9 Gr 4/13
chisq.test(cbind(c(343,3599),c(409,3599)))

###
#Chr 10
###
#Tc chr 10 Gr 9/11

chisq.test(cbind(c(397,1873),c(170,1873)))


#####
# Do a "horse race" like Freeling and Wang
#####

#use only those genes in duplicate
#rownames where coverage == 2
depth2<-Tc_depth[Tc_depth[,"depth"] == 2,]
mapping_depth2<-mapping3[Tc_depth[,"depth"] == 2,]
PacData_depth2<-PacData3[Tc_depth[,"depth"] == 2,]

names<-t(apply(PacData_depth2,1,function(x) unique(x)))
#remove the "-"
names<-t(apply(names,1,function(x) x[which(x != "-")]))

#####
##Find the names that correspond to each pair
#####

GrNames<-apply(PacData_depth2,1,function(x) which(x != "-"))

#now for each row find out what two blocks are present.
pairs<-t(apply(PacData_depth2,1,function(x) which(x != "-")))

#now for each row find whcih the the most fractionated
#for each of the pair retrun the number which corresponds to the LF one of the pair
#i.e. if first blcok is 20% retaihned, and second blcok is 60% retained return the numer 2
LF<-apply(pairs,1,function(x) mostFrac(x,fracMatOrd))

#now grab the expression for each pair.
ListExp<-apply(names,1,function(x) returnExpr(x,expression_RPKM))

#now sort the pairs according to MF and LF, 1 is MF, 2 is LF
LF<-t(LF)
colnames(LF)<-c("MF","LF","pvalue")
MFexp<-NULL
LFexp<-NULL

for ( i in seq(1,length(LF[,1]))) {
  print(i)
  MFexp<-rbind(MFexp,ListExp[[i]][LF[i,1],])
  LFexp<-rbind(LFexp,ListExp[[i]][LF[i,2],])
}#for

MFexpMean<-t(apply(MFexp,1,function(x) expPerTissue(x)))
LFexpMean<-t(apply(LFexp,1,function(x) expPerTissue(x)))
colnames(MFexpMean)<-c("petal","leaf","seed")
colnames(LFexpMean)<-c("petal","leaf","seed")

#make a df to use in ggplot

expMFLF.df<-data.frame("race"=c(length(MFexpMean[,1][MFexpMean[,1]>LFexpMean[,1]]),
                                length(MFexpMean[,2][MFexpMean[,2]>LFexpMean[,2]]),
                                  length(MFexpMean[,3][MFexpMean[,3]>LFexpMean[,3]]),
                                    length(MFexpMean[,1][LFexpMean[,1]>MFexpMean[,1]]),
                                      length(MFexpMean[,2][LFexpMean[,2]>MFexpMean[,2]]),
                                        length(MFexpMean[,3][LFexpMean[,3]>MFexpMean[,3]])),
                        "tissue"=rep(c("petal", "leaf" , "seed"),2),
                            "frac"=c(rep("MF",3), rep("LF",3))
                              )
#####
# Do some high level plotting
#####

barPlot<-ggplot(data=expMFLF.df, aes(x=frac,weight=race, fill=frac)) +
    geom_bar() +
      facet_grid(~tissue)+
        xlab("") +
            theme(strip.text.x = element_text(size = 18, colour = "black")) +
                theme(axis.title.y = element_text(size = 18)) +
                  theme(axis.text.y = element_text(size = 18, angle = 0)) +
                    theme(axis.text.x = element_blank()) +
                      theme(legend.text = element_text(size=20),legend.key.size = unit(1.2, "cm"))+
                        theme(legend.title=element_blank())+
                          scale_fill_brewer(palette="Set1")

pdf("MFLFbarPlot.pdf",width=7,height=5)
barPlot
dev.off()

####
# A chi sq test to see if there is a difference between MF and LF
####

chisq.table<-matrix(ncol=3, nrow=2, expMFLF.df$race, byrow = T)
colSums(chisq.table)
#do the chisq by column

chiSqList<-apply(chisq.table,2,function(x) chisq.test(x))

m<-glm(log(as.vector(LFexpMean)+1)~log(as.vector(MFexpMean)+1))
MFLFmean.df<-cbind(melt(MFexpMean),melt(LFexpMean))
colnames(MFLFmean.df)<-c("geneMF","tissueMF","expMF","geneLF","tissueMF","expLF")

#1:1 line
one2one<-data.frame(x=c(-20:15), y = c(-20:15))

#a color palette
customColor<-colorRampPalette(c("yellow","orange", "red"))( 100 ) ## (n)

MFplot<-ggplot(data=MFLFmean.df, aes(x=log(expLF), y=log(expMF))) +
  geom_tile() +
    facet_grid(~tissueMF) +
        stat_density2d(aes(x=log(expLF), y=log(expMF), fill = ..density..), geom="tile", 
                  contour = FALSE)+
                    scale_fill_gradientn(colours=rainbow(100)) +
                      geom_line(data=one2one, aes(x=x, y = y), colour = "white") +
                        theme(strip.text.x = element_text(size = 20, colour = "black")) +
                          ylab("log(RPKM MF)")+
                          xlab("log(RPKM LF)")+
                            #scale_x_continuous(breaks=c(-5,0,5))+
                          theme(axis.title.y = element_text(size = 18)) +
                            theme(axis.title.x = element_text(size = 18)) +
                              theme(axis.text.y = element_text(size = 18, angle = 0)) +
                                theme(axis.text.x = element_text(size = 18, angle = 0)) +
                                  theme(legend.text = element_text(size=18),legend.key.size = unit(1.2, "cm"))+
                                    theme(legend.title=element_text(size=18)) +
                                      ylim(-3,7)+ 
                                        #xlim(-5,7.5) +
                                          scale_x_continuous(breaks = c(-5,0,5), limits=c(-5,7.5))

pdf("SampleGraph.pdf",width=7,height=5)
MFplot
dev.off()


####
#make a huge function to compare pairs of reconstructed chromosomes
#this prints a collection of figures for each chromosome comparison
####

#where x is the blcok numbers for one of the chromosomes, 
#and y is the block numbers of the others, 
#REMEMBER x IS LF AND y IS MF

compare_reconstruction<-function(x,y) {
  
###
# The horse race between MF and LF should be at the whole chromosome level
#i.e. those used in the chi sq test 
###

#assemble Tc chr 2 from cotton blocks

#using group 5 from cotton THIS IS THE LEAST FRACTIONATED OF THE CHROMOSOMES
LF.df<-as.matrix(PacData3[,x])
#print(dim(LF.df))
#using group 8 from cotton THIS IS THE MOST FRACTIONATED OF THE CHROMOSOMES
MF.df<-as.matrix(PacData3[,y])
#print(dim(MF.df))
#convert LF.df and MF.df into 1 column matrices according to some rule
#first sub in NA fir "-"
LF.df[LF.df == "-"]<-NA
MF.df[MF.df == "-"]<-NA
print
print(dim(MF.df))
print(dim(MF.df))

#if the dataframes have multiple columns
 if ( dim(LF.df)[2] > 1 ) {
      LF.df<-compress(as.matrix(LF.df))
      LF.df<-LF.df[,1]
}#if
else {
  LF.df<-as.matrix(LF.df[,1])
}
  if ( dim(MF.df)[2] >1 ) {
      MF.df<-compress(as.matrix(MF.df))
      MF.df<-MF.df[,1]
}#if
else {
  MF.df<-as.matrix(MF.df[,1])
}
#print(dim(MF.df))
#print(dim(LF.df))
#join the two chromosome copies to compare
TcChr2<-cbind(LF.df,MF.df)
TcChr2all<-cbind(LF.df,MF.df)
print(dim(TcChr2))
print(head(TcChr2))
#now return rows with all genes ( i.e. we can now compare the espression of these genes)
indexa<-which(!is.na(TcChr2[,1]))
indexb<-which(!is.na(TcChr2[,2]))
TcChr2allLF<-TcChr2all[indexa,1]
TcChr2allMF<-TcChr2all[indexb,2]
#find the intersection and extract the rows that intersect
if ( length(intersect(indexa,indexb)) < 10 ) {
  print(length(intersect(indexa,indexb)))
  stop ("fewerr than 10 genes overlap")
}#if
inter_indx<-intersect(indexa,indexb)
TcChr2<-TcChr2[inter_indx,]
#grab the coprresponding expression

#print(TcChr2[,1])
exprA<-expression_RPKM[as.character(TcChr2[,1]),]
exprB<-expression_RPKM[as.character(TcChr2[,2]),]

#binf the genes names into LF and MF

LF.genes<-LF.df[inter_indx]
MF.genes<-MF.df[inter_indx]


#calcualte average expressoion
chr2_5mean<-t(apply(exprA,1,function(x) expPerTissue(x)))
chr_2_8_mean<-t(apply(exprB,1,function(x) expPerTissue(x)))

colnames(chr2_5mean)<-c("petal","leaf","seed")
colnames(chr_2_8_mean)<-c("petal","leaf","seed")

exp2_8.df<-data.frame("race"=c(length(chr_2_8_mean[,1][chr_2_8_mean[,1]>chr2_5mean[,1]]),
                                length(chr_2_8_mean[,2][chr_2_8_mean[,2]>chr2_5mean[,2]]),
                                length(chr_2_8_mean[,3][chr_2_8_mean[,3]>chr2_5mean[,3]]),
                                length(chr_2_8_mean[,1][chr2_5mean[,1]>chr_2_8_mean[,1]]),
                                length(chr_2_8_mean[,2][chr2_5mean[,2]>chr_2_8_mean[,2]]),
                                length(chr_2_8_mean[,3][chr2_5mean[,3]>chr_2_8_mean[,3]])),
                       "tissue"=rep(c("petal", "leaf" , "seed"),2),
                       "frac"=c(rep("MF",3), rep("LF",3))
)

barPlot2<-ggplot(data=exp2_8.df, aes(x=frac,weight=race, fill=frac)) +
  geom_bar() +
  facet_grid(~tissue)+
  xlab("") +
  theme(strip.text.x = element_text(size = 18, colour = "black")) +
  theme(axis.title.y = element_text(size = 18)) +
  theme(axis.text.y = element_text(size = 18, angle = 0)) +
  theme(axis.text.x = element_blank()) +
  theme(legend.text = element_text(size=20),legend.key.size = unit(1.2, "cm"))+
  theme(legend.title=element_blank())+
  scale_fill_brewer(palette="Set1")

  myList<-list("plot"=barPlot2,"table"=exp2_8.df,"LF.genes"=LF.df[inter_indx],"MF.genes"=MF.df[inter_indx],"LF.exp"=exprA,"MF.exp"=exprB, "ALL.MF"=TcChr2allMF,"ALL.LF"=TcChr2allLF)
  #return ( barPlot2 )
  #return ( exp2_8.df )
  #return (LF.df[indexa])
  #return (MF.df[indexb])
  return(myList)
}#compare_reconstruction

####
#Actually use the above function
####
#asses and store the data
Tc2Gr5_8<-compare_reconstruction(c(137,138,139),c(179,184,185))# LF>MF

Tc6Gr6_9<-compare_reconstruction(c(190),c(149,150))# LF>MF for 2/3 tissues

Tc6Gr6_10<-compare_reconstruction(c(33,34,36),c(149,150))# MF>LF 2/3 tissue

Tc6Gr9_10<-compare_reconstruction(c(190),c(33,34,36))# LF>MF

Tc7Gr2_13<-compare_reconstruction(c(88,86,89),c(76,75))# LF>MF

Tc8Gr5_9<-compare_reconstruction(c(191),c(133,132))# LF>MF

Tc9Gr9_4<-compare_reconstruction(c(188,189),c(113,130))# LF>MF

Tc9Gr9_13<-compare_reconstruction(c(188,189),c(79,80,81,82))# LF>MF

Tc9Gr4_13<-compare_reconstruction(c(113,130),c(79,80,81,82))# MF>LF

Tc10Gr9_11<-compare_reconstruction(c(195),c(44))#MF>LF 2/3

#
#A complete list of comparable genes in LF and MF, for DGE analysis
#

MF.dge.names<-c(Tc2Gr5_8$MF.genes,Tc6Gr6_9$MF.genes,Tc6Gr6_10$MF.genes,Tc6Gr9_10$MF.genes,Tc7Gr2_13$MF.genes
                ,Tc8Gr5_9$MF.genes,Tc9Gr9_4$MF.genes,Tc9Gr9_13$MF.genes,Tc9Gr4_13$MF.genes,Tc10Gr9_11$MF.genes)

LF.dge.names<-c(Tc2Gr5_8$LF.genes,Tc6Gr6_9$LF.genes,Tc6Gr6_10$LF.genes,Tc6Gr9_10$LF.genes,Tc7Gr2_13$LF.genes
                ,Tc8Gr5_9$LF.genes,Tc9Gr9_4$LF.genes,Tc9Gr9_13$LF.genes,Tc9Gr4_13$LF.genes,Tc10Gr9_11$LF.genes)

LF.exp.dge<-expression_RPKM[LF.dge.names,]
MF.exp.dge<-expression_RPKM[MF.dge.names,]

exp.dge<-cbind(LF.exp.dge,MF.exp.dge)

save(exp.dge, file="LFvsMF.dge.df.R")

#list the matrices first

listMt<-list("Tc2Gr5_8"=Tc2Gr5_8$table,"Tc6Gr6_9"=Tc6Gr6_9$table,"Tc6Gr6_10"=Tc6Gr6_10$table,"Tc6Gr9_10"=Tc6Gr9_10$table,
             "Tc7Gr2_13"=Tc7Gr2_13$table,"Tc8Gr5_9"=Tc8Gr5_9$table,"Tc9Gr9_4"=Tc9Gr9_4$table,"Tc9Gr9_13"=Tc9Gr9_13$table,"Tc9Gr4_13"=Tc9Gr4_13$table,"Tc10Gr9_11"=Tc10Gr9_11$table)
totalTable<-sapply(listMt,function(x) print(x[,1]))
sumTable<-data.frame("race"=as.numeric(rowSums(totalTable)),"tissue"=as.character(Tc2Gr5_8$table[,2]),"frac"=as.character(Tc2Gr5_8$table[,3]))

#make a data.frame, so as to make a facet plot by chromosome comparison
completeRaceTable<-rbind(Tc2Gr5_8$table,Tc6Gr6_9$table,Tc6Gr6_10$table,Tc6Gr9_10$table,
             Tc7Gr2_13$table,Tc8Gr5_9$table,Tc9Gr9_4$table,Tc9Gr9_13$table,
              Tc9Gr4_13$table,Tc10Gr9_11$table)

completeRaceTable<-cbind(completeRaceTable,"comparison"=c(rep("2:Gr5/Gr8",6),rep("6:Gr6/Gr9",6),rep("6:Gr6/Gr10",6),
                                             rep("6:Gr9/Gr10",6),rep("7:Gr2/Gr13",6),rep("8:Gr5/Gr9",6),
                                             rep("9:Gr9/Gr4",6),rep("9:Gr9/Gr13",6),rep("9:Gr4/Gr13",6),
                                             rep("10:Gr9/Gr11",6)))

facetBarPlot<-ggplot(data=completeRaceTable, aes(x=frac,weight=race, fill=frac)) +
  geom_bar() +
  #facet_wrap(~comparison) +
  facet_wrap(comparison~tissue,scales="free_y")+
  xlab("") +
  ylab("count (number of wins)") +
  theme(strip.text.x = element_text(size = 8, colour = "black")) +
  theme(axis.title.y = element_text(size = 18)) +
  theme(axis.text.y = element_text(size = 11, angle = 0)) +
  theme(axis.text.x = element_blank()) +
  theme(legend.text = element_text(size=18),legend.key.size = unit(1.2, "cm"))+
  theme(legend.title=element_blank())+
  scale_fill_brewer(palette="Set1")

pdf("sep_chr_race.pdf", width = 10 , height = 7)
facetBarPlot
dev.off()

#make a barplot of all the MF and LF races
totalBarPlot<-ggplot(data=sumTable, aes(x=frac,weight=race, fill=frac)) +
  geom_bar() +
  facet_grid(~tissue)+
  xlab("") +
  ylab("count (number of wins)") +
  theme(strip.text.x = element_text(size = 20, colour = "black")) +
  theme(axis.title.y = element_text(size = 26)) +
  theme(axis.text.y = element_text(size = 22, angle = 90, colour = "black")) +
  theme(axis.text.x = element_blank()) +
  theme(legend.text = element_text(size=20),legend.key.size = unit(1.2, "cm"))+
  theme(legend.title=element_blank())+
  scale_fill_brewer(palette="Set1")

pdf("total_bar_plot.pdf", height = 5, width =5)
totalBarPlot
dev.off()

###
# test using a cumulative binomial
###

pbinom(559, size=1064, p = 0.5)

####

all.LF.genes<-c("Tc2Gr5_8"=Tc2Gr5_8$ALL.LF,"Tc6Gr6_9"=Tc6Gr6_9$ALL.LF,"Tc7Gr2_13"=Tc7Gr2_13$ALL.LF,
                "Tc8Gr5_9"=Tc8Gr5_9$ALL.LF,"Tc9Gr9_4"=Tc9Gr9_4$ALL.LF,"Tc10Gr9_11"=Tc10Gr9_11$ALL.LF)

all.MF.genes<-c("Tc2Gr5_8"=Tc2Gr5_8$ALL.MF,"Tc6Gr6_9"=Tc6Gr6_9$ALL.MF,"Tc6Gr6_10"=Tc6Gr6_10$ALL.MF,
                "Tc7Gr2_13"=Tc7Gr2_13$ALL.MF,"Tc8Gr5_9"=Tc8Gr5_9$ALL.MF,"Tc9Gr9_4"=Tc9Gr9_4$ALL.MF,"Tc9Gr9_13"=Tc9Gr9_13$ALL.MF,"Tc10Gr4_11"=Tc10Gr9_11$ALL.MF)

#remove the middle fractionated block when there are mnore than two
all.LF.exp<-rbind("Tc2Gr5_8"=Tc2Gr5_8$LF.exp,"Tc6Gr6_9"=Tc6Gr6_9$LF.exp,"Tc7Gr2_13"=Tc7Gr2_13$LF.exp,
                "Tc8Gr5_9"=Tc8Gr5_9$LF.exp,"Tc9Gr9_4"=Tc9Gr9_4$LF.exp,"Tc10Gr9_11"=Tc10Gr9_11$LF.exp)

all.MF.exp<-rbind("Tc2Gr5_8"=Tc2Gr5_8$MF.exp,"Tc6Gr6_9"=Tc6Gr6_9$MF.exp, "Tc7Gr2_13"=Tc7Gr2_13$MF.exp,
              "Tc8Gr5_9"=Tc8Gr5_9$MF.exp,"Tc9Gr9_4"=Tc9Gr9_4$MF.exp,"Tc10Gr9_11"=Tc10Gr9_11$MF.exp)
#calculate expression per tissue
LF.exp.mean<-t(apply(all.LF.exp,1,function(x) expPerTissue(x)))
MF.exp.mean<-t(apply(all.MF.exp,1,function(x) expPerTissue(x)))
colnames(MF.exp.mean)<-c("petal","leaf","seed")
colnames(LF.exp.mean)<-c("petal","leaf","seed")

#melt the data
MFLFexp.df<-cbind(melt(MF.exp.mean),melt(LF.exp.mean))
colnames(MFLFexp.df)<-c("geneMF","tissueMF","expMF","geneLF","tissueMF","expLF")

MFplot2<-ggplot(data=MFLFexp.df, aes(x=log(expLF), y=log(expMF))) +
  geom_tile() +
  facet_grid(~tissueMF) +
  stat_density2d(aes(x=log(expMF), y=log(expLF), fill = ..density..), geom="tile", 
                 contour = FALSE)+
  scale_fill_gradientn(colours=rainbow(100)) +
  geom_line(data=one2one, aes(x=x, y = y), colour = "white") +
  theme(strip.text.x = element_text(size = 20, colour = "black")) +
  ylab("log(RPKM LF)")+
  xlab("log(RPKM MF)")+
  #scale_x_continuous(breaks=c(-5,0,5))+
  theme(axis.title.y = element_text(size = 18)) +
  theme(axis.title.x = element_text(size = 18)) +
  theme(axis.text.y = element_text(size = 18, angle = 0)) +
  theme(axis.text.x = element_text(size = 18, angle = 0)) +
  theme(legend.text = element_text(size=18),legend.key.size = unit(1.2, "cm"))+
  theme(legend.title=element_text(size=18)) +
  ylim(-3,7)+ 
  #xlim(-5,7.5) +
  scale_x_continuous(breaks = c(-5,0,5), limits=c(-5,7.5))
pdf("SampleGraph.pdf",width=7,height=5)
MFplot2
dev.off()

###
###############################################################################
#siRNA analysis
###############################################################################
###

####
#re-arrange the data to be in a logical order

sRNA_data2<-rbind(sRNA_data[502:1002,],sRNA_data[1:501,])
sRNA_data2[sRNA_data2$dir == "up",2]<-rev(-abs(sRNA_data2[sRNA_data2$dir == "up",2]))

#draw plot for LF and then overlay LF
LF_sRNA_Ave2<-apply(sRNA_data2[,all.LF.genes], 1 , function(x) mean(x))

LF_sRNA_Ave2<-data.frame("Ave"=LF_sRNA_Ave2,"window"=sRNA_data2[,2]*10, 
                            "ori" =sRNA_data2[,1])

sRNA_overall1000<-apply(sRNA_data2[c(450:550),-c(1:2)],2,function(x) mean(x))

MF_sRNA_Ave2<-apply(sRNA_data2[,all.MF.genes], 1 , function(x) mean(x))

MF_sRNA_Ave2<-data.frame("Ave"=MF_sRNA_Ave2,"window"=sRNA_data2[,2]*10, 
                        "ori" =sRNA_data2[,1])

global_sRNA_Ave2<-apply(sRNA_data2[,3:length(sRNA_data2[1,])], 1 , function(x) mean(x))

global_sRNA_Ave2<-data.frame("Ave"=global_sRNA_Ave2,"window"=sRNA_data2[,2]*10, 
                            "ori"=sRNA_data2[,1])

#modify the columns names in the exp data

rownames(MF.exp.mean)<-sub(".*\\.Gorai", "Gorai" , rownames(MF.exp.mean), perl = TRUE )
rownames(LF.exp.mean)<-sub(".*\\.Gorai", "Gorai" , rownames(LF.exp.mean), perl = TRUE )

#add the average 24 nt coverage to the dataframe

ntExpLF<-data.frame("exp"=LF.exp.mean+1 , "24nt" = sRNA_overall1000[rownames(LF.exp.mean)]+1 , "Frac" = "LF" )
ntExpMF<-data.frame("exp"=MF.exp.mean+1 , "24nt" = sRNA_overall1000[rownames(MF.exp.mean)]+1, "Frac" = "MF" )

ntExp<-rbind(ntExpLF,ntExpMF)

ntVsExpLF<-ggplot(data=ntExp, aes(x=log2(X24nt), y=log2(exp.leaf), colour = factor(Frac))) +
  geom_point(size=3) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x, size = 2)

ntVsExpLF<-ggplot(data=ntExp, aes(x=log2(X24nt), y=log2(exp.leaf), colour = factor(Frac))) +
  geom_point(size=3) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x, size = 2)

histPlot<-ggplot(data=ntExp, aes(x=log10(X24nt), fill = factor(Frac))) +
  geom_density(alpha=0.5)

ggplot(data=ntExp, aes(x=log(X24nt[Frac=="LF"]), y = log(X24nt[Frac=="MF"]))) + 
  geom_point()+
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x, size = 2)

###
# A global data set of expression and sRNA abundance
###
#get the mean expression of all genes
globalExpMean<-t(apply(expression_RPKM,1,function(x) expPerTissue(x)))
colnames(globalExpMean)<-c("petal","leaf","seed")
globalsiRNAExp<-cbind(globalExpMean+1,"siRNA"=sRNA_overall1000[rownames(globalExpMean)]+1)

#make the data frame to plot sRNAs

LF_sRNA_Ave2<-cbind(LF_sRNA_Ave2,"frac"="LF")
MF_sRNA_Ave2<-cbind(MF_sRNA_Ave2,"frac"="MF")

sRNA.dt<-rbind(LF_sRNA_Ave2,MF_sRNA_Ave2)
Palette1 <- c('red','blue')

smallRNAplotLFMF<-ggplot(data=sRNA.dt, aes(x=window,y=Ave, factor(frac))) +
    geom_line(aes(y=Ave, x = window,color=factor(frac)),size = 1.1) +
        geom_vline(xintercept = 0, linetype = "longdash" , size=1.2) +
          scale_colour_manual(values=Palette1) +
            ylab("mean number of mapped reads") +
              theme(axis.title.x = element_text(size = 18)) +
                theme(axis.title.y = element_text(size = 18)) +
                  xlab("distance from transcription start/stop site (bp)") +
                    theme(axis.text.y = element_text(size = 18, angle = 0)) +
                      theme(axis.text.x = element_text(size = 18, angle = 0)) +
                        scale_x_continuous(limits=c(-1250,1250)) +
                          theme(legend.title=element_blank()) +
                            theme(legend.key.size = unit(1.3, "cm")) +
                              theme(legend.text=element_text(size = 14))

pdf("sRNAmapping_no_masking_zoom.pdf", height = 5 , width = 8)
smallRNAplotLFMF
dev.off()

###
#Anova for siRNA, LF/MF and distance to start/stop site
###

summary(lm(sRNA.dt$Ave~sRNA.dt$frac*abs(sRNA.dt$window)))

anova(lm(sRNA.dt$Ave~sRNA.dt$frac*abs(sRNA.dt$window)))

wilcox.test(sRNA.dt[,1][sRNA.dt$frac=="MF"],sRNA.dt[,1][sRNA.dt$frac=="LF"], paired = TRUE, exact = TRUE)

#####
#examine the distribution of TE insertions in LF and MF
#####

rownames(TEpos)<-TEpos[,1]

hist(TEpos[as.character(all.LF.genes),2], breaks = 100)
hist(TEpos[as.character(all.MF.genes),2], breaks = 100)

TE.df<-data.frame("pos"=c(TEpos[as.character(all.MF.genes),2],TEpos[as.character(all.LF.genes),2]), 
                    "frac" = c(rep("MF",length(TEpos[as.character(all.MF.genes),2])),rep("LF",length(TEpos[as.character(all.LF.genes),2]))))

TEdensityplot<-ggplot(data=TE.df, aes(x=pos, col = frac)) +
  scale_colour_manual(values=Palette1) +
  geom_density(subset = .(frac == 'MF'),aes(fill=frac),size = 1.2, alpha=0.5) +
    geom_density(subset = .(frac == 'LF'),aes(fill=frac),size = 1.2, alpha=0.5) +
      xlim(-2000,2000) +
        xlab("distance from transcription stop/start site (bp)") +
          ylab("kernal density - nearest TE") +
            theme(axis.title.x = element_text(size = 18)) +
              theme(axis.title.y = element_text(size = 18)) +
theme(axis.text.y = element_text(size = 18, angle = 0)) +
  theme(axis.text.x = element_text(size = 18, angle = 0)) +
  theme(legend.title=element_blank()) +
  theme(legend.key.size = unit(1.3, "cm")) +
  theme(legend.text=element_text(size = 14)) +
  
pdf("nearest_TE.pdf", width = 8 , height = 5)
TEdensityplot
dev.off()

###
# Scale the RNA seq data like Hollister et al., 2011
###

meanExp<-cbind(petalAve,leafAve,seedAve)
#does not look very normal
LogExpression<-as.matrix(log(meanExp))
scaleExpression<-scale(LogExpression,center = TRUE, scale = TRUE)
TEpos<-TEpos[sort.list(TEpos[,1]),]
#make the breaks
#breaks<-pretty(c(-1000:5000))
breaks<-c(-1000, 0 , 1 , 1000, seq(2000,17000, by =1000))
f<-cut(abs(TEpos$distance), breaks=breaks,include.lowest=FALSE)

scaled_exp_by_bin<-as.data.frame(cbind(scaleExpression, as.factor(f)))
scaled_df<-melt(scaled_exp_by_bin[,1:3])
scaled_df<-cbind(scaled_df,f)

###
# check expression by nearest TE
###

TEplot<-ggplot(data=scaled_df, aes(y=value, factor(f)))+
      geom_boxplot()+
        facet_wrap(~variable)

##make a plot TEs for within genes and targeted sRNAs

as.data.frame(withinTE.dt)
withinTE.dt<-data.frame(withinTE.dt)
rownames(withinTE.dt)<-withinTE.dt$gene

within_SRNA.df<-data.frame("RNA"=c(withinTE.dt[as.character(all.MF.genes),2],withinTE.dt[as.character(all.LF.genes),2]), 
                  "frac" = c(rep("MF",length(withinTE.dt[as.character(all.MF.genes),2])),rep("LF",length(withinTE.dt[as.character(all.LF.genes),2]))))

ggplot(data=within_SRNA.df , aes(y=log(within_SRNA.df$RNA+1) , factor(frac), fill= frac)) + 
      geom_boxplot() +
        geom_point(position=position_jitter(width=0.3), alpha=0.2, size = 4) 

ggplot(data=within_SRNA.df , aes(x=log(within_SRNA.df$RNA+1), fill=factor(frac))) +
      geom_density(alpha = 0.4)

###
#plot TE density over 100bp window
###

#first grab all the MF and LF genes from the TE density 

LFindex<-which( TEcoverage[,1] %in% as.character(all.LF.genes) )
TEdensityLF<-TEcoverage[LFindex,]

MFindex<-which( TEcoverage[,1] %in% as.character(all.MF.genes) )
TEdensityMF<-TEcoverage[MFindex,]

TEdensityMF2<-cbind(TEdensityMF,"MF")
TEdensityLF2<-cbind(TEdensityLF,"LF")
colnames(TEdensityMF2)<-NULL
colnames(TEdensityLF2)<-NULL
rownames(TEdensityMF2)<-NULL
rownames(TEdensityLF2)<-NULL
TEdensityALL<-rbind(data.frame(TEdensityMF2),data.frame(TEdensityLF2))
TEdensityALL<-TEdensityALL[,c(1,4,8,9)]

colnames(TEdensityALL)<-c("gene","window","density","frac")
MFLF.lm<-lm(TEdensityALL$density~TEdensityALL$frac*abs(TEdensityALL$window+1))
summary(MFLF.lm)
anova(MFLF.lm)

wilcox.test(TEdensityAve[,1][TEdensityAve$frac=="MF"],TEdensityAve[,1][TEdensityAve$frac=="LF"], paired = TRUE, eaxct = TRUE)


#for each window calculate the average
TEdensityAveMF<-tapply(TEdensityMF[,8], TEdensityMF[,4], mean)
TEdensityAveLF<-tapply(TEdensityLF[,8], TEdensityLF[,4], mean)
TEdensityAve<-data.frame("density"=c(as.numeric(TEdensityAveLF),as.numeric(TEdensityAveMF)),"frac"=c(rep("LF",length(as.numeric(TEdensityAveMF))),rep("MF",length(as.numeric(TEdensityAveMF)))),"window"=as.numeric(rownames(TEdensityAveLF)))

Palette1 <- c('red','blue')

TEdensityPlot<-ggplot(data=TEdensityAve, aes(x=window*10, y=density, factor(frac), color=frac)) +
  geom_line(size = 1.1)+
  scale_colour_manual(values=Palette1) +
  ylab("TE denstiy (proportion of bp)") +
  theme(axis.title.x = element_text(size = 18)) +
  theme(axis.title.y = element_text(size = 18)) +
  xlab("distance from transcription start/stop site (bp)") +
  theme(axis.text.y = element_text(size = 18, angle = 0)) +
  theme(axis.text.x = element_text(size = 18, angle = 0)) +
  theme(legend.title=element_blank()) +
  theme(legend.key.size = unit(1.3, "cm")) +
  theme(legend.text=element_text(size = 14)) +
  geom_vline(xintercept = 0, linetype = "longdash" , size=1.2)

pdf("TE_density.pdf", height = 5 , width = 8)
TEdensityPlot
dev.off()

###
#An anova of TEdesnity in LF and Mf
###

lm(TEdensityAve$density~TEdensityAve$frac*abs(TEdensityAve$window))

#Call:
#  lm(formula = TEdensityAve$density ~ TEdensityAve$frac * abs(TEdensityAve$window))

#Residuals:
#  Min       1Q   Median       3Q      Max 
#-0.18874 -0.03111 -0.01667  0.04326  0.07411 

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)                                  1.911e-01  2.774e-03  68.875  < 2e-16 ***
#  TEdensityAve$fracMF                          7.611e-02  3.924e-03  19.397  < 2e-16 ***
#  abs(TEdensityAve$window)                     3.178e-04  9.578e-06  33.177  < 2e-16 ***
#  TEdensityAve$fracMF:abs(TEdensityAve$window) 7.939e-05  1.354e-05   5.861 5.37e-09 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 0.04385 on 2000 degrees of freedom
#Multiple R-squared:  0.7231,  Adjusted R-squared:  0.7227 
#F-statistic:  1741 on 3 and 2000 DF,  p-value: < 2.2e-16


dim()
 
##
#make a plot of GC
##
 
LF_gc<-gcTable[all.LF.genes,]
MF_gc<-gcTable[all.MF.genes,]
MF_gc<-na.omit(MF_gc)
LF_gc<-na.omit(LF_gc)
 
#transpose table and calculate row means ( i.e. mean GC over all genes in the list)
 LF_gc<-t(LF_gc)
 MF_gc<-t(MF_gc)
 rownames(MF_gc)<-seq(-5000, 4990, by = 10)
rownames(LF_gc)<-seq(-5000, 4990, by = 10)

all_MF_gc<-melt(MF_gc)
all_LF_gc<-melt(LF_gc)

all_LF_gc<-cbind(all_LF_gc,"LF")
all_MF_gc<-cbind(all_MF_gc,"MF")

colnames(all_LF_gc)<-c("window","gene","gc","frac")
colnames(all_MF_gc)<-c("window","gene","gc","frac")

all_gc_lm<-rbind(all_LF_gc,all_MF_gc)

###
#Anova on gc, predicted by LF/MF and window
###

test.lm<-lm(all_gc_lm$gc~all_gc_lm$frac*log(abs(all_gc_lm$window)+1))

anova(test.lm)

#Analysis of Variance Table
#
#Response: all_gc_lm$gc
#                                          Df Sum Sq Mean Sq F value    Pr(>F)    
#all_gc_lm$frac                             1     91   90.89  9074.3 < 2.2e-16 ***
#  abs(all_gc_lm$window)                      1    573  572.88 57196.7 < 2.2e-16 ***
#  all_gc_lm$frac:abs(all_gc_lm$window)       1     10   10.24  1021.9 < 2.2e-16 ***
#  Residuals                            6242996  62530    0.01                      
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


LF_gc_mean<-rowMeans(LF_gc) 
MF_gc_mean<-rowMeans(MF_gc)
 
 
gc_df<-data.frame("gc"=c(LF_gc_mean, MF_gc_mean),"frac" = c(rep("LF", length(LF_gc_mean)), rep("MF",length(MF_gc_mean) )), "window" = seq(-5000, 4990, by = 10)) 
 
gc_plot<-ggplot(data=gc_df, aes( x = window, y = gc , factor(frac) , color = frac )) +
  geom_line() +
   scale_colour_manual(values=Palette1) +
   ylab("GC content") +
   theme(axis.title.x = element_text(size = 18)) +
   theme(axis.title.y = element_text(size = 18)) +
   xlab("distance from transcription start/stop site (bp)") +
   theme(axis.text.y = element_text(size = 18, angle = 0)) +
   theme(axis.text.x = element_text(size = 18, angle = 0)) +
   theme(legend.title=element_blank()) +
   theme(legend.key.size = unit(1.3, "cm")) +
   theme(legend.text=element_text(size = 14)) +
   geom_vline(xintercept = 0, linetype = "longdash" , size=1.2)

 pdf("gc_plot.pdf", height = 5 , width = 8)
 gc_plot
 dev.off()

###
# Anova of distance to start/stop and LF/MF on gc content
###

summary(lm(gc_df$gc~gc_df$frac*abs(gc_df$window)))

anova(lm(gc_df$gc~gc_df$frac*abs(gc_df$window)))

wilcox.test(gc_df[,1][gc_df$frac=="MF"],gc_df[,1][gc_df$frac=="LF"], paired = TRUE, eaxct = TRUE)
 
#########
#Now look at the position of the TEs within each gene region
#########
 
region<-read.table("/Users/simonrenny-byfield/cotton/TE_Gene_interaction/TE_by_region/TE_by_region.txt", header = T , sep = "\t" )
 
#modify the table a littel so it makes more sense.
 
#name the rows
#rownames(region)<-region[,1]
#remove row one (now rownames)
#region<-region[,-1]
#remove the bits we don't need
cols<-colnames(region)
 region<-region[,-8]
 colnames(region)<-cols[2:13]
#calculate percentages
percent_Region<-(region[,1:6]/region[,7:12])*100 
#remove NaNs
 percent_Region[is.na(percent_Region)]<-0
 
 #exclude genes that overlap with other genes
percent_Region<-as.matrix(percent_Region[rowMax(as.matrix(percent_Region)) <= 100,])
#melt the matrix into a data.frame for ggplot
percent_region.df<-melt(percent_Region)
 percent_region.df<-data.frame("region" = percent_region.df$X2 , "percent" = as.numeric(percent_region.df$value) )

regiopPlot<-ggplot(data=percent_region.df, aes(y=percent, x = region, col = region, fill =region)) +
   geom_boxplot(alpha = 0.5,outlier.size = 0.3, outlier.colour=NULL,size=1.5)
 
 regiopDensity<-ggplot(data=percent_region.df, aes(x=percent, col = region)) +
   geom_density(alpha = 0.4) +
    xlim(0,2.5)
 
 #now what about LF and MF
 
 perLFMF.df<-melt(percent_Region)
 
 perMF.df<-perLFMF.df[perLFMF.df$X1 %in% all.MF.genes,]
 perLF.df<-perLFMF.df[perLFMF.df$X1 %in% all.LF.genes,]
 
 perMF.df<-data.frame("region" = perMF.df$X2 , "percent" = as.numeric(perMF.df$value) )
 perLF.df<-data.frame("region" = perLF.df$X2 , "percent" = as.numeric(perLF.df$value) )
 
 ggplot(data=perMF.df, aes(y=percent, x = region, col = region, fill =region)) +
   geom_boxplot(alpha = 0.5,outlier.size = 0.3, outlier.colour=NULL,size=1.5)
 ggplot(data=perLF.df, aes(y=percent, x = region, col = region, fill =region)) +
   geom_boxplot(alpha = 0.5,outlier.size = 0.3, outlier.colour=NULL,size=1.5)
 
######
#End
######
