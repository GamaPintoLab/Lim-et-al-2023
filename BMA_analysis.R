# -- BMA ANALYSIS -------------------------------------------------------------------------------

library(BMA)

#Import table with expression values for the 16 differential expressed miRNA

de_miRs=read.table("de_miRs.txt",sep="\t",header=T, row.names = 1)

# create vector with responders or non-responders (1/0)
resp=de_miRs$GR 

#removing the column for Good/Non-responders
mir=data.frame(de_miRs[1:15]) 

res_sum=summary(bic.glm(mir,resp,glm.family = "binomial"))
par(mar = c(2,2,2,2))
imageplot.bma(bic.glm(mir,resp,glm.family="binomial"))