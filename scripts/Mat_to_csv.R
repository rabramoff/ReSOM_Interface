#Load interface data from matlab and merge into long-format dataset
#Columns: Year, total SOC, total CO2 production rate, mineral associated DOC and enzymes, free DOC and enzymes, polymers, microbial structural biomass, microbial reserve biomass
library(R.matlab)

path <- "/Users/rzabramoff/Dropbox (Climate)/External Collaborations/Interface/code_scripts/output/"
surfinit = c(6029.120, 7332.271, 8445.413)
litterQ = round(c(2.4133-2.4133*0.43, 2.4133-2.4133*0.1), digits = 4)
expt = c(0, 1, 2, 3, 4, 5, 6, 7, 8) 

#import each pathname as an item in a list
pathname = NULL
df = list(NULL)
lenmat <- 548
a <- rep("low clay",lenmat)
b <- rep("mid clay",lenmat)
c <- rep("low quality litter",lenmat)
d <- rep("high quality litter",lenmat)
e <- rep("Control",lenmat)
f <- rep("+2C",lenmat)
g <- rep("+5C",lenmat)
h <- rep("+30% total input",lenmat)
m <- rep("+100% total input",lenmat)
n <- rep("+30% labile input",lenmat)
o <- rep("+2C and +30% total input",lenmat)
p <- rep("+2C and +30% labile input",lenmat)
q <- rep("No input",lenmat)
r <- rep("high clay",lenmat) ##need to add

#import matfiles and add labels
for (i in 1:length(surfinit)){
  for (j in 1:length(litterQ)){
    for (k in 1:length(expt)){
      pathname[length(expt)*2*(i-1)+length(expt)*(j-1)+k] <- file.path(path, paste0("mbms_noIso_",surfinit[i],"surfinit_",litterQ[j],"litterQ_",expt[k],"expt.mat"))
      df[length(expt)*2*(i-1)+length(expt)*(j-1)+k] = readMat(pathname[length(expt)*2*(i-1)+length(expt)*(j-1)+k])[2]
      
      df[length(expt)*2*(i-1)+length(expt)*(j-1)+k][[1]][,9] =  c(diff(df[length(expt)*2*(i-1)+length(expt)*(j-1)+k][[1]][,9]/40),NA)
      
      if (surfinit[i] == 6029.120){
        df[length(expt)*2*(i-1)+length(expt)*(j-1)+k][[1]] = data.frame(df[length(expt)*2*(i-1)+length(expt)*(j-1)+k][[1]],a)
      } else if (surfinit[i] == 7516.784) {
        df[length(expt)*2*(i-1)+length(expt)*(j-1)+k][[1]] = data.frame(df[length(expt)*2*(i-1)+length(expt)*(j-1)+k][[1]],b)
      } else {  
        df[length(expt)*2*(i-1)+length(expt)*(j-1)+k][[1]] = data.frame(df[length(expt)*2*(i-1)+length(expt)*(j-1)+k][[1]],r)
      }
      
      if (litterQ[j] == 1.3756){
        df[length(expt)*2*(i-1)+length(expt)*(j-1)+k][[1]] = data.frame(df[length(expt)*2*(i-1)+length(expt)*(j-1)+k][[1]],c) #CHECK DIS
      } else {
        df[length(expt)*2*(i-1)+length(expt)*(j-1)+k][[1]] = data.frame(df[length(expt)*2*(i-1)+length(expt)*(j-1)+k][[1]],d)
      }
      
      if (expt[k] == 0) { 
        df[length(expt)*2*(i-1)+length(expt)*(j-1)+k][[1]] = data.frame(df[length(expt)*2*(i-1)+length(expt)*(j-1)+k][[1]],e)
      } else if (expt[k] == 1) {
        df[length(expt)*2*(i-1)+length(expt)*(j-1)+k][[1]] = data.frame(df[length(expt)*2*(i-1)+length(expt)*(j-1)+k][[1]],f)
      } else if (expt[k] == 2) {
        df[length(expt)*2*(i-1)+length(expt)*(j-1)+k][[1]] = data.frame(df[length(expt)*2*(i-1)+length(expt)*(j-1)+k][[1]],g)
      } else if (expt[k] == 3) {
        df[length(expt)*2*(i-1)+length(expt)*(j-1)+k][[1]] = data.frame(df[length(expt)*2*(i-1)+length(expt)*(j-1)+k][[1]],h)
      } else if (expt[k] == 4) {
        df[length(expt)*2*(i-1)+length(expt)*(j-1)+k][[1]] = data.frame(df[length(expt)*2*(i-1)+length(expt)*(j-1)+k][[1]],m)
      } else if (expt[k] == 5) {
        df[length(expt)*2*(i-1)+length(expt)*(j-1)+k][[1]] = data.frame(df[length(expt)*2*(i-1)+length(expt)*(j-1)+k][[1]],n)
      } else if (expt[k] == 6) {
        df[length(expt)*2*(i-1)+length(expt)*(j-1)+k][[1]] = data.frame(df[length(expt)*2*(i-1)+length(expt)*(j-1)+k][[1]],o)
      } else if (expt[k] == 7) {
        df[length(expt)*2*(i-1)+length(expt)*(j-1)+k][[1]] = data.frame(df[length(expt)*2*(i-1)+length(expt)*(j-1)+k][[1]],p)
      } else {
        df[length(expt)*2*(i-1)+length(expt)*(j-1)+k][[1]] = data.frame(df[length(expt)*2*(i-1)+length(expt)*(j-1)+k][[1]],q)
      } 
    }
  }
}

thenames <- names(as.data.frame(df[1]))
newdf = list(NULL)
for (i in 1:length(df)){
  newdf[[i]] <- df[[i]]
  newdf[[i]] <- setNames(newdf[[i]], thenames)
}

bigdf <- rbind(newdf[[1]],newdf[[2]],newdf[[3]],newdf[[4]],newdf[[5]],newdf[[6]],newdf[[7]],newdf[[8]],newdf[[9]],newdf[[10]],newdf[[11]],newdf[[12]],newdf[[13]],newdf[[14]],newdf[[15]],newdf[[16]],newdf[[17]],newdf[[18]],newdf[[19]],newdf[[20]],newdf[[21]],newdf[[22]],newdf[[23]],newdf[[24]],newdf[[25]],newdf[[26]],newdf[[27]],newdf[[28]],newdf[[29]],newdf[[30]],newdf[[31]],newdf[[32]],newdf[[33]],newdf[[34]],newdf[[35]],newdf[[36]],newdf[[37]],newdf[[38]],newdf[[39]],newdf[[40]],newdf[[41]],newdf[[42]],newdf[[43]],newdf[[44]],newdf[[45]],newdf[[46]],newdf[[47]],newdf[[48]],newdf[[49]],newdf[[50]],newdf[[51]],newdf[[52]],newdf[[53]],newdf[[54]])

bigdf$Year <- rep(seq(0,60,length.out=lenmat),54)
conv = 0.4; #gC m-3 to gC m-2 to 20cm ##proportion of SOC in top 20cm
colnames(bigdf) <- c("Microbial structural biomass (gC m-2)", "Microbial reserve biomass (gC m-2)", "Available mineral sites (gC eqv m-2)", "Free monomers  (gC m-2)", "Adsorbed monomers (gC m-2)", "Polymers (gC m-2)", "Free enzymes (gC m-2)", "Adsorbed enzymes (gC m-2)", "CO2 production (gC m-2 d-1)","CUE","defactoTurnover","Clay content","Litter quality","Experiment","Year")

bigdf[,c("Microbial structural biomass (gC m-2)", "Microbial reserve biomass (gC m-2)","Free monomers  (gC m-2)", "Adsorbed monomers (gC m-2)", "Polymers (gC m-2)", "Free enzymes (gC m-2)", "Adsorbed enzymes (gC m-2)")] <- bigdf[,c("Microbial structural biomass (gC m-2)", "Microbial reserve biomass (gC m-2)","Free monomers  (gC m-2)", "Adsorbed monomers (gC m-2)", "Polymers (gC m-2)", "Free enzymes (gC m-2)", "Adsorbed enzymes (gC m-2)")]*conv

bigdf$SOC <- bigdf$`Microbial structural biomass (gC m-2)` + bigdf$`Microbial reserve biomass (gC m-2)` + bigdf$`Free monomers  (gC m-2)` + bigdf$`Adsorbed monomers (gC m-2)` + bigdf$`Polymers (gC m-2)` + bigdf$`Free enzymes (gC m-2)` + bigdf$`Adsorbed enzymes (gC m-2)`
drops <- c("CUE","defactoTurnover")
writedf <- bigdf[ , !(names(bigdf) %in% drops)]

write.csv(writedf, file = "/Users/rzabramoff/Dropbox (Climate)/External Collaborations/Interface/code_scripts/LBL_model_Interface_output.csv")