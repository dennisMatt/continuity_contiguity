
#Script for reproducing analysis and figures in Dennis & Huck "The Continuity-Contiguity Problem in Fragmentation-Biodiversity Studies"

#load necessary libraries
library(terra)
library(matrixStats)
library(MatrixGenerics)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

#load the classified Random Forest Image
classRF<-rast("rfClasses.tif")



# load predictions (means) for each cell
allVals<-read.csv("rf_means.csv")

#get max prediction for each cell
v<-MatrixGenerics::rowRanges(as.matrix(allVals[,2:6]),useNames=FALSE)[,2]

#add to data frame
allVals$max=v

# add column for max prediction
allVals$class=NA

for(i in 1:nrow(allVals)){
  print(i)
  
  if(allVals$Woodland[i]==allVals$max[i]){
  allVals$class[i]=1}else{
    if(allVals$Grassland[i]==allVals$max[i]){
      allVals$class[i]=2}else{
        if(allVals$Shrub[i]==allVals$max[i]){
          allVals$class[i]=3}else{
            if(allVals$Water[i]==allVals$max[i]){
              allVals$class[i]=4}else{
                if(allVals$Urban[i]==allVals$max[i]){
                  allVals$class[i]=5}
}
}
      }
  }
}





prClass<-classRF


prClass[]<-allVals$class


########################Plotting habitat amount against alpha-cut

#isolate woodland class (Boolean) from RF results

prClass1<-prClass

prClass1[prClass1>1]=NA


#set up sequence of alpha-cuts

alpha<-seq(from=0.1,to=0.9,by=0.1)

#inititate lists to hold results
areaContin<-list() #continuous amount
areaContig<-list() #contiguous amount
areaContigA_Max<-list() #contig(a-max) amount
mpaContig<-list() # mean patch area of contigous amount
mpaContigA_Max<-list() # mean patch area of contigous(a-max) amount


woodMean<-classRF
woodMean[]<-allVals$Woodland

for(i in alpha){
  
  contin.i<-woodMean
  contin.i[contin.i<i]<-NA
  
  sContin<-sum(terra::values(contin.i),na.rm=T)*res(classRF)[2]^2
  
  
  contig.i<-contin.i
  
  contig.i[contig.i>0]<-1
  
  sContig<-sum(terra::values(contig.i),na.rm=T)*res(classRF)[2]^2
  
  
  contigAlphaMax<-prClass1*contin.i
  contigAlphaMax[contigAlphaMax>0]<-1
  
  
  sContigAlphaMax<-sum(terra::values(contigAlphaMax),na.rm=T)*res(classRF)[2]^2
  
  
  contigPoly<-patches(contig.i,directions=8,zeroAsNA=T)
  contigPoly<-as.polygons(contigPoly,dissolve=T)
  
  contigMPA<-mean(expanse(contigPoly))
  
  contigPolyMax<-patches(contigAlphaMax,directions=8,zeroAsNA=T)
  contigPolyMax<-as.polygons(contigPolyMax,dissolve=T)
  
  contigMPA_Max<-mean(expanse(contigPolyMax))
  
  
  areaContin[[i*10]]<-sContin
  areaContig[[i*10]]<-sContig
  areaContigA_Max[[i*10]]<-sContigAlphaMax
  mpaContig[[i*10]]<-contigMPA
  mpaContigMax[[i*10]]<-contigMPA_Max
  
  print(i)  

}



# get polygons for mean patch area calculation
contiguous_Max_patch<-patches(prClass1,directions=8,zeroAsNA=T)


contiguous_Max_patch<-as.polygons(contiguous_Max_patch,dissolve=T)



#build data frame for amount and MPA plots
dataArea<-data.frame(continuousArea=unlist(areaContin)/10000,
                     contiguousArea=unlist(areaContig)/10000,
                     contiguousArea_Amax=unlist(areaContigA_Max)/10000,
                     mpaContig=unlist(mpaContig)/10000,
                     mpaContig_AMax=unlist(mpaContigMax)/10000,
                     
                     Alpha=alpha,
                     Contiguous_Max_Amount=sum(values(prClass1)*res(prClass1)[2]^2,na.rm=TRUE)/10000,
                     Contiguous_Max_MPA=mean(expanse(contiguous_Max_patch),na.rm=T)/10000)



#set column nakes for plotting

colnames(dataArea)<-c("Continuous Amount","Contiguous(α) Amount","Contiguous (α-Max) Amount","Contiguous(α) MPA",
                      "Contiguous(α-Max) MPA","Alpha","Contiguous(Max) Amount", "Contiguous(Max) MPA")


#subset relevant columns
dataAreaPlot<-dataArea[,c(1,2,3,6,7)]

#convert format
dfArea <- melt(dataAreaPlot, id.vars='Alpha')



#ggplot
pArea<-ggplot(dfArea, aes(x=Alpha, y=value, fill=variable)) +
  geom_bar(stat='identity', position='dodge')+
  ylab("Area (ha)")+xlab("Alpha")+
  scale_x_continuous(breaks = seq(0.1, 0.9, by = 0.1))+
  #geom_bar(aes(Alpha,Boolean))+#,lwd=1,linetype=5)+
  theme_pubr( base_size = 30,
              base_family = "",
              border = FALSE,
              margin = TRUE,
              legend = c("top", "bottom", "left", "right", "none"),
              
              x.text.angle = 0)

pArea

#set visualisation parameters
ggpar(pArea,
      legend = "right", legend.title = " ",
      font.legend = c(30, "black"))+
      scale_fill_manual(values=c("dark green", 
                                 "orange",
                                 "blue",
                                 "red")) 
      


# for MPA plot
dataMPA<-dataArea[,c(4,5,6,8)]
#re-format
dfMPA <- melt(dataMPA, id.vars='Alpha')


#ggplot
pMPA<-ggplot(dfMPA,aes(x=Alpha,y=value,fill=variable))+
  geom_bar(stat='identity', position='dodge')+
  ylab("Area (ha)")+
  scale_x_continuous(breaks = seq(0.1, 0.9, by = 0.1))+

  
  theme_pubr( base_size = 30,
              base_family = "",
              border = FALSE,
              margin = TRUE,
              legend = c("top", "bottom", "left", "right", "none"),
              
              x.text.angle = 0)

pMPA

#visualisation params
ggpar(pMPA,
      legend = "right", legend.title = " ",
      font.legend = c(30, "black"))+
  scale_fill_manual(values=c("blue",
                             "dark green",
                              
                             "red")) 




#################################################

#########################################Plotting MinMax Values


#function to estimate parameters for Beta distribution

getAlphaBeta<-function(mu, sigma){
  alpha = mu**2 * ((1 - mu) / sigma**2 - 1 / mu)
  
  beta = alpha * (1 / mu - 1)
  
  return(c(alpha, beta))
}

# variance around prediction means for each cell in the landscape
woodVarVals<-read.csv("woodland_var.csv")

#raster holder
woodVar<-woodMean

#add values to raster 
woodVar[]<-wood_vals_var$variance

# compute beta distribution parameters for each cell
p<-getAlphaBeta(mu=as.numeric(matrix(woodMean)),sigma=as.numeric(matrix(woodlandVar)))

# convert to matrix
m<-matrix(p,ncol=2)

#get number of rows in matrix
nrows<-nrow(m)

#set length of Beta distribution
n=1000

#create a matrix containing beta distribution for each cell in the landscape 
rb <- matrix(rbeta(nrows * n, shape1 = m[,1], shape2 = m[,2]),
             ncol = nrows, byrow = TRUE)

#function to compute 1st and 99th percentile of values for each cell/location
funMinMax<-function(x){
  
  maxC=ecdf(rb[,x])
  q<-quantile(maxC,c(0.01,0.99))
  as.numeric(q)
  print(x)
  return(as.numeric(q))
}


# get min and max (1 and 99 percentile of values from beta distribution from each cell in the landscape)

n<-1:ncol(rb) # vector of cell numbers

#get min/max
minMaxCols<-lapply(n,FUN=funMinMax)

# collect to data frame
minMaxBeta<-do.call("rbind",minMaxCols)

# holder for min values
rMin<-woodMean

#add min values
rMin[]<-minMaxBeta[,1]

#holder for max values
rMax<-woodlandMean

#add max values
rMax[]<-minMaxBeta[,2]


###############################Contin/Contig min/max

#holders for continuous min/max >= alpha

#Min
minContin05<-rMin

#Max
maxContin05<-rMax

# set 0.5 alpha - min
minContin05[minContin05<0.5]=NA


#set 0.5 alpha - max
maxContin05[maxContin05<0.5]=NA


#plot

#build palette
cols <- brewer.pal(7, "Greens")
colsSD <- brewer.pal(7, "Blues")

#build colour ramp
pal <- colorRampPalette(cols)



plot(minContin05,col = pal(20),axes=F)
plot(maxContin05,col = pal(20),axes=F)


#for Boolean/contiguous calculation all cells >= alpha==1 

# Min
minContig05<-minContin05
minContig05[!is.na(minContig05)]<-1

#Max
maxContig05<-maxContin05
maxContig05[!is.na(maxContig05)]<-1

#get total area -min contin
minContinSum<-sum(values(minContin05),na.rm=T)*res(minContin)[2]^2

#Hectares
minContinSum<-minContinSum/10000

#get total area - max contin
maxContinSum<-sum(values(maxContin05),na.rm=T)*res(maxContin)[2]^2

maxContinSum<-maxContinSum/10000

##################################Contig

#get total area min Contig
minContigSum<-sum(values(minContig05),na.rm=T)*res(minContig)[2]^2

#hectares
minContigSum<-minContigSum/10000

#get total area max Contig
maxContigSum<-sum(values(maxContig05),na.rm=T)*res(maxContig)[2]^2

#hectares
maxContigSum<-maxContigSum/10000


################################### build data frame for plotting

dfMinMax<-data.frame(value=c(minContinSum,maxContinSum,minContigSum,maxContigSum),
                     minMax=c(0,1,0,1),Area=c("Continuous Amount","Continuous Amount","Contiguous(α) Amount",
                                              "Contiguous(α) Amount"))



dfMinMax$Area <- factor(dfMinMax$Area, levels = c("Continuous Amount", "Contiguous(α) Amount"))

head(dfMinMax)

pMinMax<-ggplot(dfMinMax, aes(x=minMax, y=value, fill=Area)) +
  geom_bar(stat='identity', position='dodge')+
  #geom_line(aes(Alpha,Boolean),lwd=1,linetype=5)+
  theme_pubr(base_size = 30,
              base_family = "",
              border = FALSE,
              margin = TRUE,
              #legend = c("top", "bottom", "left", "right", "none"),
              
              x.text.angle = 0)




pMinMax +ylab("Area (ha)")+xlab("level")
  

  
ggpar(pMinMax,legend = "top", legend.title = " ",
        font.legend = c(30, "black"))+
  xlab("Level")+
  ylab("Area (ha)")+
  scale_x_continuous(breaks = seq(0, 1, by = 1),labels=c("Minimum","Maximum"))+
  #scale_x_continuous("",breaks=c(2.5,6,9),)
  
  
  scale_fill_manual(values=c( 
                             "blue", 
                             "red")) 


