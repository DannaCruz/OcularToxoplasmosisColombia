
setwd("~/OneDrive - Universidad del rosario/UROSARIO/ProyectoIDEAM/Datos Procesados")
require(readxl)
require(tidyverse)

library(maps)
library(plyr)
library(fields)
library(leaflet)
library(RColorBrewer)
library(maptools)
library(spdep)
library(rgdal)
library(mvtnorm)
library(readr)
library(MASS)
library(GGally)
library(ggnet)
library(network)
library(sna)
library(ggplot2)
library(invgamma)
library(bnlearn)
library(BNSL)
library(xtable)
library(corrplot)
library(expm)
library(truncnorm)
library(matrixcalc)
library(pracma)
library(LaplacesDemon)
library(MCMCpack)
library(viridis)
library(ggrepel)
library(purrr)
library(dplyr)
library(ggplot2)
library(gganimate)
#theme_set(theme_bw())
library(ape)
library(readxl)
library(flextable)
library(officer)
library(magrittr)
library("Cairo")


library("dlnm")
library("splines")
library(ggplot2)

library(compareGroups)
library(qwraps2)
library(ggmosaic)





#set the working directory from which the files will be read from
setwd("/Users/dannacruz/OneDrive - Universidad del rosario/UROSARIO/ProyectoIDEAM/Datos Procesados")

#create a list of the files from your target directory
file_list <- list.files(path="/Users/dannacruz/OneDrive - Universidad del rosario/UROSARIO/ProyectoIDEAM/Datos Procesados/ProcessFiles")

#initiate a blank data frame, each iteration of the loop will append the data from the given file to this variable
dataset <- data.frame()

setwd("/Users/dannacruz/OneDrive - Universidad del rosario/UROSARIO/ProyectoIDEAM/Datos Procesados/ProcessFiles")

for (i in 1:length(file_list)){
  temp_data <- read.csv(file_list[i])
  temp  <- temp_data[,1]
  prep<-read.csv(file_list[i])[,3]
  dataPreTm<-as.data.frame(cbind(temp, prep))
  dataPreTm$clase<- rep(gsub("[^[:digit:]]", "", file_list[i]), dim(dataPreTm)[1])
  dataset <- rbind(dataset, dataPreTm) 
}


setwd("/Users/dannacruz/OneDrive - Universidad del rosario/UROSARIO/ProyectoIDEAM/Datos Procesados")

Codigos_por_departamento<-read.csv("codigosDep.csv", sep=";")

departamento<-NULL


for(i in 1:dim(dataset)[1]){
  if(length(which(Codigos_por_departamento$CODIGO==dataset$clase[i]))>0){
    departamento[i]<-Codigos_por_departamento$DEPARTAMENTO[
      which(Codigos_por_departamento$CODIGO==dataset$clase[i])]}else{
        departamento[i]<-"Choco"
      }
}




options(digits=2)
#Función para leer todas las hojas del excel
read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}


#carga de datos

toxo<- read_excel_allsheets("TOXO.xlsx")

toxoDat<-toxo[[3]]


DatoTotalTox<-toxoDat[,c(1, 2, 3,4)]

DatoTotalTox<-na.omit(DatoTotalTox)


colnames(DatoTotalTox)[1]<-"Ano"



Dep<-unique(dataset$departamento)

dataset$departamento<-departamento

##############


##############



my_doc <- read_docx() %>%
  body_add_par(value = "Table of content", style = "heading 1") %>%
  body_add_toc(level = 2)


for(de in 1:length(Dep))
{
  
d<-Dep[de]

my_doc <- body_add_par(my_doc, d, style = "heading 1") 

prep<- dataset[which(dataset$departamento==d),]

tox<-DatoTotalTox[which(DatoTotalTox$Departamento==d),c(1,5)]

prep$ano<- substr(prep$temp,1, 4)

for(i in 1:dim(prep)[1])
{
  prep$tox[i]<-tox$Total[which(tox[,1]==prep$ano[i])]
}
  

cb1.pm <- crossbasis(as.numeric(prep$prep), lag=15, argvar=list(fun="lin"),
                     arglag=list(fun="poly",degree=7))

prep$tiempo<-1:length(prep$temp)

prep$dia <- weekdays(as.Date(prep$temp))

model1 <- glm(tox ~ cb1.pm + ns(tiempo, 7*14) + dia, family=quasipoisson(), prep)

pred1.pm <- crosspred(cb1.pm, model1, at=0:20, bylag=0.5, cumul=TRUE)

my_doc <- body_add_plot(my_doc, value = plot_instr(
  code = { plot(pred1.pm, "slices", var=20, col=3, type="p", ylab="RR", ci="bars", main="Curva de respuesta de retraso precipitación")
  }), style = "centered" )

#my_doc <- body_add_plot(my_doc, value = plot_instr(
 # code = { hist(cb1.pm[,5],  col = 2, xlab=" Precipitation (ºc)")}), style = "centered" )

my_doc <- body_add_plot(my_doc, value = plot_instr(
  code = { plot(pred1.pm, "slices",cumul=TRUE, var=20, col=3, ylab="Accumulated RR",  main="Curva de respuesta Acumulada de retraso  precipitación")
  }), style = "centered" )


my_doc <- body_add_par(my_doc, "Porcentajes de cambio (PC)", style = "heading 2") 
# Porcentajes de cambio (PC)

betas<-coef(model1)
Ses<-coef(summary(model1))[, "Std. Error"]

li<-function(betaS, SE){
  li<-(exp((betaS+1.96*SE))-1)
  return(round(li,3))
}

ls<-function(betas, SE){
  ls<-(exp((betas-(1.96*SE)))-1)
  return(round(ls,2))
}

m<-function(betas, SE){
  m<-ls(betas, SE)-SE
  return(round(m,2))
}


PM10_0<-cbind(m(betas[2], Ses[2]), li(betas[2], Ses[2]), ls(betas[2], Ses[2]))
PM10_1<-cbind(m(betas[3], Ses[3]), li(betas[3], Ses[3]), ls(betas[3], Ses[3]))
PM10_3<-cbind(m(betas[4], Ses[4]), li(betas[4], Ses[4]), ls(betas[4], Ses[4]))
PM10_4<-cbind(m(betas[5], Ses[5]), li(betas[5], Ses[5]), ls(betas[5], Ses[5]))

PM10<-rbind(PM10_0, PM10_1, PM10_3, PM10_4)
colnames(PM10)<-c("media", "LI", "LS")
rownames(PM10)<-c("lag0", "lag1", "lag2", "lag3")
PM10<-cbind(rownames(PM10), PM10)
colnames(PM10)[1]<-"lag"
PM10<-as.data.frame(PM10)

#gg<-ggplot(PM10, aes(x  = lag, y  = media))  +
 # geom_point()+
#  geom_errorbar(aes(ymin=LS, ymax=LI), width=.2,
         #       position=position_dodge(.9)) 


# <-  body_add_gg(my_doc, value = gg, style = "centered" )%>%
#  body_add_par(value = d, style = "Image Caption") %>%
##  shortcuts$slip_in_plotref(depth = 2)%>% body_end_section_portrait()


#my_doc <-  body_add_table(my_doc, value = PM10, style = "table_template")%>%
 # body_add_par(value = d, style = "Table Caption")%>%
#  shortcuts$slip_in_tableref(depth = 2)
}

print(my_doc, target = "DESCRIPTIVASFINAL.docx")



########





for(de in 1:length(Dep))
{
  d<-Dep[de]
  
  my_doc <- body_add_par(my_doc, d, style = "heading 1") 
  
  prep<- dataset[which(dataset$departamento==d),]
  
  tox<-DatoTotalTox[which(DatoTotalTox$Departamento==d),c(1,5)]
  
  prep$ano<- substr(prep$temp,1, 4)
  
  for(i in 1:dim(prep)[1])
  {
    prep$tox[i]<-tox$Total[which(tox[,1]==prep$ano[i])]
  }
  
  
  cb1.pm <- crossbasis(as.numeric(prep$prep), lag=15, argvar=list(fun="lin"),
                       arglag=list(fun="poly",degree=7))
  
  prep$tiempo<-1:length(prep$temp)
  
  prep$dia <- weekdays(as.Date(prep$temp))
  

}

print(my_doc, target = "DESCRIPTIVASFINAL.docx")

##### TABLAS RIESGO RELATIVOOO


tablaRR<-NULL

nombrestablaRR<-NULL
for(de in 1:length(Dep))
{
  d<-Dep[de]
  nombrestablaRR[de]<-d
  
  prep<- dataset[which(dataset$departamento==d),]
  tox<-DatoTotalTox[which(DatoTotalTox$Departamento==d),c(1,5)]
  prep$ano<- substr(prep$temp,1, 4)
  for(i in 1:dim(prep)[1])
  {
    prep$tox[i]<-tox$Total[which(tox[,1]==prep$ano[i])]
  }
  cb1.pm <- crossbasis(as.numeric(prep$prep), lag=15, argvar=list(fun="lin"),
                       arglag=list(fun="poly",degree=7))
  
  prep$tiempo<-1:length(prep$temp)
  prep$dia <- weekdays(as.Date(prep$temp))
  model1 <- glm(tox ~ cb1.pm + ns(tiempo, 7*14) + dia, family=quasipoisson(), prep)
  pred1.pm <- crosspred(cb1.pm, model1, at=0:20, bylag=0.5, cumul=TRUE)
  options(digits = 5)
  
  tabla<- pred1.pm[c("matRRfit","matRRlow","matRRhigh")]
  
  intervalosRR<-lapply(tabla, colMeans)
  
  tablaRlag0<-paste(paste(signif(intervalosRR$matRRfit["lag0"], 4), paste(signif(intervalosRR$matRRlow["lag0"],3), signif(intervalosRR$matRRhigh["lag0"],3), sep= ","), sep = " ("), ")")
  tablaRlag6<-paste(paste(signif(intervalosRR$matRRfit["lag6"], 4), paste(signif(intervalosRR$matRRlow["lag6"],3), signif(intervalosRR$matRRhigh["lag6"],3), sep= ","), sep = " ("), ")")
  tablaRlag10<-paste(paste(signif(intervalosRR$matRRfit["lag10"], 4), paste(signif(intervalosRR$matRRlow["lag10"],3), signif(intervalosRR$matRRhigh["lag10"],3), sep= ","), sep = " ("), ")")
  tablaRlag13<-paste(paste(signif(intervalosRR$matRRfit["lag13"], 4), paste(signif(intervalosRR$matRRlow["lag13"],3), signif(intervalosRR$matRRhigh["lag13"],3), sep= ","), sep = " ("), ")")
  tablaRlag15<-paste(paste(signif(intervalosRR$matRRfit["lag15"], 4), paste(signif(intervalosRR$matRRlow["lag15"],3), signif(intervalosRR$matRRhigh["lag15"],3), sep= ","), sep = " ("), ")")
  
  tablaRResul<-cbind(tablaRlag0, tablaRlag6, tablaRlag10, tablaRlag13,tablaRlag15)
  
  tablaRR<-rbind(tablaRR,tablaRResul)
}


my_doc <- read_docx() %>%
  body_add_par(value = "Table of content", style = "heading 1") %>%
  body_add_toc(level = 2)

tablaRR<-as.data.frame(tablaRR)
rownames(tablaRR)<-nombrestablaRR


tablaRR<-cbind(rownames(tablaRR), tablaRR)

my_doc <-  body_add_table(my_doc, value =  tablaRR, style = "table_template")%>%
  body_add_par(value = d, style = "Table Caption")%>%
  shortcuts$slip_in_tableref(depth = 2)

print(my_doc, target = "tablaRRL.docx")


###################################
#############################

tablaprep<-NULL
for(de in 1:length(Dep))
{
  d<-Dep[de]
  nombrestablaRR[de]<-d
  
  prep<- dataset[which(dataset$departamento==d),]
  tox<-DatoTotalTox[which(DatoTotalTox$Departamento==d),c(1,5)]
  prep$ano<- substr(prep$temp,1, 4)
  i=1
  for(i in 1:dim(prep)[1])
  {
    prep$tox[i]<-tox$Total[which(tox[,1]==prep$ano[i])]
  }
  tablaprep<-rbind(tablaprep,prep)
}

class(tablaprep$departamento)
tablaprep$prep<-as.numeric(tablaprep$prep)
tablaprepgRUPO <-  tablaprep %>% dplyr::group_by(departamento) %>%  dplyr::summarize(mean_prep=mean(prep), mean_tox=mean(tox))

cor(tablaprepgRUPO$mean_prep, tablaprepgRUPO$mean_tox)


f<-"departamento ~ mean_prep"

res <- compareGroups(f, data=tablaprepgRUPO, max.ylev = 40, max.xlev = 40)

restab<- createTable(res)

#####
f<-"departamento ~ prep + "

res <- compareGroups(f, data=tablaprep, max.ylev = 40, max.xlev = 40)

restab<- createTable(res)

export2xls(restab, file="resumenprep.xlsx")
#####

export2xls(restab, file="resumenDepTox.xlsx")

f<-"ano ~tox + prep"

res <- compareGroups(f, data=tablaprep)

restab<- createTable(res)

export2xls(restab, file="resumenderestabTox.xlsx")

#####

f<-"Total ~ Femenino +Masculino + Departamento"

res <- createTable(compareGroups(f , data=DatoTotalTox, include.miss = TRUE,max.ylev = 40, max.xlev = 40), show.n = TRUE)


export2xls(res, file="resumendescripti.xlsx")




#####
####
####



my_doc <- read_docx() %>%
  body_add_par(value = "Table of content", style = "heading 1") %>%
  body_add_toc(level = 2)


for(de in 1:length(Dep))
{
  
  d<-Dep[de]
  
  my_doc <- body_add_par(my_doc, d, style = "heading 1") 
  
  prep<- dataset[which(dataset$departamento==d),]
  
  tox<-DatoTotalTox[which(DatoTotalTox$Departamento==d),c(1,5)]
  
  prep$ano<- substr(prep$temp,1, 4)
  
  for(i in 1:dim(prep)[1])
  {
    prep$tox[i]<-tox$Total[which(tox[,1]==prep$ano[i])]
  }
  
  
  prep$prep<- as.numeric(prep$prep)
  prep$tox<-as.numeric(prep$tox)

  correla<-cor.test(prep$prep, prep$tox)
  
  prepG<-prep %>%
    dplyr::group_by(ano) %>%
    dplyr::summarize(mean_prep=mean(prep), mean_tox=mean(tox))
  
  correla<-cor.test(as.numeric(prepG$mean_prep), as.numeric(prepG$mean_tox))
  tpvalue<-t.test(as.numeric(prepG$mean_prep), as.numeric(prepG$mean_tox))
  tablaFi<-data.frame(c(round(correla$estimate,2), format.pval(tpvalue$p.value), round(mean(as.numeric(prep$prep)),2)))
  tablaFi<-cbind(c("Corr", "p-value", "prep"), tablaFi)
  colnames(tablaFi)[1]<-"Result"
  TablaF<-rbind(TablaF, t(tablaFi)[2,])
  rownames(TablaF)[dim(TablaF)[1]]<-d
  
}

  
print(my_doc, target = "correla2.docx")



#####
####respuesta de manuscrito
####

deparSelec<-c("Guajira",  
"Atlántico", 
"Norte de Santander",  
"Santander",  
"Caquetá",
"Quindío")



for(de in 2:length(deparSelec))
{
  
  d<-deparSelec[de]
  prep<- dataset[which(dataset$departamento==d),]
  tox<-DatoTotalTox[which(DatoTotalTox$Departamento==d),c(1,5)]
  
  prep$ano<- substr(prep$temp,1, 4)
  
  for(i in 1:dim(prep)[1])
  {
    prep$tox[i]<-tox$Total[which(tox[,1]==prep$ano[i])]
  }
  prep$prep<- as.numeric(prep$prep)
  prep$tox<-as.numeric(prep$tox)
  prepG<-prep %>%
    dplyr::group_by(ano) %>%
    dplyr::summarize(mean_prep=mean(prep), mean_tox=mean(tox))
  
  tablaFi<-cbind(d, prepG[5,])
  TablaF<-rbind(TablaF, tablaFi)
}


deparSelec<-c("Guajira",  
              "Atlántico", 
              "Norte de Santander",  
              "Santander",  
              "Caquetá",
              "Quindío")

dataCorr<-cbind(deparSelec,c(1.3, 2.3, 4.2, 5.5,3.99, 2), c(11,1, 18,40,22.2, 13))



#dataCorr<-cbind(deparSelec,c(1.61, 2.67, 4.45, 5.39,3.99, 2), c(4,6.4, 8.2,28.2,22.2, 13))

dataCorr<-as.data.frame(dataCorr)

cor.test(c(1.3, 2.3, 4.2, 5.5,3.99, 2), c(11,1, 18,40,22.2, 13))

dataCorr$V2<-as.numeric(dataCorr$V2)
dataCorr$V3<-as.numeric(dataCorr$V3)
library(ggplot2)
#v2 precipitacion y v3 datos


dataCorr<-cbind(deparSelec[1:4],c(1.61, 2.67, 4.45, 5.39), c(4,6.4, 8.2,28.2))
dataCorr<-as.data.frame(dataCorr)

dataCorr<-dataCorr[order(dataCorr$V2), ]
dataCorr$V1<-factor(dataCorr$V1, ordered = TRUE)
levels(dataCorr$V1)<-dataCorr$V1

dataCorr$V2<-as.numeric(dataCorr$V2)
dataCorr$V3<-as.numeric(dataCorr$V3)

ggplot(dataCorr, aes(x =V1)) +
  geom_col(aes( y = V2, fill="redfill")) +
  geom_text(aes(y = V2, label = V3), fontface = "bold", vjust = 1.4, color = "black", size = 4) +
  geom_line(aes(y = V3, group = 1, color = 'blueline'))+
  geom_text(aes(y = V3, label = round(V3, 2)), vjust = 1.4, color = "blue", size = 3) +
  scale_fill_manual('', labels = 'Occurance', values = "#C00000") +
  scale_color_manual('', labels = 'Time Reshop', values = 'black') +
  theme_minimal()



