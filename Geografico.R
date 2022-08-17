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

options(digits=2)
#Función para leer todas las hojas del excel
read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}


########arreglar base 
toxo<- read_excel_allsheets("TOXO.xlsx")
toxoDat<-toxo[[3]]
DatoTotalTox<-toxoDat[,c(1, 2, 3,4)]
DatoTotalTox<-na.omit(DatoTotalTox)
colnames(DatoTotalTox)[1]<-"Ano"
DatoTotalTox$total<-rowSums(DatoTotalTox[,c(3,4)])
DatoTotalTox <- DatoTotalTox %>%
  group_by(Departamento) %>%
  summarize(mean_toxo=mean(total))
write.csv(DatoTotalTox, file="DatoTotalTox.csv")

prep<-  read_excel_allsheets("pr_timeseries_annual_cru_1901-2020_COL 2.xlsx")
prepDat<-prep[[1]]
prepDat<-prepDat[c(2,117:121),]
colnames(prepDat)<-prepDat[1,]
prepDat<-prepDat[-1,]
prepDat<-prepDat[,-c(1,2)]
for(i in 1:33){
  prepDat[,i]<-as.numeric(prepDat[,i])
}
prepDat<-colMeans(prepDat)
write.csv(prepDat, file="prepDat.csv")
################3
########3########


library("dplyr")
load("coL_Dep.RData")

toxoPlot<-  read_excel_allsheets("prep_toxo_mapa.xlsx")
toxoPlot<-  toxoPlot[[1]]

#coL_Dep1  <- merge(x =col_Dep, y = toxoPlot)

coL_Dep1  <- col_Dep %>% left_join(toxoPlot)

library("ggspatial")

flcities <-coL_Dep1[coL_Dep1$NOMBRE_DPT %in% c("LA GUAJIRA", "ATLANTICO","NORTE DE SANTANDER", "SANTANDER","CAQUETA", "QUINDIO"), ]


pdf(file = "mean_prep.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 7)
ggplot(data = coL_Dep1) +
  geom_sf(fill = "antiquewhite1") +
  labs(fill = "Precipitation") +
  geom_sf(data = coL_Dep1, aes(fill = mean_prep/1000))+scale_fill_distiller(palette = "Oranges", direction = 1)+
  annotation_scale(location = "bl", width_hint = 0.4)+
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.1, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering)+  xlab("Longitude") + ylab("Latitude") +
  ggtitle("Colombia precipitation map", subtitle = "Observation Sites")+
  theme(panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", 
                                        size = 0.5), panel.background = element_rect(fill = "aliceblue"))+
  geom_text_repel(data = flcities, aes(coords_x, coords_y, label = NOMBRE_DPT), 
                  fontface = "bold",  arrow = arrow(length = unit(0.02, "npc")),
                  box.padding = 2, size=2)
dev.off()

####

  
library(spdep)
library(INLA)

toxoTtotal<- read_excel_allsheets("ToxoCARMODEL.xlsx")

ET<-toxoTtotal[[3]]
DR<-toxoTtotal[[1]]

quitar<-grep(pattern = 'Dep_nodef', x = ET$Departamento) 
ET<-ET[-quitar,]
quitar<-grep(pattern = 'TOXO SISPRO', x = ET$Departamento) 
ET<-ET[-quitar,]
ET<-ET[1:165,]
quitar<-grep(pattern = 'Dep_nodef', x = DR$Departamento) 
DR<-DR[-quitar,]
quitar<-grep(pattern = 'TOXO SISPRO', x = DR$Departamento) 
DR<-DR[-quitar,]

ranmu<-grep(pattern = 'M_[0-9]', x = colnames(ET))  #rango mujeres
ranho<-grep(pattern = 'H_[0-9]', x = colnames(ET))#rango hombres
rantotal<-c(ranmu,ranho) #total conteo
mapPreva1<-as.data.frame(cbind(ET$...1, ET[,rantotal]))

ranmu<-grep(pattern = 'M_[0-9]', x = colnames(DR))  #rango mujeres
ranho<-grep(pattern = 'H_[0-9]', x = colnames(DR))#rango hombres
rantotal<-c(ranmu,ranho) #total conteo

mapPreva2<-as.data.frame(cbind(DR$Año, DR[,rantotal]))


ranmu<-grep(pattern = 'M_[0-9]', x = colnames(ET))  #rango mujeres
ranho<-grep(pattern = 'H_[0-9]', x = colnames(ET))#rango hombres
rantotal<-c(ranmu,ranho) #total conteo

Y<-round(rowSums(ET[,rantotal]))

Oij<-mapPreva1
Pij<-mapPreva2
Pij[Pij==0]<-1

rij<-Oij[,-c(1)]/Pij[,-c(1)]
rij[is.na(rij)==TRUE]<-0
rj<-colSums(rij)

i=1
Rj<-Oij[,-c(1)]
for(i in 1:length(rj)){
  Rj[,i]<-rep(rj[i], dim(Rj)[1])	
}

Eij<-Pij[,-c(1)]*Rj

Ei<-rowSums(Eij)

E<-as.numeric(Ei)
E[E==0]<-0.0001


shp <- readOGR("depto.shp")
col_Dep <- st_as_sf(shp)
coords <- st_centroid(st_geometry(col_Dep), of_largest_polygon=TRUE) 
ind <- row.names(col_Dep)
col.tri.nb <- tri2nb(coords , row.names=ind)
shpnb.mat <- nb2mat(col.tri.nb, style="B",zero.policy=TRUE)
A<-shpnb.mat


nombre_Y<-as.data.frame(cbind(NOMBRE_DPT=ET$Departamento, Y))
nombre_Y$Y<-as.numeric(nombre_Y$Y)

nombre_Y <- nombre_Y %>%
  dplyr::group_by(NOMBRE_DPT) %>%
  dplyr::summarize(mean_Y=mean(Y))

nombre_E<-as.data.frame(cbind(NOMBRE_DPT=ET$Departamento, E))
nombre_E$E<-as.numeric(nombre_E$E)

nombre_E <- nombre_E %>%
  dplyr::group_by(NOMBRE_DPT) %>%
  dplyr::summarize(mean_E=mean(E))

shp<-merge(shp, nombre_E, by.x ="NOMBRE_DPT")
shp<-merge(shp, nombre_Y, by.x ="NOMBRE_DPT")
Xprep<-toxoPlot$mean_prep



data<-cbind(shp$mean_E, Xprep, round(shp$mean_Y))
data<-as.data.frame(data)
colnames(data)<-c("E", "X", "Y")



library(CARBayes)

formula <- Y  ~  X+offset(log(E))

model_spatial<-S.CARbym(formula=formula, data=data,family="poisson", W=A, burnin=10000, n.sample=16000)

model<-glm(Y  ~  X, data=data,family="poisson")

# Compute posterior SIR
y.fit <- model_spatial$samples$fitted
SIR <- t(t(y.fit) / data$E)

coL_Dep1$RR<-apply(SIR, 2, median)
coL_Dep1$E<-data$Y#/data$E

pdf(file = "RR.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 7)
ggplot(data = coL_Dep1) +
  geom_sf(fill = "antiquewhite1") +
  labs(fill = "RR") +
  geom_sf(data = coL_Dep1, aes(fill =RR))+scale_fill_distiller(palette = "Oranges", direction = 1)+
  annotation_scale(location = "bl", width_hint = 0.4)+
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.1, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering)+  xlab("Longitude") + ylab("Latitude") +
  ggtitle("Colombia relative risk map", subtitle = "Observation Sites")+
  theme(panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", 
                                        size = 0.5), panel.background = element_rect(fill = "aliceblue"))+
  geom_text_repel(data = flcities, aes(coords_x, coords_y, label = NOMBRE_DPT), 
                  fontface = "bold",  arrow = arrow(length = unit(0.02, "npc")),
                  box.padding = 2, size=2)
dev.off()










###probar en el compu





my_doc <- read_docx() %>%
  body_add_par(value = "Table of content", style = "heading 1") %>%
  body_add_toc(level = 2)

my_doc <- body_add_par(my_doc, "Mapas Colombia", style = "heading 1")

my_doc <- body_add_par(my_doc, "Precipitacones", style = "heading 2")

my_doc <-  body_add_gg(my_doc, value = gg, style = "centered")%>%
  body_add_par(value ="Precipitaciones", style = "Image Caption") %>%
  shortcuts$slip_in_plotref(depth = 2)


toxoPlot <-toxo %>% dplyr::group_by(departamentos) %>% dplyr::summarize(total=mean(Total))

Y<-toxoPlot$total

toxo<- read_excel_allsheets("TOXO.xlsx")
toxo<-toxo[[3]]
colnames(toxo)[1]<-"ano"

Oij<-toxo[1:33,c(9:26, 28:45, 47:64)]

Pij<-toxo[1:33,c(3:5)]

Oij<-na.omit(Oij)
Pij<-na.omit(Pij)

rij<-colSums(Oij/colSums(Pij))

i=1
Rj<-Pij

for(i in 1:length(rij)){
  Rj[,i]<-rep(rij[i], length(Pij[,1]))	
}

Eij<-t(Pij)*Rj

Ei<-rowSums(Eij)

#E[E==0]<-0.0001
E<-as.numeric(Ei)

tasa<-ifelse(E==0, 0, Y/E)

toxoTasa<-cbind(nombresDep, tasa)

toxoTasa<-as.data.frame(toxoTasa)
colnames(toxoTasa)<-c("NOMBRE_DPT", "Tasa")

coL_Dep1  <- coL_Dep1 %>% left_join(toxoTasa)

gg<-ggplot(data = coL_Dep1) + geom_sf(aes(fill = as.numeric(Tasa)))+ geom_text_repel(mapping = aes(coords_x, coords_y, label = NOMBRE_DPT), size = 2, min.segment.length = 0) +scale_fill_gradient(low = "white",high = "Red")+labs(fill = "Precipitación") 


my_doc <- body_add_par(my_doc, "Tasa incidencia Toxo", style = "heading 2")

my_doc <-  body_add_gg(my_doc, value = gg, style = "centered")%>%
  body_add_par(value ="Tasa incidencia Toxo", style = "Image Caption") %>%
  shortcuts$slip_in_plotref(depth = 2)


gg<-ggplot(data = coL_Dep1) + geom_sf(aes(fill = as.numeric(Tasa)))+scale_fill_gradient(low = "white",high = "Red")+labs(fill = "Tasa incidencia") 

my_doc <- body_add_par(my_doc, "Tasa incidencia", style = "heading 2")

my_doc <-  body_add_gg(my_doc, value = gg, style = "centered")%>%
  body_add_par(value ="Tasa incidencia Toxo", style = "Image Caption") %>%
  shortcuts$slip_in_plotref(depth = 2)


print(my_doc, target = "Geograficos.docx")




