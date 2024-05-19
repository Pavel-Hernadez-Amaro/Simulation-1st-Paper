# Creation date: 29-07-2022

library(refund)
library(fda)
library(mgcv)
library(SOP)
library(expm)
library(lattice)
library(psych)
library(ggplot2)
library(readxl)
library(writexl)
library(Tplyr)
library(ggpubr)
library(tidyverse)
library(gridExtra)
# library(ggpubr)
library(remotes)
# install_github("Pavel-Hernadez-Amaro/VDPO")
library(VDPO)

options(error=NULL)

# time_data=read_csv(file.choose())
# 
# plot(table(time_data$fecha))
# 
# time_table=table(time_data$tel_p)
# 
# first_date=NULL
# first_date[1]=paste((time_data$fecha[1]))
# 
# for (i in 2:length(time_table)) {
#   
#    first_date[i]=paste((time_data$fecha[(i*time_table[i])+1]))
#   
# }

# load("C:/Users/user/Desktop/Trabajo/Escuela/Doctorado/Pavel/Tesis/Código/Mix models SOP/Propios/Datos_EPOC.RData")

load("C:/Users/user/Desktop/Trabajo/Escuela/Doctorado/Pavel/Tesis/Código/Mix models SOP/Propios/Datos_EPOC_SALVA.RData")

## HAY QUE CORRER LAS FUNCIONES EN FUNCIONES Finales.R

total_lcadl_baseline[109]=round(mean(total_lcadl_baseline[-109]))

o=order(rowSums(is.na(X)),decreasing = 1) # ESTO ME DA EL ORDER DE MENOS PASOS A MÁS

X=X[o,] # AQUÍ ESTOY ORDENANDO A LAS PERSONAS DE MENOS OBS A MÁS.

N=dim(X)[1]

T=dim(X)[2]

X_suave=matrix(nrow = N,ncol = T)

M_EPOC = dim(X)[2]-rowSums(is.na(X)) # LO MISMO QUE M_EPOC = ta[o]

datos$Sex=as.factor(datos$Sex)


c1=25
c2=25
c3=25

####################
#Sex+previous_adm_COPD+ansiedad_baseline+depresion_baseline

aa=model.matrix(prog_adm_COPD~Sex+previous_adm_COPD+ansiedad_baseline+depresion_baseline)
# aa=model.matrix(prog_adm_COPD~Edad+Sex+Sitl+AP1+AP2+AP3+AP4+AP5+PF2+PF4+bmi+bmicat+Indice_Charlson+indxcharl_cat+bode+previous_adm_COPD+ansiedad_baseline+depresion_baseline+total_stg_basal+w22_baseline+total_lcadl_baseline)
# aa=aa[,-1]

# ESTO ES PARA CREAR LA BASE INCLUYENDO LAS VARIABLES NO FUNCIONALES 

BB_EPOC=Data2B_simpson(X, M_EPOC, nbasis=c(c1,c2,c3),sub = 50)
E_EPOC=B2XZG(BB_EPOC$B,c=c(c2,c3)) # AQUÍ UTILIZO LA FUNCIÓN B2XZG() PARA GENERAR LAS MATRICES NECESARIAS EN UN MODELO MIXTO

X_all=cbind(aa[,-1],E_EPOC$X) # ESTE ES LA X DE LA PARTE FIJA QUE INCLUYE VARIABLES NO FUNCIONALES.

res_EPOC_all=XZG2theta(X = X_all, Z = E_EPOC$Z, G =  E_EPOC$G, T =  E_EPOC$T, y =  prog_adm_COPD, family = poisson(), offset = log((M_EPOC)/365)) # AQUÍ AJUSTO EL MODELO MIXTO Y RECUPERO LOS COEFICIENTES ORIGINALES
#

j=110 #ind_mortal[21]

prod=as.matrix(kronecker(BB_EPOC$B_Phi[[j]]$B,t(BB_EPOC$B_T$B[j,])))

Beta_h_all=prod %*% res_EPOC_all$theta # Aquí estamos multiplicando la base bidimesional por los theta

plot(Beta_h_all,type = "l",lwd=2,col=23, lty=3, ylim=c(-0.001,0.001),ylab = expression(paste(beta(t,T))),xlab = "Days")
# 

leyenda=NULL

iter_legend=1

for (j in c(12,24,45,66,110)) { #ind_mortal[seq(1,21,2)] 12,24,66,75
  
  # ind=which(unique(M_EPOC)== M_EPOC[j])
  # ind_t=max(M_EPOC) # OJO AQUÍ ESTAMOS SUPONIENDO QUE LAS OBSERVACIONES SON CONSECUTIVAS DESDE 1:MAX(M)
  # Beta_refund_all=Beta_G_all$value[(ind_t*(ind-1)+1):(ind*ind_t)]
  # Beta_refund_all=Beta_refund_all[1:M_EPOC[j]]
  
  ### SIN SUAVIZAR
  
  
  prod=as.matrix(kronecker(BB_EPOC$B_Phi[[j]]$B,t(BB_EPOC$B_T$B[j,])))
  
  # Beta_h=prod %*% res_EPOC$theta # Aquí estamos multiplicando la base bidimesional por los theta
  
  Beta_h_all=prod %*% res_EPOC_all$theta # Aquí estamos multiplicando la base bidimesional por los theta
  # Beta_h_mortalidad=prod %*% res_EPOC_mortalidad$theta 
  
  # Beta_h_inf = prod %*% lim_inf # Aquí estamos multiplicando la base bidimesional por los theta
  # Beta_h_sup = prod %*% lim_sup
  
  lines(Beta_h_all,type = "l",lwd=2,lty=iter_legend,col=j+1)
  
  # lines(Beta_h_mortalidad,type = "l",col=j+1,lwd=2)
  
  # lines(Beta_h_inf,lty=2,lwd=2,col=j+1)#,ylim=c(min(Beta[j,1:M[j]]),max(Beta[j,1:M[j]])))
  # lines(Beta_h_sup,lty=2,lwd=2,col=j+1)#,ylim=c(min(Beta[j,1:M[j]]),max(Beta[j,1:M[j]])))
  
  # lines(Beta_refund_all,col=j,type="l",lwd=2)
  
  leyenda[iter_legend]=paste(M_EPOC[j], "days")
  
  iter_legend=iter_legend+1
  
}


T_dat=IND=NULL
t_dat=rep(1:T,N)
Beta_estimated=array(dim=T*N)

for (ind in 1:N) {
  
  prod=as.matrix(kronecker(BB_EPOC$B_Phi[[ind]]$B,t(BB_EPOC$B_T$B[ind,])))
  Beta_estimated[(T*(ind-1)+1):((T*(ind-1))+M_EPOC[ind])]=prod %*% res_EPOC_all$theta # Aquí estamos multiplicando la base bidimesional por los theta
  T_dat=c(T_dat,rep(M_EPOC[ind],T))
  IND=c(IND,rep(ind,M_EPOC[ind]))
  
}


Heat_map_data=data.frame(t=t_dat,M=T_dat,Beta=Beta_estimated)
Heat_map_data=Heat_map_data[t_dat <= T_dat,]
Heat_map_data$IND=IND

lims <- range(Heat_map_data$Beta)

heat_map=if (requireNamespace("ggplot2", quietly = TRUE) &
             requireNamespace("RColorBrewer", quietly = TRUE)) {
  est <- Heat_map_data
  ggplot2::ggplot(est, ggplot2::aes(t, IND)) +
    ggplot2::geom_tile(ggplot2::aes(colour=Beta, fill=Beta)) +
    ggplot2::scale_fill_gradientn(  name="", limits=lims,
                                    colours=rev(RColorBrewer::brewer.pal(11,"Spectral"))) +
    ggplot2::scale_colour_gradientn(name="", limits=lims,
                                    colours=rev(RColorBrewer::brewer.pal(11,"Spectral"))) +
    ggplot2::scale_y_continuous(expand = c(0,0)) +
    ggplot2::scale_x_continuous(expand = c(0,0)) +
    geom_hline(yintercept = 12,linetype="dotted",color = "blue", linewidth=0.8, alpha=0.4) +
    geom_hline(yintercept = 61,linetype="dotted",color = "blue", linewidth=0.8, alpha=0.4) +
    geom_hline(yintercept = 102,linetype="dotted",color = "blue", linewidth=0.8, alpha=0.4) +
    ggplot2::theme_bw()+
    labs(y="")
}

heat_map

########################### HERE WE ARE GOING TO PERFORM A BOOSTRAP-LIKE INFERENCE FOR THE DATA

N=length(prog_adm_COPD)
B=200

theta_new=vector(mode = "list", length = B)
Beta_All=vector(mode = "list", length = B)

Beta_B=matrix(nrow = N,ncol = max(M_EPOC))

for (ind in 1:B) {

set.seed(1000+ind)
  
  print(paste0("ind = ",ind))
  
  y_new=rpois(N,res_EPOC_all$fit$fitted.values)
  
  res_EPOC_boostrap=XZG2theta(X = X_all, Z = E_EPOC$Z, G =  E_EPOC$G, T =  E_EPOC$T, y =  y_new, family = poisson(), offset = log((M_EPOC)/365)) # AQUÍ AJUSTO EL MODELO MIXTO Y RECUPERO LOS COEFICIENTES ORIGINALES
  
  theta_new[[ind]]=res_EPOC_boostrap$theta
  
  for (j in 1:N) {

    prod=as.matrix(kronecker(BB_EPOC$B_Phi[[j]]$B,t(BB_EPOC$B_T$B[j,])))
    
    Beta_B[j,1:M_EPOC[j]]=prod %*% res_EPOC_boostrap$theta # Aquí estamos multiplicando la base bidimesional por los theta
    
  }
  
  Beta_All[[ind]]=Beta_B
}

Beta_aux=Curva_inf=Curva_sup=vector(mode = "list",length = N)
aux=matrix(nrow = B,ncol=max(M_EPOC))

for (j in 1:N) {
for (i in 1:B) {
  aux[i,]=Beta_All[[i]][j,]
  
}
  Beta_aux[[j]]=aux
  }

for (ind in 1:N) {

  Curva_inf[[ind]]=apply(Beta_aux[[ind]], MARGIN = 2,FUN =  quantile,probs = c(0.05),na.rm=TRUE)[1:M_EPOC[ind]]
  Curva_sup[[ind]]=apply(Beta_aux[[ind]], MARGIN = 2,FUN =  quantile,probs = c(0.95),na.rm=TRUE)[1:M_EPOC[ind]]
  
}

j=77 #22 40 50 77

prod=as.matrix(kronecker(BB_EPOC$B_Phi[[j]]$B,t(BB_EPOC$B_T$B[j,])))

Beta_h_all=prod %*% res_EPOC_all$theta # Aquí estamos multiplicando la base bidimesional por los theta

plot(Beta_h_all,type = "l",lwd=2,col=23, lty=3, ylim=c(-0.001,0.001),ylab = expression(paste(beta(t,T))),xlab = "Days")
lines(Curva_inf[[j]])
lines(Curva_sup[[j]])
abline(h=0)

leyenda=NULL

iter_legend=1

for (j in c(12,24,45,66,110)) { #ind_mortal[seq(1,21,2)] 12,24,66,75

  prod=as.matrix(kronecker(BB_EPOC$B_Phi[[j]]$B,t(BB_EPOC$B_T$B[j,])))
  

  Beta_h_all=prod %*% res_EPOC_all$theta # Aquí estamos multiplicando la base bidimesional por los theta

  lines(Beta_h_all,type = "l",lwd=2,lty=iter_legend,col=j+1)

  lines(Curva_inf[[j]],type = "l",lwd=2,lty=iter_legend,col=j+1)
  lines(Curva_sup[[j]],type = "l",lwd=2,lty=iter_legend,col=j+1)
  
  leyenda[iter_legend]=paste(M_EPOC[j], "days")
  
  iter_legend=iter_legend+1
  
}


T_dat=IND=NULL
t_dat=rep(1:T,N)
Beta_estimated=array(dim=T*N)

for (ind in 1:N) {
  
  prod=as.matrix(kronecker(BB_EPOC$B_Phi[[ind]]$B,t(BB_EPOC$B_T$B[ind,])))
  Beta_estimated[(T*(ind-1)+1):((T*(ind-1))+M_EPOC[ind])]=prod %*% res_EPOC_boostrap$theta # Aquí estamos multiplicando la base bidimesional por los theta
  T_dat=c(T_dat,rep(M_EPOC[ind],T))
  IND=c(IND,rep(ind,M_EPOC[ind]))
  
}


Heat_map_data=data.frame(t=t_dat,M=T_dat,Beta=Beta_estimated)
Heat_map_data=Heat_map_data[t_dat <= T_dat,]
Heat_map_data$IND=IND

lims <- range(Heat_map_data$Beta)

heat_map_boostrap=if (requireNamespace("ggplot2", quietly = TRUE) &
             requireNamespace("RColorBrewer", quietly = TRUE)) {
  est <- Heat_map_data
  ggplot2::ggplot(est, ggplot2::aes(t, IND)) +
    ggplot2::geom_tile(ggplot2::aes(colour=Beta, fill=Beta)) +
    ggplot2::scale_fill_gradientn(  name="", limits=lims,
                                    colours=rev(RColorBrewer::brewer.pal(11,"Spectral"))) +
    ggplot2::scale_colour_gradientn(name="", limits=lims,
                                    colours=rev(RColorBrewer::brewer.pal(11,"Spectral"))) +
    ggplot2::scale_y_continuous(expand = c(0,0)) +
    ggplot2::scale_x_continuous(expand = c(0,0)) +
    geom_hline(yintercept = 12,linetype="dotted",color = "blue", linewidth=0.8, alpha=0.4) +
    geom_hline(yintercept = 61,linetype="dotted",color = "blue", linewidth=0.8, alpha=0.4) +
    geom_hline(yintercept = 102,linetype="dotted",color = "blue", linewidth=0.8, alpha=0.4) +
    ggplot2::theme_bw()+
    labs(y="")
}

heat_map_boostrap

