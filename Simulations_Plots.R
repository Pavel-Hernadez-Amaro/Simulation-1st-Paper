# load("C:/Users/user/Desktop/Trabajo/Escuela/Doctorado/Pavel/Tesis/1er Paper/Codigo/Simulation-1st-Paper/Simulations_N_100_X_noise_Poisson NegBin.RData")

###### mean
err_B_FFVD=apply(X = error_Beta_FFVD, MARGIN = 2, FUN = mean)
err_B_VD=apply(X = error_Beta_VD, MARGIN = 2, FUN = mean)

err_Y_FFVD=apply(X = Y_ERROR_2_sop, MARGIN = 2, FUN = mean)
err_Y_VD=apply(X = Y_ERROR_2_ellos, MARGIN = 2, FUN = mean)
err_Y_SOF=apply(X = Y_ERROR_2_te_ellos, MARGIN = 2, FUN = mean)

###### var
err_B_FFVD_var=apply(X = error_Beta_FFVD, MARGIN = 2, FUN = var)
err_B_VD_var=apply(X = error_Beta_VD, MARGIN = 2, FUN = var)

err_Y_FFVD_var=apply(X = Y_ERROR_2_sop, MARGIN = 2, FUN = var)
err_Y_VD_var=apply(X = Y_ERROR_2_ellos, MARGIN = 2, FUN = var)
err_Y_SOF_var=apply(X = Y_ERROR_2_te_ellos, MARGIN = 2, FUN = var)



B_ERRORES=data.frame(error_Beta_FFVD,error_Beta_VD) #B_ERROR_2,B_ERROR_2_te,

Y_ERRORES=data.frame(Y_ERROR_2_ellos,Y_ERROR_2_te_ellos,Y_ERROR_2_sop) # #Y_ERROR_2,Y_ERROR_2_te,

# write_xlsx(B_ERRORES,path = "C:/Users/user/Desktop/Trabajo/Escuela/Doctorado/Pavel/Tesis/Código/Mix models SOP/Propios/Simulaciones K-fold good version N_200/Poisson_Uniform_N_200/B_ERRORES.xlsx")
# 
# write_xlsx(Y_ERRORES,path = "C:/Users/user/Desktop/Trabajo/Escuela/Doctorado/Pavel/Tesis/Código/Mix models SOP/Propios/Simulaciones K-fold good version N_200/Poisson_Uniform_N_200/Y_ERRORES.xlsx")

Y_values=cbind(c(Y_ERRORES[,1],Y_ERRORES[,2],Y_ERRORES[,3]))
# Y_values_1=cbind(c(Y_ERRORES[,1],Y_ERRORES[,3],Y_ERRORES[,5]))
Y_values_2=cbind(c(Y_ERRORES[,2],Y_ERRORES[,10],Y_ERRORES[,18]))
Y_ERRORES_DF=data.frame(values=Y_values_2,method=as.factor(c(rep("VDFR",R),rep("SOF",R),rep("FF-VDFR",R))))

aux=Y_values[which(Y_values<50)] #Y_values[which(Y_values<quantile(Y_values,probs = 0.95))] #quantile(Y_values,probs = 0.95)
# aux=Y_values_2[which(Y_values_2<quantile(Y_values_2,probs = 0.95))] #quantile(Y_values,probs = 0.95)
Y_ERRORES_DF=data.frame(values=aux,method=as.factor(c(rep("VDFR",R),rep("SOF",R-35),rep("FF-VDFR",R-1)))) # POISSON UNIFORM R-7,R,R Y_1 # POISSON NegBin R-12,R-1,R-4 Y_1 # POISSON NegBin R-7,R-2,R-2 Y_2


 #quantile(Y_values_1,probs = 0.95)

Y_plot=ggplot(Y_ERRORES_DF,aes(x=method,y=values,fill=method))+
  geom_violin()+
  ylab("Error") +
  scale_x_discrete(labels = NULL, breaks = NULL)+
  xlab("")+
  theme_bw() +
  stat_summary(fun=median, geom="point", size=2, color="black")

Y_plot

B_values=cbind(c(B_ERRORES[,3],B_ERRORES[,11]))
# B_values_1=cbind(c(B_ERRORES[,1],B_ERRORES[,3]))
B_values_2=cbind(c(B_ERRORES[,2],B_ERRORES[,10]))
B_ERRORES_DF=data.frame(values=B_values_2,method=as.factor(c(rep("FF-VDFR",R),rep("VDFR",R))))

aux=B_values_2[which(B_values_2<quantile(B_values_2,probs = 0.95))] #B_values_2[-c(76,97)]  #B_values_2[which(B_values_2<quantile(B_values_2,probs = 0.95))] #quantile(B_values_2,probs = 0.95)
B_ERRORES_DF=data.frame(values=aux,method=as.factor(c(rep("FF-VDFR",R-6),rep("VDFR",R-4)))) # POISSON NegBin R-6,R-4 B_2 # POISSON NegBin round3 R-1,R-4 B_2

mean(aux[1:46]) #VDPO
mean(tail(aux,49)) #VD-Gellar

var(aux[1:46]) #VDPO
var(tail(aux,49)) #VD-Gellar

B_plot=ggplot(B_ERRORES_DF,aes(x=method,y=values,fill=method))+
  geom_violin()+
  ylab("Error") +
  scale_x_discrete(labels = NULL, breaks = NULL)+
  xlab("")+
  theme_bw() +
  stat_summary(fun=median, geom="point", size=2, color="black")

B_plot

plot(X_s[30,])
lines(X_se[30,])

case=96
case_iter=which.min(B_ERRORES_DF[1:R,1])

plot(Beta[case,,iter_out],ylim=c(-1,1))
lines(Beta_FFVD[case,,case_iter,iter_out],col="lightgreen",lwd=2,ylim=c(-2.5,2.5))
lines(Beta_FFVD_inf_ind[[case_iter]][case,],col="blue",lty=2,lwd=2)
lines(Beta_FFVD_sup_ind[[case_iter]][case,],col="red",lty=2,lwd=2)
abline(h=0)
