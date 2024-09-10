library(refund)
library(fda)
library(mgcv)
library(SOP)
library(expm)
library(writexl)
library(Tplyr)
library(ggpubr)
library(tidyverse)
library(gridExtra)
library(dismo)
options(error=NULL)
library(remotes)
# install_github("Pavel-Hernadez-Amaro/VDPO") 
library(VDPO)

########### Here we generate the data and set some fixed parameters

Rsq=0.95

N=100 # NUMBER OF SUBJECTS

J=100 # NUMBER OF MAXIMUM OBSERVATIONS PER SUBJECTS

k=4  # NUMBER OF GROUPS IN THE K-fold

fold=N/k

R=100 # NUMBER OF ITERATIONS FOR THE SIMULATION STUDY

c1=25
c2=50
c3=30

case=8 # THIS IS THE TOTAL NUMBER OF SCENARIOS SIMULATED

registrated_points=function(x) {
  J.i <- sum(!is.na(x))
  y <- approx(x=seq(0,1,length=J.i), y=x[1:J.i],
              xout=seq(0,1,length=J))$y
  y
}

time_no_vc=time_SOP=time_Gellar=time_Goldsmith=array(dim = c(R,case))

# # HERE WE GENERATE TO EMPTY ARRAYS THAT WILL BE THE ESTIMATED FUNCTIONAL COEFFICIENTS
Beta_VD=Beta_FFVD=array(dim =c(N,J,R,case))


error_Beta_FFVD=error_Beta_VD=array(dim =c(R,case))
# #

Beta_FFVD_sup_ind=Beta_FFVD_inf_ind=vector(mode = "list",length = R)

#ndb=15 # THIS IS THE NUMBER OF LAMBDAS FOR THE ADPATIVE CASE

# # HERE WE GENERATE TO EMPTY ARRAYS THAT WILL BE THE TRUE FUNCTIONAL COEFFICIENTS
Beta_h=Beta_h_no_vc=array(dim =c((k-1)*fold,J,R))
Beta_h_test=array(dim =c(fold,J,R))
# #

# EMPTY LISTS THAT WILL BE THE ESTIMATED FUNCTIONAL COEFFICIENTS USING THE Gellar APPROACH
Beta_G=list()
Beta_G_te=list()
#

# EMPTY ARRAYS WHERE THE DOMAINS OF EVERY SUBJECTS IN EVERY ITERATION WILL BE STORE
M_it=array(dim=c(R,N,case))
M_it_test=array(dim=c(R,fold,case))
M_it_train=array(dim=c(R,(k-1)*fold,case))
#

# EMPTY ARRAYS WHERE THE ESTIMATED RESPOND VARIABLE WILL BE STORE FOR EVERY ITERATION IN EVERY ESCANARIO
y_h_ellos=array(dim=c(R,fold,case))
y_h_ellos_te=array(dim=c(R,fold,case))
y_h_sop=array(dim=c(R,fold,case))
y_h_no_vc=array(dim=c(R,fold,case))
y_h_adaptive=array(dim=c(R,fold,case))
#

# EMPTY ARRAYS WHERE THE ERROR OF THE ESTIMATED COEFFICIENT FUNCTION WILL BE STORE FOR EVERY ITERATION IN EVERY ESCANARIO
B_ERROR_2_ellos=B_ERROR_2_te_ellos=array(dim = c(R,case))
B_ERROR_2_sop=B_ERROR_2_no_vc=B_ERROR_2_sop_ad=array(dim = c(R,case))
B_ERROR_2_ellos_test=B_ERROR_2_te_ellos_test=array(dim = c(R,case))
B_ERROR_2_sop_test=B_ERROR_2_sop_ad_test=array(dim = c(R,case))
#


# THIS ARE THE N/K ERRORS IN EACH GROUPS
error_group_FF_VDFR=error_group_SB=error_group_Carmen=error_group_VDFR=array(dim=c(R,k,case))

# EMPTY ARRAYS WHERE THE ERROR OF THE ESTIMATED RESPOND VARIABLE WILL BE STORE FOR EVERY ITERATION IN EVERY ESCANARIO
Y_ERROR_2_ellos=Y_ERROR_2_te_ellos=array(dim = c(R,case))
Y_ERROR_2_sop=Y_ERROR_2_no_vc=Y_ERROR_2_sop_ad=Y_ERROR_2_sop_manual=array(dim = c(R,case))
#

start=proc.time()

for (iter_out in 1) { # HERE WE CAN SELECTED WHICH SCENARIO(S) SIMULATE
  
  print(c("case = ",iter_out))
  
  for (iter in 1:R) {
    
    print(c("iter = ",iter))
    
    set.seed(1000+iter)
    
    M = round(runif(N,10,J),digits = 0) # HERE WE GENERATE THE DOMAIN FOR ALL SUBJECTS WITH A MINIMUM OF A 10 OBSERVATIONS
    
    # M = rnegbin(N,24,1) # Para 1000 poner (240,2)
    
    if (max(M)>J) {
      # print("si")
      M[which(M>J)]=J
    }
    
    
    if (min(M)<=10) {
      # print("si")
      # M[which(M<=10)]=round(runif(length(which(M<=10)),31,max(M)))
      M[which(M<=10)]=10
    }
    
    M=sort(M) # WE SORT THE DATA WITHOUT LOSS OF GENERALITY
    
    M_max=max(M)
    
    t=seq(from = 0, to = 1, length.out=M_max)
    
    # t=1:M_max
    
    M_it[iter,,iter_out]=M # WE STORE THE DOMAINS FOR EVERY ITERATION AND SCENARIO
    
    
    ############ HERE WE GENERATE THE FUNCTIONAL DATA
    
    X_se=matrix(NA,N,M_max) # NOISY
    
    X_s=matrix(NA,N,M_max) # NOT NOISY
    
    for (i in 1:N) {
      
      u=rnorm(1)
      
      temp=matrix(NA,10,M_max)
      
      for (e in 1:10) {
        
        v_i1=rnorm(1,0,4/e^2)
        v_i2=rnorm(1,0,4/e^2)
        
        temp[e,1:M[i]]=v_i1*sin(2*pi*e*(1:M[i])/100)+v_i2*cos(2*pi*e*(1:M[i])/100)
      }
      
      B=apply(temp,2,sum)
      
      B=B+u
      
      X_s[i,]=B
      X_se[i,]=(B)+rnorm(M_max,0,1.5) # WE ADD NOISE
      
    }
    
    X_reg = t(apply(X_se, 1,registrated_points)) # THESE ARE THE REGISTRATED POINTS
    
    # SOME SAVE CHECKS FOR UNWANTED NAs
    for (i in 1:M_max) {
      if (length(which(is.na(X_s[,i])))==N) {
        print(c("iteracion",l,"columna",i))
        stop("Hay columnas con NA")
      }}
    
    
    
    
    ###### HERE WE GENERATE THE TRUE Beta COEFFICIENT AND THE RESPONSE VARIABLE
    
    Beta=array(dim = c(N,M_max,4))
    nu=y=rep(0,N)
    
    for (i in 1:N) {
      
      # TRUE FUNCTIONAL COEFFICIENTS
      
      Beta[i,1:(M[i]),1]=((10*t[1:(M[i])]/M[i])-5)/10
      Beta[i,1:(M[i]),2]=((1-(2*M[i]/M_max))*(5-40*((t[1:(M[i])]/M[i])-0.5)^2))/10
      Beta[i,1:(M[i]),3]=(5-10*((M[i]-t[1:(M[i])])/M_max))/10
      Beta[i,1:(M[i]),4]=(sin(2*pi*M[i]/M_max)*(5-10*((M[i]-t[1:(M[i])])/M_max)))/10
      #
      
      
      # HERE WE GENERATE THE RESPONSE VARIABLE FOR EVERY BETA FROM THE 2 DIFFERENT FUNCTIONAL DATA (NOISY AND OTHERWISE)
      
      if (iter_out==8) {
        nu[i]=sum(X_se[i,]*Beta[i,,4],na.rm = 1)/(M[i]) # NOISY
      }
      if (iter_out<=4) {
        nu[i]=sum(X_s[i,]*Beta[i,,iter_out],na.rm = 1)/(M[i]) #NOT NOISY
      }
      if(iter_out>4 & iter_out<8){
        nu[i]=sum(X_se[i,]*Beta[i,,iter_out%%4],na.rm = 1)/M[i] # NOISY
      }
      
    }
    
    
    var_e <- (1 / Rsq - 1) * stats::var(nu) # (1-Rsq)*var(nu[ind,])
    
    y=nu+rnorm(N,sd = var_e) # ADDING NOISE TO THE GAUSSIAN MODEL
    
    # y=rpois(N,exp(nu)) # POISSON MODEL ##### CAHNGE THE ERROR EQUATIONS IN LINE 477
    
    # y=rbinom(N,1,(exp(nu)/(1+exp(nu)))) # BINOMIAL MODEL
    
    ############## HERE ENDS THE GENERATION OF THE DATA
    
    ############ HERE BEGINS THE ESTIMATION OF THE MODEL
    data=data.frame(y=y)
    data[["X_se"]]=X_se
    data[["X_s"]]=X_s
    data[["Beta"]]=Beta
    
    ############## HERE ENDS THE GENERATION OF THE DATA  
    
    ############ HERE BEGINS THE ESTIMATION OF THE MODEL
    
    # ASSIGNING THE TRAIN AND TEST SETS 
    
    groups=kfold(data,k=k) # creating the K-fold groups    
    
    for (group in 1:k) {
      
      print(c("group = ",group))
      
      X_test=X_se[which(groups==group),]
      X_train=X_se[-which(groups==group),]
      
      X_reg_train=X_reg[-which(groups==group),]
      X_reg_test=X_reg[which(groups==group),]
      
      M_test=M[which(groups==group)]
      M_train=M[-which(groups==group)]
      
      M_it_test[iter,,iter_out]=M_test
      M_it_train[iter,,iter_out]=M_train
      
      y_test=y[which(groups==group)]
      y_train=y[-which(groups==group)]
      
      # MODEL ESTIMATION BY OUR APPROACH

      data_train=data.frame(y_train)
      data_train[["X_train"]]=X_train
      
      data_test=data.frame(y_test)
      data_test[["X_test"]]=X_test
      
      start_SOP=proc.time()
      formula <- y_train ~ ffvd(X_train, t, nbasis = c(c1,c2,c3))
      BB=ffvd(X_train, t, nbasis = c(c1,c2,c3))
      res <- VDPO(formula = formula, data = data_train)#,family = poisson())
      
      
      end_SOP=proc.time()
      time_SOP[iter,iter_out]=end_SOP[1]-start_SOP[1]
      #
      
      # MODEL ESTIMATION USING GELLAR APPROACH
      
      start_Gellar=proc.time()
      
      fit_Gellar <- pfr(y_train ~ lf.vd(X_train,k=89))#,family = poisson())
      
      end_Gellar=proc.time()
      time_Gellar[iter,iter_out]=end_Gellar[1]-start_Gellar[1]
      #
      
      # MODEL ESTIMATION USING GOLDSMITH APPROACH
      
      start_Goldsmith=proc.time()
      
      fit_Gellar_te <- pfr(y_train ~ lf(X_reg_train, bs = "ps",k=25,presmooth = "bspline",presmooth.opts = list(nbasis=25)))#,family = poisson())
      
      end_Goldsmith=proc.time()
      time_Goldsmith[iter,iter_out]=end_Goldsmith[1]-start_Goldsmith[1]
      
      #################### HERE ENDS THE ESTIMATION OF THE MODELS
      
      #################### HERE WE CALCULATE THE ESTIMATION ERRORS
      
      error_2_ellos=error_2_te_ellos=matrix(0,nrow = (k-1)*fold,ncol = 1)
      error_2_sop=matrix(0,nrow = (k-1)*fold,ncol = 1)
      
      error_2_ellos_test=error_2_te_ellos_test=matrix(0,nrow = fold,ncol = 1)
      error_2_sop_test=matrix(0,nrow = fold,ncol = 1)
      
      ERROR_2_sop=0
      ERROR_2_ellos=0
      
      ############# HERE WE CALCULATE THE ERROR OF THE ESTIMATED RESPONSE VARIABLE
      
      BB_test=ffvd(X_test, t, nbasis = c(c1,c2,c3))
      
      y_h_sop[iter,,iter_out] = BB_test$B_ffvd %*% res$theta_ffvd # ESTIMATED REPSONSE VARIABLE USING OUR APPROACH
      
      # ESTIMATED REPSONSE VARIABLE USING HE GELLAR AND GOLDSMITH APPROACHES
      
      Beta_G[[iter]]=coef(fit_Gellar, n=c(length(t),length(unique(M))))$value
      Beta_G_te[[iter]]=coef(fit_Gellar_te, n=M_max)$value
      
      for (j in 1:fold) {
        
        ind=which(unique(M_it[iter,,iter_out])== M_it_test[iter,j,iter_out])
        ind_t=max(M_it[iter,,iter_out]) # OJO AQU? ESTAMOS SUPONIENDO QUE LAS OBSERVACIONES SON CONSECUTIVAS DESDE 1:MAX(M)
        
        Beta_refund=Beta_G[[iter]][(ind_t*(ind-1)+1):(ind*ind_t)]
        Beta_refund=Beta_refund[1:M_it_test[iter,j,iter_out]]
        
        Beta_refund_te=Beta_G_te[[iter]][1:M_it_test[iter,j,iter_out]]
        
        y_h_ellos[iter,j,iter_out] = sum(X_test[j,1:M_it_test[iter,j,iter_out]]*Beta_refund,na.rm = 1)/M_it_test[iter,j,iter_out]
        
        y_h_ellos_te[iter,j,iter_out] = sum(X_test[j,1:M_it_test[iter,j,iter_out]]*Beta_refund_te,na.rm = 1)/M_it_test[iter,j,iter_out]
        
        
      }
      
      # ESTIMATION ERRORS FOR THE TEST DATA SET IN THE CASE OF NORMAL RESPONSE
      
      error_group_VDFR[iter,group,iter_out]=sqrt(sum((y_test-y_h_ellos[iter,,iter_out])^2)/fold)
      error_group_FF_VDFR[iter,group,iter_out]=sqrt(sum((y_test-y_h_sop[iter,,iter_out])^2)/fold)
      error_group_Carmen[iter,group,iter_out]=sqrt(sum((y_test-y_h_ellos_te[iter,,iter_out])^2)/fold)
      
      # ESTIMATION ERRORS FOR THE TEST DATA SET IN THE CASE OF POISSON RESPONSE
      
      # error_group_VDFR[iter,group,iter_out]=sqrt(sum((y_test-(exp(y_h_ellos[iter,,iter_out])))^2)/fold)
      # error_group_Carmen[iter,group,iter_out]=sqrt(sum((y_test-(exp(y_h_ellos_te[iter,,iter_out])))^2)/fold)
      # error_group_FF_VDFR[iter,group,iter_out]=sqrt(sum((y_test-(exp(y_h_sop[iter,,iter_out])))^2)/fold)
      
      
    } # HERE END THE FOR OF group = 1:fold
    
    Y_ERROR_2_ellos[iter,iter_out]= mean(error_group_VDFR[iter,,iter_out])
    Y_ERROR_2_te_ellos[iter,iter_out]= mean(error_group_Carmen[iter,,iter_out])
    Y_ERROR_2_sop[iter,iter_out]= mean(error_group_FF_VDFR[iter,,iter_out])
    
    # HERE WE ESTIMATE THE FUNCTIONAL COEFFICIENT 
    
    formula <- y ~ ffvd(X_se, t, nbasis = c(c1,c2,c3))
    res <- VDPO(formula = formula, data = data)#,family = poisson())

  #   B_data=ffvd(X_se, t, nbasis = c(c1,c2,c3))
  #   
  #   B_all = B_data$B_ffvd
  #   L_Phi <- B_data$L_Phi
  #   B_T <- B_data$B_T
  #   
  #   foo <- VDPO:::B2XZG(B_all, deglist)
  #   X <- foo$X_ffvd
  #   Z <- foo$Z_ffvd
  #   G <- foo$G_ffvd
  #   
  #   fit <- sop.fit(X = X, Z = Z, G = G, y = y, control = list(trace = FALSE))
  #   
  #   theta=c(fit$b.fixed,fit$b.random)
  #   
  #   for (loop_ind in 1:N) {
  # 
  # Beta_FFVD[loop_ind,1:M[loop_ind],iter,iter_out] <- as.matrix(kronecker(L_Phi$B[1:M[loop_ind],], t(B_T$B[i,]))) %*% theta
  #     
  #   }
  #   
    Beta_FFVD[,1:M_max,iter,iter_out]=res$Beta_ffvd[[1]]
    
    error_Beta_FFVD[iter,iter_out]=sum(((Beta[,,iter_out]-Beta_FFVD[,1:M_max,iter,iter_out])^2)/(J*(J+1)), na.rm=TRUE)
    
    Beta_FFVD_inf=Beta_FFVD_sup=matrix(nrow = N,ncol = M_max)
    
    aux_base=ffvd(X_se, t, nbasis = c(c1,c2,c3))
    
    for (aux_ind in 1:N) {
      
      prod=as.matrix(kronecker(aux_base$L_Phi$B[1:M[aux_ind],],t(aux_base$B_T$B[aux_ind,])))
      
      var_curve = diag(prod %*% res$covar_theta %*% t(prod))
      
      std_curve=sqrt(var_curve)
      
      Beta_FFVD_inf[aux_ind,1:M[aux_ind]] =  Beta_FFVD[aux_ind,1:M[aux_ind],iter,iter_out] - 1.96*std_curve
      Beta_FFVD_sup[aux_ind,1:M[aux_ind]] =  Beta_FFVD[aux_ind,1:M[aux_ind],iter,iter_out] + 1.96*std_curve
      
    }
    
    Beta_FFVD_inf_ind[[iter]]=Beta_FFVD_inf
    Beta_FFVD_sup_ind[[iter]]=Beta_FFVD_sup
    
    #####
    
    fit_Gellar_all <- pfr(y ~ lf.vd(X_se,k=89))#,family = poisson())
    
    Beta_G[[iter]]=coef(fit_Gellar_all, n=c(length(t),length(unique(M))))$value
    
    for (j in 1:M_max) {
      
      ind=which(unique(M_it[iter,,iter_out])==M[j])
      ind_t=max(M_it[iter,,iter_out]) # OJO AQU? ESTAMOS SUPONIENDO QUE LAS OBSERVACIONES SON CONSECUTIVAS DESDE 1:MAX(M)
      
      Beta_refund=Beta_G[[iter]][(ind_t*(ind-1)+1):(ind*ind_t)]
      Beta_VD[j,1:M[j],iter,iter_out]=Beta_refund[1:M_it[iter,j,iter_out]]
      
    }
    
    error_Beta_VD[iter,iter_out]=sum(((Beta[,,iter_out]-Beta_VD[,1:M_max,iter,iter_out])^2)/(J*(J+1)), na.rm=TRUE)
    
    
    
  } # HERE ENDS THE MIDDLE FOR ITERATION (RUNS ON R) 
  
  
} # HERE ENDS THE OUTTER FOR ITERATION (RUNS ON iter_out)

end=proc.time()
time=end[1]-start[1]

time/60/60


################# BOXPLOTS
err_B_FFVD=apply(X = error_Beta_FFVD, MARGIN = 2, FUN = mean)
err_B_VD=apply(X = error_Beta_VD, MARGIN = 2, FUN = mean)

err_Y_FFVD=apply(X = Y_ERROR_2_sop, MARGIN = 2, FUN = mean)
err_Y_VD=apply(X = Y_ERROR_2_ellos, MARGIN = 2, FUN = mean)
err_Y_SOF=apply(X = Y_ERROR_2_te_ellos, MARGIN = 2, FUN = mean)

# save.image("C:/Users/user/Desktop/Trabajo/Escuela/Doctorado/Pavel/Tesis/1er Paper/Codigo/Simulation-1st-Paper/Results new method/vacations.RData")

# aux_save=data.frame(Y_ERROR_2_sop[,1],Y_ERROR_2_ellos[,1],Y_ERROR_2_te_ellos[,1],error_Beta_FFVD[1:R,1],error_Beta_VD[1:R,1])

# new_method=data.frame(Y_ERROR_2_sop[,1],Y_ERROR_2_ellos[,1],Y_ERROR_2_te_ellos[,1],error_Beta_FFVD[1:R,1],error_Beta_VD[1:R,1])

B_ERRORES=data.frame(error_Beta_FFVD,error_Beta_VD) #B_ERROR_2,B_ERROR_2_te,

Y_ERRORES=data.frame(Y_ERROR_2_ellos,Y_ERROR_2_te_ellos,Y_ERROR_2_sop) # #Y_ERROR_2,Y_ERROR_2_te,

# write_xlsx(B_ERRORES,path = "C:/Users/user/Desktop/Trabajo/Escuela/Doctorado/Pavel/Tesis/Código/Mix models SOP/Propios/Simulaciones K-fold good version N_200/Poisson_Uniform_N_200/B_ERRORES.xlsx")
# 
# write_xlsx(Y_ERRORES,path = "C:/Users/user/Desktop/Trabajo/Escuela/Doctorado/Pavel/Tesis/Código/Mix models SOP/Propios/Simulaciones K-fold good version N_200/Poisson_Uniform_N_200/Y_ERRORES.xlsx")

Y_values=cbind(c(Y_ERRORES[,2],Y_ERRORES[,10],Y_ERRORES[,18]))
Y_values_1=cbind(c(Y_ERRORES[,1],Y_ERRORES[,3],Y_ERRORES[,5]))
Y_values_2=cbind(c(Y_ERRORES[,2],Y_ERRORES[,4],Y_ERRORES[,6]))
Y_ERRORES_DF=data.frame(values=Y_values,method=as.factor(c(rep("VDFR",R),rep("SOF",R),rep("FF-VDFR",R))))

Y_plot=ggplot(Y_ERRORES_DF,aes(x=method,y=values,fill=method))+
  geom_violin()+
  ylab("Error") +
  scale_x_discrete(labels = NULL, breaks = NULL)+
  xlab("")+
  theme_bw() +
  stat_summary(fun=median, geom="point", size=2, color="black")

Y_plot

B_values=cbind(c(B_ERRORES[,2],B_ERRORES[,10]))
B_values_1=cbind(c(B_ERRORES[,1],B_ERRORES[,3]))
B_values_2=cbind(c(B_ERRORES[,2],B_ERRORES[,4]))
B_ERRORES_DF=data.frame(values=B_values,method=as.factor(c(rep("FF-VDFR",R),rep("VDFR",R))))

B_plot=ggplot(B_ERRORES_DF,aes(x=method,y=values,fill=method))+
  geom_violin()+
  ylab("Error") +
  scale_x_discrete(labels = NULL, breaks = NULL)+
  xlab("")+
  theme_bw() +
  stat_summary(fun=median, geom="point", size=2, color="black")

B_plot

apply(X_s, 1, range,na.rm=TRUE)

case=100
case_iter=which.min(B_ERRORES_DF[1:R,1])

plot(Beta_FFVD[case,,case_iter,iter_out],ylim=c(-2,2))
lines(Beta_FFVD_inf_ind[[case_iter]][case,],col="blue",lty=2,lwd=2)
lines(Beta_FFVD_sup_ind[[case_iter]][case,],col="red",lty=2,lwd=2)
abline(h=0)

# B_ERRORES_test=data.frame(B_ERROR_2_ellos_test,B_ERROR_2_te_ellos_test,B_ERROR_2_sop_test,B_ERROR_2_sop_ad_test) #B_ERROR_2,B_ERROR_2_te,
# 
# B_ERRORES=data.frame(B_ERROR_2_ellos,B_ERROR_2_te_ellos,B_ERROR_2_sop,B_ERROR_2_sop_ad) #B_ERROR_2,B_ERROR_2_te,
# 
# Y_ERRORES=data.frame(Y_ERROR_2_ellos,Y_ERROR_2_te_ellos,Y_ERROR_2_sop,Y_ERROR_2_no_vc) #,Y_ERROR_2_sop_manual) #Y_ERROR_2,Y_ERROR_2_te,
# 
# write_xlsx(B_ERRORES,path = "C:/Users/user/Desktop/Trabajo/Escuela/Doctorado/Pavel/Tesis/Código/Mix models SOP/Propios/Simulaciones K-fold good version N_200/Normal_NegBin_N_200/B_ERRORES.xlsx")
# 
# write_xlsx(Y_ERRORES,path = "C:/Users/user/Desktop/Trabajo/Escuela/Doctorado/Pavel/Tesis/Código/Mix models SOP/Propios/Simulaciones K-fold good version N_200/Normal_NegBin_N_200/Y_ERRORES.xlsx")
# 
# save.image("Normal_NegBin_N_200.Rdata")
# 
# # save(list = ls()[40],file = "Normal_Uniform_N_200_results.Rdata")
# 
# # SELECT THE SCENARIO TO PLOT 
# 
# scenario=8
# 
# VDFR=scenario
# Goldsmith=1*8+scenario
# FF_VDFR=2*8+scenario
# SB_VDFR=3*8+scenario
# 
# Y_values=rbind(as.matrix(Normal_NegBin_N_200_results$Y_ERRORES[,VDFR]),as.matrix(Normal_NegBin_N_200_results$Y_ERRORES[,Goldsmith]),as.matrix(Normal_NegBin_N_200_results$Y_ERRORES[,FF_VDFR]),as.matrix(Normal_NegBin_N_200_results$Y_ERRORES[,SB_VDFR]))
# Y_ERRORES_DF=data.frame(values=Y_values,method=as.factor(c(rep("VDFR",Normal_NegBin_N_200_results$R),rep("Goldsmith",Normal_NegBin_N_200_results$R),rep("FF_VDFR",Normal_NegBin_N_200_results$R),rep("SB_VDFR",Normal_NegBin_N_200_results$R))))
# 
# Y_plot=ggplot(Y_ERRORES_DF,aes(x=method,y=values,fill=method))+
#   geom_violin()+
#   ylab("Error") +
#   scale_x_discrete(labels = NULL, breaks = NULL)+
#   xlab("")+
#   theme_bw() +
#   stat_summary(fun=median, geom="point", size=2, color="black")
# 
# Y_plot
# 
# B_values=rbind(as.matrix(Normal_NegBin_N_200_results$B_ERRORES[,VDFR]),as.matrix(Normal_NegBin_N_200_results$B_ERRORES[,FF_VDFR]))
# B_ERRORES_DF=data.frame(values=B_values,method=as.factor(c(rep("VDFR",Normal_NegBin_N_200_results$R),rep("FF-VDFR",Normal_NegBin_N_200_results$R))))
# 
# B_plot=ggplot(B_ERRORES_DF,aes(x=method,y=values,fill=method))+
#   geom_violin()+
#   ylab("Error") +
#   scale_x_discrete(labels = NULL, breaks = NULL)+
#   xlab("")+
#   theme_bw() +
#   stat_summary(fun=median, geom="point", size=2, color="black")
# 
# B_plot
# 
# ######################
# #############
# ######
# 
# 
# # scenario=8
# # 
# # Gellar=scenario
# # Gellar_te=1*8+scenario
# # New_approach=2*8+scenario
# # NO_VC=3*8+scenario
# # 
# # # aux_b_er=which(B_ERRORES[,New_approach]>0.075)
# # 
# # boxplot(B_ERRORES[,c(Gellar,New_approach)],pars  =  list(xaxt = "n"),xlab="")
# # axis(1, at=c(1,2),gap.axis = 0.75, labels = c("Gellar","New Approach"))
# # 
# # 
# # # aux_er_Gellar=which(Y_ERRORES[,c(Gellar)]>15)
# # # aux_er_Gellar_te=which(Y_ERRORES[,c(Gellar_te)]>15)
# # # aux_er_SOP=which(Y_ERRORES[,c(New_approach)]>15)
# # 
# # # aux=c(aux_er_Gellar_te,aux_er_SOP)
# # 
# # boxplot(Y_ERRORES[,c(Gellar,Gellar_te,New_approach,NO_VC)],pars  =  list(xaxt = "n"),xlab="")
# # axis(1, at=c(1,2,3,4),gap.axis = 0.75, labels = c("VDFR","Goldsmith","FF-VDFR","SB-VDFR"))
# # 
# # 
# # 
# # ############# IN CASE WE NEED SOME METRICS OF THE ERRORS
# # 
# # #  case=8
# # # 
# # # Gellar=case
# # # Gellar_te=1*8+case
# # # New_approach=2*8+case
# # # # Adaptive=3*8+case
# # # 
# # # c(mean(Y_ERRORES[,c(Gellar)]),mean(Y_ERRORES[,c(Gellar_te)]),mean(Y_ERRORES[,c(New_approach)]))
# # # # mean(Y_ERRORES[,c(Adaptive)])
# # # 
# # # c(median(Y_ERRORES[,c(Gellar)]),median(Y_ERRORES[,c(Gellar_te)]), median(Y_ERRORES[,c(New_approach)]))
# # # # median(Y_ERRORES[,c(Adaptive)])
# # # 
# # # c(sd(Y_ERRORES[,c(Gellar)]), sd(Y_ERRORES[,c(Gellar_te)]), sd(Y_ERRORES[,c(New_approach)]))
# # # # sd(Y_ERRORES[,c(Adaptive)])
# # # 
# # # c(mean(B_ERRORES[,c(Gellar)]), mean(B_ERRORES[,c(New_approach)]))
# # # 
# # # # mean(B_ERRORES[,c(Adaptive)])
# # # 
# # # c(median(B_ERRORES[,c(Gellar)]), median(B_ERRORES[,c(New_approach)]))
# # # 
# # # #median(B_ERRORES[,c(Adaptive)])
# # # 
# # # c(sd(B_ERRORES[,c(Gellar)]), sd(B_ERRORES[,c(New_approach)]))
# # #   
# # # #sd(B_ERRORES[,c(Adaptive)])
# # # #
# # 
# # #### MEANS TO AN EXCEL FILE
# # 
# # Y_means=colMeans(Y_ERRORES)
# # B_means=colMeans(B_ERRORES)[1:24]
# # b_matrix=matrix(data = B_means,nrow = 3,byrow = 1)
# # 
# # Y_matrix=matrix(data = Y_means,nrow = 4,byrow = 1)
# # 
# # B_matrix=b_matrix[-2,]
# # 
# # # write_xlsx(as.data.frame(Y_matrix),path = "C:/Users/user/Desktop/Trabajo/Escuela/Doctorado/Pavel/Tesis/Código/Mix models SOP/Propios/Simulaciones/Adaptive/Curvas más suaves/K-fold/R=100/1 y 1/Normal/T_i NegBin/N=100/Despues de IWSM/Y_means.xlsx")
# # 
# # # write_xlsx(as.data.frame(B_matrix),path = "C:/Users/user/Desktop/Trabajo/Escuela/Doctorado/Pavel/Tesis/Código/Mix models SOP/Propios/Simulaciones/Adaptive/Curvas más suaves/K-fold/R=100/1 y 1/Normal/T_i NegBin/N=100/Despues de IWSM/B_means.xlsx")
# # 
# # write_xlsx(as.data.frame(Y_matrix[4,]),path = "C:/Users/user/Desktop/Trabajo/Escuela/Doctorado/Pavel/Tesis/Código/Mix models SOP/Propios/Simulaciones/Adaptive/Curvas más suaves/K-fold/R=100/1 y 1/Poisson/T_i Uniform/N = 500/Y_means_no_vc.xlsx")
# # 
# # 
# # ############ SD TO AN EXCEL FILE
# # 
# # Y_std=matrix(0,nrow = 3,ncol=32)
# # B_std=matrix(0,nrow = 3,ncol=24)
# # q_ric_B=q_ric_Y=NULL
# # N_cases=c(100,200,500)
# # outliers_B=outliers_Y=list()
# # 
# # for (i_ind in 1:3) {
# #   
# #   print(c("i_ind =",i_ind))
# #   
# #   if (i_ind%%3==1) {
# #     
# #     texto="N = 100/Data.RData"
# #     texto_env="N = 100/Data_no_vc.RData"
# #   }
# #   if (i_ind%%3==2) {
# #     
# #     texto="N = 200/Data.RData"
# #     texto_env="N = 200/Data_no_vc.RData"
# #   }
# #   if (i_ind%%3==0) {
# #     
# #     texto="N = 500/Data.RData"
# #     texto_env="N = 500/Data_no_vc.RData"
# #   }
# #   
# #   nam_X <- paste("C:/Users/user/Desktop/Trabajo/Escuela/Doctorado/Pavel/Tesis/Código/Mix models SOP/Propios/Simulaciones/Adaptive/Curvas más suaves/K-fold/R=100/1 y 1/Poisson/T_i Uniform/", texto ,sep="")
# #   
# #   load(nam_X)
# #   
# #   nam_env=paste("env","NegBin_POISSON_N",N_cases[i_ind],sep = "_")
# #   
# #   env_global=assign(nam_env,new.env())
# #   
# #   nam_load=paste("C:/Users/user/Desktop/Trabajo/Escuela/Doctorado/Pavel/Tesis/Código/Mix models SOP/Propios/Simulaciones/Adaptive/Curvas más suaves/K-fold/R=100/1 y 1/Poisson/T_i Uniform/", texto_env, sep="")
# #   
# #   load(nam_load, envir = get(nam_env))
# #   
# #   Y_ERRORES[,25:32]=env_global$Y_ERRORES[,25:32]
# #   
# #   for (iind in 1:(dim(Y_ERRORES)[2])) {
# #     
# #     nam_RIC_Y=paste("RIC_case_Y",iind,sep = "_")
# #     
# #     q_Y=quantile(Y_ERRORES[,iind])
# #     
# #     q_ric_Y[iind]=assign(nam_RIC_Y,q_Y[4]-q_Y[2]) #sqrt(diag(var(Y_ERRORES)))
# #     
# #     outliers_Y[[iind]]=which(Y_ERRORES[,iind]>q_Y[4]+1.5*q_ric_Y[iind] | Y_ERRORES[,iind]<q_Y[2]-1.5*q_ric_Y[iind])
# #     
# #     if (iind<25) {
# #       
# #       nam_RIC_B=paste("RIC_case_B",iind,sep = "_")
# #       
# #       q_B=quantile(B_ERRORES[,iind])
# #       
# #       q_ric_B[iind]=assign(nam_RIC_B,q_B[4]-q_B[2]) #sqrt(diag(var(Y_ERRORES)))
# #       
# #       outliers_B[[iind]]=which(B_ERRORES[,iind]>q_B[4]+1.5*q_ric_B[iind] | B_ERRORES[,iind]<q_B[2]-1.5*q_ric_B[iind])
# #       
# #     }
# #     
# #   }
# #   
# #   for (j_ind in 1:32) {
# #     
# #     print(c("j_ind =",j_ind))
# #     
# #     # case=j_ind
# #     # 
# #     # Gellar=case
# #     # Gellar_te=1*8+case
# #     # New_approach=2*8+case
# #     # NO_VC=3*8+case
# #     
# #     if (length(as.double(outliers_Y[[j_ind]]))==0) {
# #       
# #       Y_std[i_ind,j_ind]=sqrt(var(Y_ERRORES[,j_ind]))
# #       
# #     }else{
# #       
# #       Y_std[i_ind,j_ind]=sqrt(var(Y_ERRORES[-outliers_Y[[j_ind]],j_ind]))
# #     }
# #     
# #     if (j_ind<25) {
# #       
# #       if (length(as.double(outliers_B[[j_ind]]))==0) {
# #         
# #         B_std[i_ind,j_ind]=sqrt(var(B_ERRORES[,j_ind]))
# #         
# #       }else{
# #         B_std[i_ind,j_ind]=sqrt(var(B_ERRORES[-outliers_B[[j_ind]],j_ind]))
# #         
# #       }
# #       
# #     }
# #     
# #   }
# # }
# # 
# # # Y  Gellar Goldsmith SS-VDFR SB-VDFR
# # # B  Gellar Goldsmith SS-VDFR
# # 
# # aux_1=matrix(Y_std[1,],nrow = 4, ncol = 8, byrow = 1)
# # aux_2=matrix(Y_std[2,],nrow = 4, ncol = 8, byrow = 1)
# # aux_3=matrix(Y_std[3,],nrow = 4, ncol = 8, byrow = 1)
# # 
# # Y_matrix=rbind(aux_1, aux_2, aux_3)
# # 
# # aux_1=matrix(B_std[1,],nrow = 3, ncol = 8, byrow = 1)
# # aux_2=matrix(B_std[3,],nrow = 3, ncol = 8, byrow = 1)
# # 
# # B_matrix=rbind(aux_1, aux_2)
# # 
# # write_xlsx(as.data.frame(Y_matrix),path = "C:/Users/user/Desktop/Trabajo/Escuela/Doctorado/Pavel/Tesis/Código/Mix models SOP/Propios/Simulaciones/Adaptive/Curvas más suaves/K-fold/R=100/1 y 1/Poisson/T_i Uniform/Y_sd.xlsx")
# # 
# # write_xlsx(as.data.frame(B_matrix),path = "C:/Users/user/Desktop/Trabajo/Escuela/Doctorado/Pavel/Tesis/Código/Mix models SOP/Propios/Simulaciones/Adaptive/Curvas más suaves/K-fold/R=100/1 y 1/Poisson/T_i Uniform/B_sd.xlsx")
# # 
# # # COMPUTATION TIMES
# # 
# # colMeans(time_SOP)
# # 
# # colMeans(time_no_vc)
# # 
# # colMeans(time_Gellar)
# # 
# # colMeans(time_Goldsmith)
# # 

#### Y

aux_old_Y=c(aux_save$Y_ERROR_2_sop...1.,aux_save$Y_ERROR_2_ellos...1.,aux_save$Y_ERROR_2_te_ellos...1.)

old_version_Y=data.frame(values=aux_old_Y,method=as.factor(c(rep("FF-VDFR",R),rep("VDFR",R),rep("SOF",R))))

Y_plot=ggplot(old_version_Y,aes(x=method,y=values,fill=method))+
  geom_violin()+
  ylab("Error") +
  scale_x_discrete(labels = NULL, breaks = NULL)+
  xlab("")+
  theme_bw() +
  stat_summary(fun=median, geom="point", size=2, color="black")

Y_plot


aux_new_Y=c(new_method$Y_ERROR_2_sop...1.,new_method$Y_ERROR_2_ellos...1.,new_method$Y_ERROR_2_te_ellos...1.)

new_version_Y=data.frame(values=aux_new_Y,method=as.factor(c(rep("FF-VDFR",R),rep("VDFR",R),rep("SOF",R))))

Y_plot=ggplot(new_version_Y,aes(x=method,y=values,fill=method))+
  geom_violin()+
  ylab("Error") +
  scale_x_discrete(labels = NULL, breaks = NULL)+
  xlab("")+
  theme_bw() +
  stat_summary(fun=median, geom="point", size=2, color="black")

Y_plot

#### Beta


aux_old_B=c(aux_save$error_Beta_FFVD.1.R..1.,aux_save$error_Beta_VD.1.R..1.)

old_version_B=data.frame(values=aux_old_B,method=as.factor(c(rep("FF-VDFR",R),rep("VDFR",R))))

B_plot=ggplot(old_version_B,aes(x=method,y=values,fill=method))+
  geom_violin()+
  ylab("Error") +
  scale_x_discrete(labels = NULL, breaks = NULL)+
  xlab("")+
  theme_bw() +
  stat_summary(fun=median, geom="point", size=2, color="black")

B_plot


aux_new_B=c(new_method$error_Beta_FFVD.1.R..1.,new_method$error_Beta_VD.1.R..1.)

new_version_B=data.frame(values=aux_new_B,method=as.factor(c(rep("FF-VDFR",R),rep("VDFR",R))))

B_plot=ggplot(new_version_B,aes(x=method,y=values,fill=method))+
  geom_violin()+
  ylab("Error") +
  scale_x_discrete(labels = NULL, breaks = NULL)+
  xlab("")+
  theme_bw() +
  stat_summary(fun=median, geom="point", size=2, color="black")

B_plot

case=90

plot(X_s[case,])
lines(X_se[case,])






