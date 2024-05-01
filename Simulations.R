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

N=100 # NUMBER OF SUBJECTS

J=100 # NUMBER OF MAXIMUM OBSERVATIONS PER SUBJECTS

k=10  # SIZE OF THE TEST DATA SET 

fold=N/k

R=100 # NUMBER OF ITERATIONS FOR THE SIMULATION STUDY 

case=2 # THIS IS THE TOTAL NUMBER OF SCENARIOS SIMULATED 

c1=30
c2=30
c3=30

Rsq=0.95

registrated_points=function(x) {
  J.i <- sum(!is.na(x))
  y <- approx(x=seq(0,1,length=J.i), y=x[1:J.i],
              xout=seq(0,1,length=J))$y
  y
}

time_SOP=time_Gellar=time_Goldsmith=array(dim = c(R,case))

#ndb=15 # THIS IS THE NUMBER OF LAMBDAS FOR THE ADPATIVE CASE

# # HERE WE GENERATE TO EMPTY ARRAYS THAT WILL BE THE ESTIMATED FUNCTIONAL COEFFICIENTS
Beta_VD=Beta_FFVD=array(dim =c(N,J,R,case))


error_Beta_FFVD=error_Beta_VD=array(dim =c(R,case))
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

#

# EMPTY ARRAYS WHERE THE ERROR OF THE ESTIMATED COEFFICIENT FUNCTION WILL BE STORE FOR EVERY ITERATION IN EVERY ESCANARIO
B_ERROR_2_ellos=B_ERROR_2_te_ellos=array(dim = c(R,case))
B_ERROR_2_sop=array(dim = c(R,case))
B_ERROR_2_ellos_test=B_ERROR_2_te_ellos_test=array(dim = c(R,case))
B_ERROR_2_sop_test=array(dim = c(R,case))
#


# THIS ARE THE N/K ERRORS IN EACH GROUPS
error_group_FF_VDFR=error_group_SB=error_group_Carmen=error_group_VDFR=array(dim=c(R,k,case))

# EMPTY ARRAYS WHERE THE ERROR OF THE ESTIMATED RESPOND VARIABLE WILL BE STORE FOR EVERY ITERATION IN EVERY ESCANARIO
Y_ERROR_2_ellos=Y_ERROR_2_te_ellos=array(dim = c(R,case))
Y_ERROR_2_sop=array(dim = c(R,case))
#

start=proc.time() 

for (iter_out in 1:case) { # HERE WE CAN SELECTED WHICH SCENARIO(S) SIMULATE 
  
  print(c("case = ",iter_out))
  
  for (iter in 1:R) {
    
    print(c("iter = ",iter))
    
    set.seed(1000+iter) 
    
    # M = round(runif(N,10,J),digits = 0) # HERE WE GENERATE THE DOMAIN FOR ALL SUBJECTS WITH A MINIMUM OF A 10 OBSERVATIONS
    
    M = rnegbin(N,24,1) # Para 1000 poner (240,2)
    
    if (max(M)>J) {
      # print("si")
      M[which(M>J)]=J
    }
    
    
    if (min(M)<=10) {
      # print("si")
      # M[which(M<=10)]=round(runif(length(which(M<=10)),31,max(M)))
      M[which(M<=10)]=10
    }
    
    M_max=max(M)
    
    t=seq(from = 0, to = 1, length.out=J)
    
    M=sort(M) # WE SORT THE DATA WITHOUT LOSS OF GENERALITY
    
    M_it[iter,,iter_out]=M # WE STORE THE DOMAINS FOR EVERY ITERATION AND SCENARIO
    
    
    ############ HERE WE GENERATE THE FUNCTIONAL DATA
    
    X_se=matrix(NA,N,M_max) # NOISY
    
    X_s=matrix(NA,N,M_max) # NOT NOISY
    
    
    for (i in 1:N) {
      
      u=rnorm(1)
      
      temp=matrix(NA,5,M_max)
      
      for (aux in 1:5) {
        
        v_i1=rnorm(1,0,4/aux^2)
        
        
        temp[aux,1:M[i]]=v_i1*sin(2*pi*aux*(1:M[i])/M_max)
      }
      
      B=apply(temp,2,sum)
      
      B=B+u
      
      X_s[i,]=B 
      X_se[i,]=(B)+rnorm(M_max,0,2.5) # WE ADD NOISE
      
    }
    
    X_reg = t(apply(X_se, 1,registrated_points)) # THESE ARE THE REGISTRATED POINTS
    
    # SOME SAVE CHECKS FOR UNWANTED NAs
    for (i in 1:M_max) {
      if (length(which(is.na(X_s[,i])))==N) {
        print(c("iteracion",l,"columna",i))
        stop("Hay columnas con NA")  
      }}
    
    
    
    
    ###### HERE WE GENERATE THE TRUE Beta COEFFICIENT AND THE RESPONSE VARIABLE
    
    Beta=array(dim = c(N,M_max,2))  
    nu=y=rep(0,N)
    
    for (i in 1:N) {
      
      # TRUE FUNCTIONAL COEFFICIENTS
      
      Beta[i,1:(M[i]),1]=2*sin(0.5*pi*(t[1:M[i]]-M[i])/J)+4*sin(1.5*pi*(t[1:M[i]]-M[i])/J)+5*sin(2.5*pi*t[1:M[i]])
      Beta[i,1:(M[i]),2]=exp(-2*((t[1:M[i]]-0.75))) - exp(-2*((t[1:M[i]]-0.5)))
      
      # HERE WE GENERATE THE RESPONSE VARIABLE FOR EVERY BETA FROM THE 2 DIFFERENT FUNCTIONAL DATA (NOISY AND OTHERWISE)  
      
      
      if (iter_out==2) {
        nu[i]=sum(X_s[i,]*Beta[i,,2],na.rm = 1)/(M[i]) #NOT NOISY
      }
      if (iter_out==1) {
        nu[i]=sum(X_s[i,]*Beta[i,,1],na.rm = 1)/(M[i]) #NOT NOISY
      }
      
    }
    
    
    var_e <- (1 / Rsq - 1) * stats::var(nu) # (1-Rsq)*var(nu[ind,])
    
    y=nu+rnorm(N,sd = var_e) # ADDING NOISE TO THE GAUSSIAN MODEL 
    
    # y=rpois(N,exp(nu)) # POISSON MODEL ##### CAHNGE THE ERROR EQUATIONS IN LINE 477
    
    # y=rbinom(N,1,(exp(nu)/(1+exp(nu)))) # BINOMIAL MODEL
    
    ###  HERE WE STORE THE DATA INTO A DATA FRAME FOR THE VDPO package
    
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
      
      start_SOP=proc.time()
      
      # BB=Data2B_simpson(X_train, M_train, nbasis=c(c1,c2,c3),sub = 25, lim =c(min(M),max(M)),) # HERE WE TRASFORM THE FUNCTIONAL MODEL INTO A MULTIVARIATE MODEL
      # E=B2XZG(BB$B,c=c(c2,c3)) # HERE THE FUNCTION B2XZG() GENERATE THE NECESSARY MATRICES FOR THE MIXED MODEL
      # res=XZG2theta(E$X, E$Z, E$G, E$T, y_train,family = poisson()) # HERE THE MODEL COEFICIENT ARE ESTIMATED.
      
      data_train=data.frame(y_train)
      data_train[["X_train"]]=X_train
      
      data_test=data.frame(y_test)
      data_test[["X_test"]]=X_test
      
      formula <- y_train ~ ffvd(X_train,nbasis = c(c1,c2,c3))
      BB=ffvd(X_train,nbasis = c(c1,c2,c3))
      res <- VDPO(formula = formula, data = data_train) #family = poisson()
      
      
      end_SOP=proc.time()
      time_SOP[iter,iter_out]=end_SOP[1]-start_SOP[1]
      #
      
      # MODEL ESTIMATION USING GELLAR APPROACH
      
      start_Gellar=proc.time()
      
      fit_Gellar <- pfr(y_train ~ lf.vd(X_train,k=89)) #,family = poisson()
      
      end_Gellar=proc.time()
      time_Gellar[iter,iter_out]=end_Gellar[1]-start_Gellar[1]
      #
      
      # MODEL ESTIMATION USING GOLDSMITH APPROACH
      
      start_Goldsmith=proc.time()
      
      fit_Gellar_te <- pfr(y_train ~ lf(X_reg_train, bs = "ps",k=25,presmooth = "bspline",presmooth.opts = list(nbasis=25)))#,family = poisson()
      
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
      
      BB_test=ffvd(X_test,nbasis = c(c1,c2,c3))
      
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
    
    formula <- y ~ ffvd(X_se,nbasis = c(c1,c2,c3))
    res <- VDPO(formula = formula, data = data) #family = poisson()
    
    Beta_FFVD[,1:M_max,iter,iter_out]=res$Beta_ffvd[[1]]
    
    error_Beta_FFVD[iter,iter_out]=sum(((Beta[,,iter_out]-Beta_FFVD[,1:M_max,iter,iter_out])^2)/(J*(J+1)), na.rm=TRUE)
    
    #####
    
    Beta_G[[iter]]=coef(fit_Gellar, n=c(length(t),length(unique(M))))$value

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


B_ERRORES=data.frame(B_ERROR_2_ellos,B_ERROR_2_te_ellos,B_ERROR_2_sop) #B_ERROR_2,B_ERROR_2_te,

Y_ERRORES=data.frame(Y_ERROR_2_ellos,Y_ERROR_2_te_ellos,Y_ERROR_2_sop) # #Y_ERROR_2,Y_ERROR_2_te,

write_xlsx(B_ERRORES,path = "C:/Users/user/Desktop/Trabajo/Escuela/Doctorado/Pavel/Tesis/Código/Mix models SOP/Propios/Simulaciones K-fold good version N_200/Poisson_Uniform_N_200/B_ERRORES.xlsx")

write_xlsx(Y_ERRORES,path = "C:/Users/user/Desktop/Trabajo/Escuela/Doctorado/Pavel/Tesis/Código/Mix models SOP/Propios/Simulaciones K-fold good version N_200/Poisson_Uniform_N_200/Y_ERRORES.xlsx")

save.image("Poisson_Uniform_N_200.Rdata")

# SELECT THE SCENARIO TO PLOT 

scenario=8

VDFR=scenario
Goldsmith=1*8+scenario
FF_VDFR=2*8+scenario
SB_VDFR=3*8+scenario

Y_values=rbind(as.matrix(Poisson_Uniform_N_200_results$Y_ERRORES[,VDFR]),as.matrix(Poisson_Uniform_N_200_results$Y_ERRORES[,Goldsmith]),as.matrix(Poisson_Uniform_N_200_results$Y_ERRORES[,FF_VDFR]),as.matrix(Poisson_Uniform_N_200_results$Y_ERRORES[,SB_VDFR]))
Y_ERRORES_DF=data.frame(values=Y_values,method=as.factor(c(rep("VDFR",Poisson_Uniform_N_200_results$R),rep("Goldsmith",Poisson_Uniform_N_200_results$R),rep("FF_VDFR",Poisson_Uniform_N_200_results$R),rep("SB_VDFR",Poisson_Uniform_N_200_results$R))))

Y_plot=ggplot(Y_ERRORES_DF,aes(x=method,y=values,fill=method))+
  geom_violin()+
  ylab("Error") +
  scale_x_discrete(labels = NULL, breaks = NULL)+
  xlab("")+
  theme_bw() +
  stat_summary(fun=median, geom="point", size=2, color="black")

Y_plot

B_values=rbind(as.matrix(Poisson_Uniform_N_200_results$B_ERRORES[,VDFR]),as.matrix(Poisson_Uniform_N_200_results$B_ERRORES[,FF_VDFR]))
B_ERRORES_DF=data.frame(values=B_values,method=as.factor(c(rep("VDFR",Poisson_Uniform_N_200_results$R),rep("FF-VDFR",Poisson_Uniform_N_200_results$R))))

B_plot=ggplot(B_ERRORES_DF,aes(x=method,y=values,fill=method))+
  geom_violin()+
  ylab("Error") +
  scale_x_discrete(labels = NULL, breaks = NULL)+
  xlab("")+
  theme_bw() +
  stat_summary(fun=median, geom="point", size=2, color="black")

B_plot

