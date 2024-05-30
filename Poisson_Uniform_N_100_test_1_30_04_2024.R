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
options(error=NULL)
library(remotes)
# install_github("Pavel-Hernadez-Amaro/VDPO")
library(VDPO)

########### Here we generate the data and set some fixed parameters

N=100 # NUMBER OF SUBJECTS

J=100 # NUMBER OF MAXIMUM OBSERVATIONS PER SUBJECTS

k=10  # SIZE OF THE TEST DATA SET 

fold=N/k

R=1 # NUMBER OF ITERATIONS FOR THE SIMULATION STUDY 

case=2 # THIS IS THE TOTAL NUMBER OF SCENARIOS SIMULATED 

registrated_points=function(x) {
  J.i <- sum(!is.na(x))
  y <- approx(x=seq(0,1,length=J.i), y=x[1:J.i],
              xout=seq(0,1,length=J))$y
  y
}

time_no_vc=time_SOP=time_Gellar=time_Goldsmith=array(dim = c(R,case))

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

for (iter_out in 1:case) { # HERE WE CAN SELECTED WHICH SCENARIO(S) SIMULATE 
  
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
      
      for (k in 1:5) {
        
        v_i1=rnorm(1,0,4/k^2)
        
        
        temp[k,1:M[i]]=v_i1*sin(2*pi*k*(1:M[i])/M_max)
      }
      
      B=apply(temp,2,sum)
      
      B=B+u
      
      X_s[i,]=B 
      X_se[i,]=(B)+rnorm(M_max,0,1) # WE ADD NOISE
      
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
    
    
    y=nu+rnorm(N,sd = 1) # ADDING NOISE TO THE GAUSSIAN MODEL 
    
    # y=rpois(N,exp(nu)) # POISSON MODEL ##### CAHNGE THE ERROR EQUATIONS IN LINE 477
    
    # y=rbinom(N,1,(exp(nu)/(1+exp(nu)))) # BINOMIAL MODEL
    
    ############## HERE ENDS THE GENERATION OF THE DATA  
    
    ############ HERE BEGINS THE ESTIMATION OF THE MODEL
    
    # ASSIGNING THE TRAIN AND TEST SETS 
    
    for (group in 1:k) {
      
      print(c("group = ",group))
      
      current_group=(fold*(group-1)+1):(group*fold)
      
      X_test=X_se[current_group,]
      X_train=X_se[-current_group,]
      
      X_reg_train=X_reg[-current_group,]
      X_reg_test=X_reg[current_group,]
      
      M_test=M[current_group]
      M_train=M[-current_group]
      
      M_it_test[iter,,iter_out]=M_test
      M_it_train[iter,,iter_out]=M_train
      
      y_test=y[current_group]
      y_train=y[-current_group]
      
      # NUMBER OF BASIS FOR EVERY MARGINAL FUNCTION IN OUR APPROACH
      c1=25
      c2=25
      c3=25
      
      # MODEL ESTIMATION BY OUR APPROACH
      
      start_SOP=proc.time()
      
      # BB=Data2B_simpson(X_train, M_train, nbasis=c(c1,c2,c3),sub = 25, lim =c(min(M),max(M)),) # HERE WE TRASFORM THE FUNCTIONAL MODEL INTO A MULTIVARIATE MODEL
      # E=B2XZG(BB$B,c=c(c2,c3)) # HERE THE FUNCTION B2XZG() GENERATE THE NECESSARY MATRICES FOR THE MIXED MODEL
      # res=XZG2theta(E$X, E$Z, E$G, E$T, y_train,family = poisson()) # HERE THE MODEL COEFICIENT ARE ESTIMATED.
      
      data=data.frame(y=y)
      data[["X_se"]]=X_se

      data_train=data.frame(y_train)
      data_train[["X_train"]]=X_train
      
      data_test=data.frame(y_test)
      data_test[["X_test"]]=X_test
      
      formula <- y_train ~ ffvd(X_train)
      res <- VDPO(formula = formula, data = data_train) #family = poisson()
      
      
      end_SOP=proc.time()
      time_SOP[iter,iter_out]=end_SOP[1]-start_SOP[1]
      #
      
      start_no_vc=proc.time() 
      
      BB_EPOC_no_vc=Data2B_simpson_no_vc(X_train, M_train, nbasis=c(c1,c2),sub = 25) # lim ES NULO PORQUE NO ES UNA SIMULACION
      E_EPOC_no_vc=B2XZG_1d(BB_EPOC_no_vc$B,c=c(c2)) # AQUÍ UTILIZO LA FUNCIÓN B2XZG() PARA GENERAR LAS MATRICES NECESARIAS EN UN MODELO MIXTO
      res_EPOC_no_vc=XZG2theta(X = E_EPOC_no_vc$X, Z = E_EPOC_no_vc$Z, G = E_EPOC_no_vc$G, T = E_EPOC_no_vc$T, y = y_train, family = poisson())
      
      end_no_vc=proc.time()
      time_no_vc[iter,iter_out]=end_no_vc[1]-start_no_vc[1]
      #
      
      # ADAPTIVE
      # BB_ad=Data2B_simpson_ad(X_train, M_train, nbasis=c(c1,c2,c3),sub = 25, ndb = ndb)
      # E_ad=B2XZG(BB_ad$B,c=c(c2,c3)) # AQU? UTILIZO LA FUNCI?N B2XZG() PARA GENERAR LAS MATRICES NECESARIAS EN UN MODELO MIXTO
      # res_ad=XZG2theta(E_ad$X, E_ad$Z, E_ad$G, E_ad$T, y_train) # AQU? AJUSTO EL MODELO MIXTO Y RECUPERO LOS COEFICIENTES ORIGINALES
      
      # MODEL ESTIMATION USING GELLAR APPROACH
      
      start_Gellar=proc.time()
      
      fit_Gellar <- pfr(y_train ~ lf.vd(X_train,k=89),family = poisson())
      
      end_Gellar=proc.time()
      time_Gellar[iter,iter_out]=end_Gellar[1]-start_Gellar[1]
      #
      
      # MODEL ESTIMATION USING GOLDSMITH APPROACH
      
      start_Goldsmith=proc.time()
      
      fit_Gellar_te <- pfr(y_train ~ lf(X_reg_train, bs = "ps",k=25,presmooth = "bspline",presmooth.opts = list(nbasis=25)),family = poisson())#, offset = log((M_train)/365)      
      
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
      
      # ESTIMATED COEFFICIENTS FOR THE GELLAR AND GOLDSMITH APPROACHES 
      
      Beta_G[[iter]]=coef(fit_Gellar, n=c(length(t),length(unique(M))))$value
      Beta_G_te[[iter]]=coef(fit_Gellar_te, n=M_max)$value
      
      for (j in 1:((k-1)*fold)) { # THESE ARE THE ESTIMATED COEFFICIENTS BY MY APPROACH
        
        prod=as.matrix(kronecker(BB$B_Phi[[j]]$B,t(BB$B_T$B[j,])))
        
        Beta_h[j,1:M_it_train[iter,j,iter_out],iter]=as.vector(prod %*% res$theta) # HERE WE MULTIPLY THE BIDIMENSIONAL BASIS FOR THE COEFFICIENT
        
        # HERE WE SELECT THE CORRECT INDEX FOR THE COEFFICIENT OF THE GELLAR AND GOLDSMITH APPROACHES
        
        ind=which(unique(M_it[iter,,iter_out])== M_it_train[iter,j,iter_out])
        ind_t=max(M_it[iter,,iter_out]) # HERE IS CRITICAL THAT OBSERVATIONS ARE CONSECUTIVES IN 1:MAX(M)
        Beta_refund=Beta_G[[iter]][(ind_t*(ind-1)+1):(ind*ind_t)]
        Beta_refund=Beta_refund[1:M_it_train[iter,j,iter_out]]
        # Beta_refund_te=Beta_G_te[[iter]][(ind_t*(ind-1)+1):(ind*ind_t)]
        # Beta_refund_te=Beta_refund_te[1:M_it_train[iter,j,iter_out]]
        #
        
        
        # HERE WE CALCULATE THE ERROR OF THE ESTIMATED FUNCTIONAL COEFFICIENTS
        
        if (iter_out<=4) {
          
          True_Beta=Beta[-(current_group),1:M_it_train[iter,j,iter_out],iter_out]
          
          error_2_sop[j,1]=sum((True_Beta[j,]-Beta_h[j,1:M_it_train[iter,j,iter_out],iter])^2)
          
          error_2_ellos[j,1]=sum((True_Beta[j,]-Beta_refund)^2)
          
          # error_2_te_ellos[j,1]=sum((Beta[train[j],1:M_it_train[iter,j,iter_out],iter_out]-Beta_refund_te)^2)
          
        }
        
        if (iter_out>4 & iter_out<8) {
          
          True_Beta=Beta[(-current_group),1:M_it_train[iter,j,iter_out],iter_out%%4]
          
          error_2_sop[j,1]=sum((True_Beta[j,]-Beta_h[j,1:M_it_train[iter,j,iter_out],iter])^2)
          
          error_2_ellos[j,1]=sum((True_Beta[j,]-Beta_refund)^2)
          
          # error_2_te_ellos[j,1]=sum((Beta[train[j],1:M_it_train[iter,j,iter_out],iter_out%%4]-Beta_refund_te)^2)
          
        }
        
        if (iter_out==8) {
          
          True_Beta=Beta[(-current_group),1:M_it_train[iter,j,iter_out],4]
          
          error_2_sop[j,1]=sum((True_Beta[j,]-Beta_h[j,1:M_it_train[iter,j,iter_out],iter])^2)
          
          error_2_ellos[j,1]=sum((True_Beta[j,]-Beta_refund)^2)
          
          # error_2_te_ellos[j,1]=sum((Beta[train[j],1:M_it_train[iter,j,iter_out],4]-Beta_refund_te)^2)
          
        }
        
        
        
      }
      
      ERROR_2_ellos=sum(error_2_ellos[,1])/(J*(J+1))
      
      ERROR_2_sop=sum(error_2_sop[,1])/(J*(J+1))
      
      
      #############
      
      B_ERROR_2_ellos[iter,iter_out]=ERROR_2_ellos
      
      B_ERROR_2_sop[iter,iter_out]=ERROR_2_sop
      
      ########### END OF ESTIMATION ERRORS
      
      
      ############# HERE WE CALCULATE THE ERROR OF THE ESTIMATED RESPONSE VARIABLE
      
      BB_test=Data2B_simpson(X_test, M_test, nbasis=c(c1,c2,c3),sub = 25,lim = c(min(M),max(M))) # WE GENERATE THE BASIS OF THE TEST DATA SET
      
      y_h_sop[iter,,iter_out] = BB_test$B %*% res$theta # ESTIMATED REPSONSE VARIABLE USING OUR APPROACH
      
      
      BB_test_no_vc=Data2B_simpson_no_vc(X_test, M_test, nbasis=c(c1,c2),sub = 25,lim = c(min(M),max(M))) # WE GENERATE THE BASIS OF THE TEST DATA SET
      
      y_h_no_vc[iter,,iter_out] = BB_test_no_vc$B %*% res_EPOC_no_vc$theta 
      
      
      # ESTIMATED REPSONSE VARIABLE USING HE GELLAR AND GOLDSMITH APPROACHES
      
      for (j in 1:fold) {
        
        prod=as.matrix(kronecker(BB_test$B_Phi[[j]]$B,t(BB_test$B_T$B[j,])))
        
        Beta_h_test[j,1:M_it_test[iter,j,iter_out],iter]=as.vector(prod %*% res$theta) # HERE WE MULTIPLY THE BIDIMENSIONAL BASIS FOR THE COEFFICIENT
        
        ind=which(unique(M_it[iter,,iter_out])== M_it_test[iter,j,iter_out])
        ind_t=max(M_it[iter,,iter_out]) # OJO AQU? ESTAMOS SUPONIENDO QUE LAS OBSERVACIONES SON CONSECUTIVAS DESDE 1:MAX(M)
        
        Beta_refund=Beta_G[[iter]][(ind_t*(ind-1)+1):(ind*ind_t)]
        Beta_refund=Beta_refund[1:M_it_test[iter,j,iter_out]]
        
        Beta_refund_te=Beta_G_te[[iter]][1:M_it_test[iter,j,iter_out]]
        
        y_h_ellos[iter,j,iter_out] = sum(X_test[j,1:M_it_test[iter,j,iter_out]]*Beta_refund,na.rm = 1)/M_it_test[iter,j,iter_out]
        
        y_h_ellos_te[iter,j,iter_out] = sum(X_test[j,1:M_it_test[iter,j,iter_out]]*Beta_refund_te,na.rm = 1)/M_it_test[iter,j,iter_out]
        
        
        # THIS SECTION OF CODE IS FOR THE CASE WHEN THE USER WANTS TO CALCULATE THE ERRORS FOR THE ESTIMATION OF THE FUNCTIONAL COEFFICIENT ONLY IN THE TEST SET
        
        #
        # if (iter_out<=4) {
        #
        #   error_2_sop_test[j,1]=sum((Beta[test[j],1:M_it_test[iter,j,iter_out],iter_out]-Beta_h_test[j,1:M_it_test[iter,j,iter_out],iter])^2)
        #
        #   error_2_sop_ad_test[j,1]=sum((Beta[test[j],1:M_it_test[iter,j,iter_out],iter_out]-Beta_h_ad_test[j,1:M_it_test[iter,j,iter_out],iter])^2)
        #
        #   # error_2[j,1]=sum((Beta[j,1:M[j],4]-Beta_refund_Carmen[1:M[j]])^2)
        #   #
        #   # error_2_te[j,1]=sum((Beta[j,1:M[j],4]-Beta_refund_Carmen_te[1:M[j]])^2)
        #
        #   error_2_ellos_test[j,1]=sum((Beta[test[j],1:M_it_test[iter,j,iter_out],iter_out]-Beta_refund)^2)
        #
        #   error_2_te_ellos_test[j,1]=sum((Beta[test[j],1:M_it_test[iter,j,iter_out],iter_out]-Beta_refund_te)^2)
        #
        # }
        #
        # if (iter_out>4 & iter_out<8) {
        #
        #   error_2_sop_test[j,1]=sum((Beta[test[j],1:M_it_test[iter,j,iter_out],iter_out%%4]-Beta_h_test[j,1:M_it_test[iter,j,iter_out],iter])^2)
        #
        #   error_2_sop_ad_test[j,1]=sum((Beta[test[j],1:M_it_test[iter,j,iter_out],iter_out%%4]-Beta_h_ad_test[j,1:M_it_test[iter,j,iter_out],iter])^2)
        #
        #   # error_2[j,1]=sum((Beta[j,1:M[j],4]-Beta_refund_Carmen[1:M[j]])^2)
        #   #
        #   # error_2_te[j,1]=sum((Beta[j,1:M[j],4]-Beta_refund_Carmen_te[1:M[j]])^2)
        #
        #   error_2_ellos_test[j,1]=sum((Beta[test[j],1:M_it_test[iter,j,iter_out],iter_out%%4]-Beta_refund)^2)
        #
        #   error_2_te_ellos_test[j,1]=sum((Beta[test[j],1:M_it_test[iter,j,iter_out],iter_out%%4]-Beta_refund_te)^2)
        #
        # }
        #
        # if (iter_out==8) {
        #
        #   error_2_sop_test[j,1]=sum((Beta[test[j],1:M_it_test[iter,j,iter_out],4]-Beta_h_test[j,1:M_it_test[iter,j,iter_out],iter])^2)
        #
        #   error_2_sop_ad_test[j,1]=sum((Beta[test[j],1:M_it_test[iter,j,iter_out],4]-Beta_h_ad_test[j,1:M_it_test[iter,j,iter_out],iter])^2)
        #
        #   # error_2[j,1]=sum((Beta[j,1:M[j],4]-Beta_refund_Carmen[1:M[j]])^2)
        #   #
        #   # error_2_te[j,1]=sum((Beta[j,1:M[j],4]-Beta_refund_Carmen_te[1:M[j]])^2)
        #
        #   error_2_ellos_test[j,1]=sum((Beta[test[j],1:M_it_test[iter,j,iter_out],4]-Beta_refund)^2)
        #
        #   error_2_te_ellos_test[j,1]=sum((Beta[test[j],1:M_it_test[iter,j,iter_out],4]-Beta_refund_te)^2)
        #
        # }
        #
        # ERROR_2_ellos_test=sum(error_2_ellos_test[,1])/(J*(J+1))
        #
        # ERROR_2_te_ellos_test=sum(error_2_te_ellos_test[,1])/(J*(J+1))
        #
        # ERROR_2_sop_test=sum(error_2_sop_test[,1])/(J*(J+1))
        #
        # # ERROR_2_sop_ad_test=sum(error_2_sop_ad_test[,1])/(J*(J+1))
        #
        # #############
        #
        # B_ERROR_2_ellos_test[iter,iter_out]=ERROR_2_ellos_test
        #
        # B_ERROR_2_te_ellos_test[iter,iter_out]=ERROR_2_te_ellos_test
        #
        # B_ERROR_2_sop_test[iter,iter_out]=ERROR_2_sop_test
        #
        # # B_ERROR_2_sop_ad_test[iter,iter_out]=ERROR_2_sop_ad_test
        
        
      } # HERE ENDS THE INNNER FOR ITERATION (RUNS ON K)
      
      
      # ESTIMATION ERRORS FOR THE TEST DATA SET IN THE CASE OF NORMAL RESPONSE
      
      error_group_VDFR[iter,group,iter_out]=sqrt(sum((y_test-y_h_ellos[iter,,iter_out])^2)/fold)
      error_group_FF_VDFR[iter,group,iter_out]=sqrt(sum((y_test-y_h_sop[iter,,iter_out])^2)/fold)
      error_group_SB[iter,group,iter_out]=sqrt(sum((y_test-y_h_no_vc[iter,,iter_out])^2)/fold)
      error_group_Carmen[iter,group,iter_out]=sqrt(sum((y_test-y_h_ellos_te[iter,,iter_out])^2)/fold)
      
      # ESTIMATION ERRORS FOR THE TEST DATA SET IN THE CASE OF POISSON RESPONSE
      
      # error_group_VDFR[iter,group,iter_out]=sqrt(sum((y_test-round(exp(y_h_ellos[iter,,iter_out])))^2)/fold)
      # error_group_Carmen[iter,group,iter_out]=sqrt(sum((y_test-round(exp(y_h_ellos_te[iter,,iter_out])))^2)/fold)
      # error_group_FF_VDFR[iter,group,iter_out]=sqrt(sum((y_test-round(exp(y_h_sop[iter,,iter_out])))^2)/fold)
      # error_group_SB[iter,group,iter_out]=sqrt(sum((y_test-round(exp(y_h_no_vc[iter,,iter_out])))^2)/fold)
      
      
    } # HERE END THE FOR OF group = 1:fold
    
    Y_ERROR_2_ellos[iter,iter_out]= mean(error_group_VDFR[iter,,iter_out])
    Y_ERROR_2_te_ellos[iter,iter_out]= mean(error_group_Carmen[iter,,iter_out])
    Y_ERROR_2_sop[iter,iter_out]= mean(error_group_FF_VDFR[iter,,iter_out])
    Y_ERROR_2_no_vc[iter,iter_out]= mean(error_group_SB[iter,,iter_out])
    
    
  } # HERE ENDS THE MIDDLE FOR ITERATION (RUNS ON R) 
  
} # HERE ENDS THE OUTTER FOR ITERATION (RUNS ON iter_out)

end=proc.time()
time=end[1]-start[1]

time/60/60


################# BOXPLOTS

B_ERRORES_test=data.frame(B_ERROR_2_ellos_test,B_ERROR_2_te_ellos_test,B_ERROR_2_sop_test,B_ERROR_2_sop_ad_test) #B_ERROR_2,B_ERROR_2_te,

B_ERRORES=data.frame(B_ERROR_2_ellos,B_ERROR_2_te_ellos,B_ERROR_2_sop,B_ERROR_2_sop_ad) #B_ERROR_2,B_ERROR_2_te,

Y_ERRORES=data.frame(Y_ERROR_2_ellos,Y_ERROR_2_te_ellos,Y_ERROR_2_sop,Y_ERROR_2_no_vc) #,Y_ERROR_2_sop_manual) #Y_ERROR_2,Y_ERROR_2_te,

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


######################
#############
######

case=1

Gellar=case
Gellar_te=1*8+case
New_approach=2*8+case
NO_VC=3*8+case

# aux_b_er=which(B_ERRORES[,New_approach]>0.075)

boxplot(B_ERRORES[,c(Gellar,New_approach)],pars  =  list(xaxt = "n"),xlab="")
axis(1, at=c(1,2),gap.axis = 0.75, labels = c("Gellar","New Approach"))


# aux_er_Gellar=which(Y_ERRORES[,c(Gellar)]>15)
# aux_er_Gellar_te=which(Y_ERRORES[,c(Gellar_te)]>15)
# aux_er_SOP=which(Y_ERRORES[,c(New_approach)]>15)

# aux=c(aux_er_Gellar_te,aux_er_SOP)

boxplot(Y_ERRORES[,c(Gellar,Gellar_te,New_approach,NO_VC)],pars  =  list(xaxt = "n"),xlab="")
axis(1, at=c(1,2,3,4),gap.axis = 0.75, labels = c("VDFR","Goldsmith","FF-VDFR","SB-VDFR"))
# 
# 
# 
# ############# IN CASE WE NEED SOME METRICS OF THE ERRORS
# 
# #  case=8
# # 
# # Gellar=case
# # Gellar_te=1*8+case
# # New_approach=2*8+case
# # # Adaptive=3*8+case
# # 
# # c(mean(Y_ERRORES[,c(Gellar)]),mean(Y_ERRORES[,c(Gellar_te)]),mean(Y_ERRORES[,c(New_approach)]))
# # # mean(Y_ERRORES[,c(Adaptive)])
# # 
# # c(median(Y_ERRORES[,c(Gellar)]),median(Y_ERRORES[,c(Gellar_te)]), median(Y_ERRORES[,c(New_approach)]))
# # # median(Y_ERRORES[,c(Adaptive)])
# # 
# # c(sd(Y_ERRORES[,c(Gellar)]), sd(Y_ERRORES[,c(Gellar_te)]), sd(Y_ERRORES[,c(New_approach)]))
# # # sd(Y_ERRORES[,c(Adaptive)])
# # 
# # c(mean(B_ERRORES[,c(Gellar)]), mean(B_ERRORES[,c(New_approach)]))
# # 
# # # mean(B_ERRORES[,c(Adaptive)])
# # 
# # c(median(B_ERRORES[,c(Gellar)]), median(B_ERRORES[,c(New_approach)]))
# # 
# # #median(B_ERRORES[,c(Adaptive)])
# # 
# # c(sd(B_ERRORES[,c(Gellar)]), sd(B_ERRORES[,c(New_approach)]))
# #   
# # #sd(B_ERRORES[,c(Adaptive)])
# # #
# 
# #### MEANS TO AN EXCEL FILE
# 
# Y_means=colMeans(Y_ERRORES)
# B_means=colMeans(B_ERRORES)[1:24]
# b_matrix=matrix(data = B_means,nrow = 3,byrow = 1)
# 
# Y_matrix=matrix(data = Y_means,nrow = 4,byrow = 1)
# 
# B_matrix=b_matrix[-2,]
# 
# # write_xlsx(as.data.frame(Y_matrix),path = "C:/Users/user/Desktop/Trabajo/Escuela/Doctorado/Pavel/Tesis/Código/Mix models SOP/Propios/Simulaciones/Adaptive/Curvas más suaves/K-fold/R=100/1 y 1/Normal/T_i NegBin/N=100/Despues de IWSM/Y_means.xlsx")
# 
# # write_xlsx(as.data.frame(B_matrix),path = "C:/Users/user/Desktop/Trabajo/Escuela/Doctorado/Pavel/Tesis/Código/Mix models SOP/Propios/Simulaciones/Adaptive/Curvas más suaves/K-fold/R=100/1 y 1/Normal/T_i NegBin/N=100/Despues de IWSM/B_means.xlsx")
# 
# write_xlsx(as.data.frame(Y_matrix[4,]),path = "C:/Users/user/Desktop/Trabajo/Escuela/Doctorado/Pavel/Tesis/Código/Mix models SOP/Propios/Simulaciones/Adaptive/Curvas más suaves/K-fold/R=100/1 y 1/Poisson/T_i Uniform/N = 500/Y_means_no_vc.xlsx")
# 
# 
# ############ SD TO AN EXCEL FILE
# 
# Y_std=matrix(0,nrow = 3,ncol=32)
# B_std=matrix(0,nrow = 3,ncol=24)
# q_ric_B=q_ric_Y=NULL
# N_cases=c(100,200,500)
# outliers_B=outliers_Y=list()
# 
# for (i_ind in 1:3) {
#   
#   print(c("i_ind =",i_ind))
#   
#   if (i_ind%%3==1) {
#     
#     texto="N = 100/Data.RData"
#     texto_env="N = 100/Data_no_vc.RData"
#   }
#   if (i_ind%%3==2) {
#     
#     texto="N = 200/Data.RData"
#     texto_env="N = 200/Data_no_vc.RData"
#   }
#   if (i_ind%%3==0) {
#     
#     texto="N = 500/Data.RData"
#     texto_env="N = 500/Data_no_vc.RData"
#   }
#   
#   nam_X <- paste("C:/Users/user/Desktop/Trabajo/Escuela/Doctorado/Pavel/Tesis/Código/Mix models SOP/Propios/Simulaciones/Adaptive/Curvas más suaves/K-fold/R=100/1 y 1/Poisson/T_i Uniform/", texto ,sep="")
#   
#   load(nam_X)
#   
#   nam_env=paste("env","NegBin_POISSON_N",N_cases[i_ind],sep = "_")
#   
#   env_global=assign(nam_env,new.env())
#   
#   nam_load=paste("C:/Users/user/Desktop/Trabajo/Escuela/Doctorado/Pavel/Tesis/Código/Mix models SOP/Propios/Simulaciones/Adaptive/Curvas más suaves/K-fold/R=100/1 y 1/Poisson/T_i Uniform/", texto_env, sep="")
#   
#   load(nam_load, envir = get(nam_env))
#   
#   Y_ERRORES[,25:32]=env_global$Y_ERRORES[,25:32]
#   
#   for (iind in 1:(dim(Y_ERRORES)[2])) {
#     
#     nam_RIC_Y=paste("RIC_case_Y",iind,sep = "_")
#     
#     q_Y=quantile(Y_ERRORES[,iind])
#     
#     q_ric_Y[iind]=assign(nam_RIC_Y,q_Y[4]-q_Y[2]) #sqrt(diag(var(Y_ERRORES)))
#     
#     outliers_Y[[iind]]=which(Y_ERRORES[,iind]>q_Y[4]+1.5*q_ric_Y[iind] | Y_ERRORES[,iind]<q_Y[2]-1.5*q_ric_Y[iind])
#     
#     if (iind<25) {
#       
#       nam_RIC_B=paste("RIC_case_B",iind,sep = "_")
#       
#       q_B=quantile(B_ERRORES[,iind])
#       
#       q_ric_B[iind]=assign(nam_RIC_B,q_B[4]-q_B[2]) #sqrt(diag(var(Y_ERRORES)))
#       
#       outliers_B[[iind]]=which(B_ERRORES[,iind]>q_B[4]+1.5*q_ric_B[iind] | B_ERRORES[,iind]<q_B[2]-1.5*q_ric_B[iind])
#       
#     }
#     
#   }
#   
#   for (j_ind in 1:32) {
#     
#     print(c("j_ind =",j_ind))
#     
#     # case=j_ind
#     # 
#     # Gellar=case
#     # Gellar_te=1*8+case
#     # New_approach=2*8+case
#     # NO_VC=3*8+case
#     
#     if (length(as.double(outliers_Y[[j_ind]]))==0) {
#       
#       Y_std[i_ind,j_ind]=sqrt(var(Y_ERRORES[,j_ind]))
#       
#     }else{
#       
#       Y_std[i_ind,j_ind]=sqrt(var(Y_ERRORES[-outliers_Y[[j_ind]],j_ind]))
#     }
#     
#     if (j_ind<25) {
#       
#       if (length(as.double(outliers_B[[j_ind]]))==0) {
#         
#         B_std[i_ind,j_ind]=sqrt(var(B_ERRORES[,j_ind]))
#         
#       }else{
#         B_std[i_ind,j_ind]=sqrt(var(B_ERRORES[-outliers_B[[j_ind]],j_ind]))
#         
#       }
#       
#     }
#     
#   }
# }
# 
# # Y  Gellar Goldsmith SS-VDFR SB-VDFR
# # B  Gellar Goldsmith SS-VDFR
# 
# aux_1=matrix(Y_std[1,],nrow = 4, ncol = 8, byrow = 1)
# aux_2=matrix(Y_std[2,],nrow = 4, ncol = 8, byrow = 1)
# aux_3=matrix(Y_std[3,],nrow = 4, ncol = 8, byrow = 1)
# 
# Y_matrix=rbind(aux_1, aux_2, aux_3)
# 
# aux_1=matrix(B_std[1,],nrow = 3, ncol = 8, byrow = 1)
# aux_2=matrix(B_std[3,],nrow = 3, ncol = 8, byrow = 1)
# 
# B_matrix=rbind(aux_1, aux_2)
# 
# write_xlsx(as.data.frame(Y_matrix),path = "C:/Users/user/Desktop/Trabajo/Escuela/Doctorado/Pavel/Tesis/Código/Mix models SOP/Propios/Simulaciones/Adaptive/Curvas más suaves/K-fold/R=100/1 y 1/Poisson/T_i Uniform/Y_sd.xlsx")
# 
# write_xlsx(as.data.frame(B_matrix),path = "C:/Users/user/Desktop/Trabajo/Escuela/Doctorado/Pavel/Tesis/Código/Mix models SOP/Propios/Simulaciones/Adaptive/Curvas más suaves/K-fold/R=100/1 y 1/Poisson/T_i Uniform/B_sd.xlsx")
# 
# # COMPUTATION TIMES
# 
# colMeans(time_SOP)
# 
# colMeans(time_no_vc)
# 
# colMeans(time_Gellar)
# 
# colMeans(time_Goldsmith)
# 
