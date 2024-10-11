library(refund)
library(fda)
library(mgcv)
library(SOP)
library(plotly)
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

k=10  # NUMBER OF GROUPS IN THE K-fold

fold=N/k

R=10 # NUMBER OF ITERATIONS FOR THE SIMULATION STUDY

c1=30
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
y_h_sop_all=y_h_ellos_all=array(dim=c(R,N,case))
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

err_y_ellos_all=err_y_sop_all=rep(0,R)
#

start=proc.time()

for (iter_out in 3) { # HERE WE CAN SELECTED WHICH SCENARIO(S) SIMULATE
  
  print(c("case = ",iter_out))
  
  for (iter in 1:R) {
    
    print(c("iter = ",iter))
    
    set.seed(1000+iter)
    
    M = round(runif(N,10,J),digits = 0) # HERE WE GENERATE THE DOMAIN FOR ALL SUBJECTS WITH A MINIMUM OF A 10 OBSERVATIONS
    
    # M=sort(rep((J-9):J,10)) #rep((1:5)*10,20)
    
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
    
    # t=seq(from = 0, to = 1, length.out=M_max)
    
    t=1:M_max
    
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
        
        temp[e,1:M[i]]=v_i1*sin(2*pi*e*(1:M[i])/J)+v_i2*cos(2*pi*e*(1:M[i])/J)
      }
      
      B=apply(temp,2,sum)
      
      B=B+u
      
      X_s[i,]=B
      
      aux=var(B,na.rm=TRUE)
      
      X_se[i,]=(B)+rnorm(M_max,0,aux/4) # WE ADD NOISE
      
    }
    
    # j=10
    # plot(X_s[j,])
    # lines(X_se[j,],col=3,type="o")

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
      
      nu[i]=sum(X_s[i,]*Beta[i,,iter_out],na.rm = 1)/(M[i]) #NOT NOISY
      
      # if (iter_out==8) {
      #   nu[i]=sum(X_se[i,]*Beta[i,,4],na.rm = 1)/(M[i]) # NOISY
      # }
      # if (iter_out<=4) {
      #   nu[i]=sum(X_s[i,]*Beta[i,,iter_out],na.rm = 1)/(M[i]) #NOT NOISY
      # }
      # if(iter_out>4 & iter_out<8){
      #   nu[i]=sum(X_se[i,]*Beta[i,,iter_out%%4],na.rm = 1)/M[i] # NOISY
      # }
      
    }
    
    # plot_ly(z=Beta[,,2], type="surface")
    
    # var_e <- (1 / Rsq - 1) * stats::var(nu) # (1-Rsq)*var(nu[ind,])
    # 
    # y=nu+rnorm(N,sd = sqrt(var_e)) # ADDING NOISE TO THE GAUSSIAN MODEL
    
    y=rpois(N,exp(nu)+1) # POISSON MODEL ##### CAHNGE THE ERROR EQUATIONS IN LINE 477
    
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
    
    groups=kfold(data,k=k) #rep(1:k,fold) # # creating the K-fold groups
    

    
    for (group in 1:k) {
      
      test=which(groups==group)
      train=which(groups!=group)
      
      if (M[tail(train,1)]!=M_max) {

        aux=test[length(test)]
        test[length(test)]=train[length(train)]
        train[length(train)]=aux
        test=sort(test)
      }
      
      if (M[head(train,1)]!=min(M)) {
        
        aux=test[1]
        test[1]=train[1]
        train[1]=aux
        test=sort(test)
      }
      
      print(c("group = ",group))
      
      X_test=X_se[test,]
      X_train=X_se[train,]
      
      X_reg_train=X_reg[train,]
      X_reg_test=X_reg[test,]
      
      M_test=M[test]
      M_train=M[train]
      
      M_it_test[iter,,iter_out]=M_test
      M_it_train[iter,,iter_out]=M_train
      
      y_test=y[test]
      y_train=y[train]
      
      # MODEL ESTIMATION BY OUR APPROACH

      data_train=data.frame(y_train)
      data_train[["X_train"]]=X_train
      
      data_test=data.frame(y_test)
      data_test[["X_test"]]=X_test
      
      start_SOP=proc.time()
      # formula <- y ~ ffvd(X_se, t, nbasis = c(c1,c2,c3))
      # BB=ffvd(X_se, t[1:max(M)], nbasis = c(c1,c2,c3))
      # res <- VDPO(formula = formula, data = data)#,family = poisson())
      
      formula <- y_train ~ ffvd(X_train, t, nbasis = c(c1,c2,c3))
      BB=ffvd(X_train, t[1:max(M_train)], nbasis = c(c1,c2,c3))
      res <- VDPO(formula = formula, data = data_train,family = poisson())
      
      end_SOP=proc.time()
      time_SOP[iter,iter_out]=end_SOP[1]-start_SOP[1]
      #
      
      # MODEL ESTIMATION USING GELLAR APPROACH
      
      start_Gellar=proc.time()
      
      fit_Gellar <- pfr(y_train ~ lf.vd(X_train,k=89),family = poisson())
      
      end_Gellar=proc.time()
      time_Gellar[iter,iter_out]=end_Gellar[1]-start_Gellar[1]
      #
      
      # MODEL ESTIMATION USING GOLDSMITH APPROACH
      
      start_Goldsmith=proc.time()
      
      fit_Gellar_te <- pfr(y_train ~ lf(X_reg_train, bs = "ps",k=15,presmooth = "bspline",presmooth.opts = list(nbasis=15)),family = poisson())
      
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
      
      BB_test=ffvd(X_test, t[1:max(M_test)], nbasis = c(c1,c2,c3))
      # 
      # as.matrix(kronecker(BB_test$L_Phi[[1]], t(BB_test$B_T[[1]][j,]))) %*% res$theta_ffvd
      # 
      # y_h_sop[iter,,iter_out] = BB_test$B_ffvd %*% res$theta_ffvd # ESTIMATED REPSONSE VARIABLE USING OUR APPROACH
      
      
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
        #
        
        # Beta_sop=as.matrix(kronecker(BB_test$L_Phi[[1]][1:M_it_test[iter,j,iter_out],], t(BB_test$B_T[[1]][j,]))) %*% res$theta_ffvd
        
        # Beta_sop_2=as.matrix(kronecker(BB$L_Phi[[1]][1:M_it_test[iter,j,iter_out],], t(BB$B_T[[1]][test[j],]))) %*% res$theta_ffvd
        
        
        B_T_fold=splines::spline.des(BB$B_T$knots, M_test[j], 3 + 1, 0 * M_test[j])$design
        L_Phi_fold=splines::spline.des(BB$L_Phi$knots, t[1:M_test[j]], 3 + 1, 0 * t[1:M_test[j]])$design
        
        Beta_sop=as.matrix(kronecker(L_Phi_fold,B_T_fold)) %*% res$theta_ffvd
        
        y_h_sop[iter,j,iter_out]=sum(BB_test$X_hat[j,1:M_it_test[iter,j,iter_out]]*Beta_sop,na.rm = 1)/M_it_test[iter,j,iter_out]
        
      }
      
      # j=10
      # plot(Beta[which(groups==group)[j],,2])
      # lines(Beta_sop,col=2)
      # lines(Beta_refund)
      # 
      # 
      # plot(X_s[test[j],])
      # lines(X_test[j,1:M_it_test[iter,j,iter_out]],col=3,type="o")
      # lines(BB_test$X_hat[j,1:M_it_test[iter,j,iter_out]],col=2,type="o")
      
      # ESTIMATION ERRORS FOR THE TEST DATA SET IN THE CASE OF NORMAL RESPONSE
      
      # error_group_VDFR[iter,group,iter_out]=sqrt(sum((y_test-y_h_ellos[iter,,iter_out])^2)/fold)
      # error_group_FF_VDFR[iter,group,iter_out]=sqrt(sum((y_test-y_h_sop[iter,,iter_out])^2)/fold)
      # error_group_Carmen[iter,group,iter_out]=sqrt(sum((y_test-y_h_ellos_te[iter,,iter_out])^2)/fold)
      
      # ESTIMATION ERRORS FOR THE TEST DATA SET IN THE CASE OF POISSON RESPONSE
      
      error_group_VDFR[iter,group,iter_out]=sqrt(sum((y_test-(exp(y_h_ellos[iter,,iter_out])))^2)/fold)
      error_group_Carmen[iter,group,iter_out]=sqrt(sum((y_test-(exp(y_h_ellos_te[iter,,iter_out])))^2)/fold)
      error_group_FF_VDFR[iter,group,iter_out]=sqrt(sum((y_test-(exp(y_h_sop[iter,,iter_out])))^2)/fold)
      
      
    } # HERE END THE FOR OF group = 1:fold
    
    Y_ERROR_2_ellos[iter,iter_out]= mean(error_group_VDFR[iter,,iter_out])
    Y_ERROR_2_te_ellos[iter,iter_out]= mean(error_group_Carmen[iter,,iter_out])
    Y_ERROR_2_sop[iter,iter_out]= mean(error_group_FF_VDFR[iter,,iter_out])
    
    # HERE WE ESTIMATE THE FUNCTIONAL COEFFICIENT 
    
    formula <- y ~ ffvd(X_se, t, nbasis = c(c1,c2,c3))
    BB_all=ffvd(X_se, t, nbasis = c(c1,c2,c3))
    res_Beta <- VDPO(formula = formula, data = data,family = poisson())

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
    Beta_FFVD[,1:M_max,iter,iter_out]=res_Beta$Beta_ffvd[[1]]
    
    error_Beta_FFVD[iter,iter_out]=sum(((Beta[,,iter_out]-Beta_FFVD[,1:M_max,iter,iter_out])^2)/(J*(J+1)), na.rm=TRUE)
    
    Beta_FFVD_inf=Beta_FFVD_sup=matrix(nrow = N,ncol = M_max)
    
    aux_base=ffvd(X_se, t, nbasis = c(c1,c2,c3))
    
    for (aux_ind in 1:N) {
      
      prod=as.matrix(kronecker(aux_base$L_Phi$B[1:M[aux_ind],],t(aux_base$B_T$B[aux_ind,])))
      
      var_curve = diag(prod %*% res_Beta$covar_theta %*% t(prod))
      
      std_curve=sqrt(var_curve)
      
      Beta_FFVD_inf[aux_ind,1:M[aux_ind]] =  Beta_FFVD[aux_ind,1:M[aux_ind],iter,iter_out] - 1.96*std_curve
      Beta_FFVD_sup[aux_ind,1:M[aux_ind]] =  Beta_FFVD[aux_ind,1:M[aux_ind],iter,iter_out] + 1.96*std_curve
      
    }
    
    Beta_FFVD_inf_ind[[iter]]=Beta_FFVD_inf
    Beta_FFVD_sup_ind[[iter]]=Beta_FFVD_sup
    
    #####
    
    fit_Gellar_all <- pfr(y ~ lf.vd(X_se,k=89),family = poisson())
    
    Beta_G[[iter]]=coef(fit_Gellar_all, n=c(length(t),length(unique(M))))$value
    
    for (j in 1:M_max) {
      
      ind=which(unique(M_it[iter,,iter_out])==M[j])
      ind_t=max(M_it[iter,,iter_out]) # OJO AQU? ESTAMOS SUPONIENDO QUE LAS OBSERVACIONES SON CONSECUTIVAS DESDE 1:MAX(M)
      
      Beta_refund=Beta_G[[iter]][(ind_t*(ind-1)+1):(ind*ind_t)]
      Beta_VD[j,1:M[j],iter,iter_out]=Beta_refund[1:M_it[iter,j,iter_out]]
      
    }
    
    error_Beta_VD[iter,iter_out]=sum(((Beta[,,iter_out]-Beta_VD[,1:M_max,iter,iter_out])^2)/(J*(J+1)), na.rm=TRUE)
    
    for (j in 1:N) {
      
      y_h_ellos_all[iter,j,iter_out] = sum(X_se[j,1:M[j]]*Beta_VD[j,1:M[j],iter,iter_out],na.rm = 1)/M[j]        
      y_h_sop_all[iter,j,iter_out] = sum(BB_all$X_hat[j,1:M[j]]*Beta_FFVD[j,1:M[j],iter,iter_out],na.rm = 1)/M[j]      
      
      }
    
    err_y_sop_all[iter]=sqrt(sum((y-exp(y_h_sop_all[iter,,iter_out]))^2)/N)
    err_y_ellos_all[iter]=sqrt(sum((y-exp(y_h_ellos_all[iter,,iter_out]))^2)/N)
      

    
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

case_iter=3

plot_ly(z=Beta_VD[,,case_iter,iter_out], type="surface")
plot_ly(z=Beta_FFVD[,,case_iter,iter_out], type="surface")
plot_ly(z=Beta[,,iter_out], type="surface")

# save.image("C:/Users/user/Desktop/Trabajo/Escuela/Doctorado/Pavel/Tesis/1er Paper/Codigo/Simulation-1st-Paper/Results new method/vacations.RData")

# aux_save=data.frame(Y_ERROR_2_sop[,1],Y_ERROR_2_ellos[,1],Y_ERROR_2_te_ellos[,1],error_Beta_FFVD[1:R,1],error_Beta_VD[1:R,1])

# new_method=data.frame(Y_ERROR_2_sop[,1],Y_ERROR_2_ellos[,1],Y_ERROR_2_te_ellos[,1],error_Beta_FFVD[1:R,1],error_Beta_VD[1:R,1])

B_ERRORES=data.frame(error_Beta_FFVD,error_Beta_VD) #B_ERROR_2,B_ERROR_2_te,

Y_ERRORES=data.frame(Y_ERROR_2_ellos,Y_ERROR_2_te_ellos,Y_ERROR_2_sop) # #Y_ERROR_2,Y_ERROR_2_te,

# write_xlsx(B_ERRORES,path = "C:/Users/user/Desktop/Trabajo/Escuela/Doctorado/Pavel/Tesis/Código/Mix models SOP/Propios/Simulaciones K-fold good version N_200/Poisson_Uniform_N_200/B_ERRORES.xlsx")
# 
# write_xlsx(Y_ERRORES,path = "C:/Users/user/Desktop/Trabajo/Escuela/Doctorado/Pavel/Tesis/Código/Mix models SOP/Propios/Simulaciones K-fold good version N_200/Poisson_Uniform_N_200/Y_ERRORES.xlsx")

Y_values=cbind(c(Y_ERRORES[,3],Y_ERRORES[,11],Y_ERRORES[,19]))
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

B_values=cbind(c(B_ERRORES[,3],B_ERRORES[,11]))
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


#### Y ##################

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

case=100
plot(X_s[case,])
lines(X_se[case,],type="o",col=3)
lines(BB_all$X_hat[case,],type="o",col=2)

case=1
plot(X_s[test[case],])
lines(X_test[case,],type="o",col=3)
lines(BB_test$X_hat[case,],col=2,type="o")


j=1
case=which(groups==group)[j]
plot(Beta[which(groups==group)[j],,2])
lines(Beta_sop,col=2)
lines(Beta_refund)

case=10
plot(Beta[case,,2])
lines(Beta_VD[case,1:M[case],iter,iter_out])
lines(Beta_FFVD[case,1:M[case],iter,iter_out],col=2,lwd=2)
