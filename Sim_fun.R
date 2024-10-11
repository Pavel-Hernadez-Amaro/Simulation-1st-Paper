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
library(remotes)
# install_github("Pavel-Hernadez-Amaro/VDPO")
library(VDPO)

### DEFINE CONSTANSTS ###
Rsq <- 0.95

N <- 100 # NUMBER OF SUBJECTS

J <- 100 # NUMBER OF MAXIMUM OBSERVATIONS PER SUBJECTS

k <- 10 # NUMBER OF GROUPS IN THE K-fold

fold <- N / k

R <- 10 # NUMBER OF ITERATIONS FOR THE SIMULATION STUDY

c1 <- 10
c2 <- 10
c3 <- 10

registrated_points <- function(x) {
  J.i <- sum(!is.na(x))
  y <- approx(
    x = seq(0, 1, length = J.i), y = x[1:J.i],
    xout = seq(0, 1, length = J)
  )$y
  y
}

function_par <- function(iter,
                         beta_index = 1,
                         M_distribution = c("unif", "negbin"),
                         family = c("gaussian", "poisson")) {
  if (beta_index < 1 | beta_index > 4) {
    stop("'beta_index' must be inside the interval [1,4]", call. = FALSE)
  }
  
  M_distribution <- match.arg(M_distribution)
  family <- match.arg(family)
  
  print(c("iter = ", iter))
  
  set.seed(1000 + iter)
  
  
  if (M_distribution == "unif") {
    M <- round(runif(N, 10, J), digits = 0) # HERE WE GENERATE THE DOMAIN FOR ALL SUBJECTS WITH A MINIMUM OF A 10 OBSERVATIONS
  } else if (M_distribution == "negbin") {
    M <- rnegbin(N, 24, 1) # Para 1000 poner (240,2)
  }
  
  if (max(M) > J) {
    M[which(M > J)] <- J
  }
  
  
  if (min(M) <= 10) {
    M[which(M <= 10)] <- 10
  }
  
  M <- sort(M) # WE SORT THE DATA WITHOUT LOSS OF GENERALITY
  
  M_max <- max(M)
  
  # t=seq(from = 0, to = 1, length.out=M_max)
  
  t <- 1:M_max
  Beta_VD <- Beta_FFVD <- array(dim = c(N, J))
  ############ HERE WE GENERATE THE FUNCTIONAL DATA
  
  X_se <- matrix(NA, N, M_max) # NOISY
  
  X_s <- matrix(NA, N, M_max) # NOT NOISY
  
  for (i in 1:N) {
    u <- rnorm(1)
    
    temp <- matrix(NA, 10, M_max)
    
    for (e in 1:10) {
      v_i1 <- rnorm(1, 0, 4 / e^2)
      v_i2 <- rnorm(1, 0, 4 / e^2)
      
      temp[e, 1:M[i]] <- v_i1 * sin(2 * pi * e * (1:M[i]) / J) + v_i2 * cos(2 * pi * e * (1:M[i]) / J)
    }
    
    B <- apply(temp, 2, sum)
    
    B <- B + u
    
    X_s[i, ] <- B
    
    aux <- var(B, na.rm = TRUE)
    
    X_se[i, ] <- (B) + rnorm(M_max, 0, aux / 4) # WE ADD NOISE
  }
  
  X_reg <- t(apply(X_se, 1, registrated_points)) # THESE ARE THE REGISTRATED POINTS
  
  # SOME SAVE CHECKS FOR UNWANTED NAs
  for (i in 1:M_max) {
    if (length(which(is.na(X_s[, i]))) == N) {
      print(c("iteracion", l, "columna", i))
      stop("Hay columnas con NA")
    }
  }
  
  
  
  
  ###### HERE WE GENERATE THE TRUE Beta COEFFICIENT AND THE RESPONSE VARIABLE
  
  Beta <- array(dim = c(N, M_max))
  nu <- y <- rep(0, N)
  
  for (i in 1:N) {
    # TRUE FUNCTIONAL COEFFICIENTS
    
    if (beta_index == 1) {
      Beta[i, 1:(M[i])] <- ((10 * t[1:(M[i])] / M[i]) - 5) / 10
    } else if (beta_index == 2) {
      Beta[i, 1:(M[i])] <- ((1 - (2 * M[i] / M_max)) * (5 - 40 * ((t[1:(M[i])] / M[i]) - 0.5)^2)) / 10
    } else if (beta_index == 3) {
      Beta[i, 1:(M[i])] <- (5 - 10 * ((M[i] - t[1:(M[i])]) / M_max)) / 10
    } else if (beta_index == 4) {
      Beta[i, 1:(M[i])] <- (sin(2 * pi * M[i] / M_max) * (5 - 10 * ((M[i] - t[1:(M[i])]) / M_max))) / 10
    }
    
    
    # HERE WE GENERATE THE RESPONSE VARIABLE FOR EVERY BETA FROM THE 2 DIFFERENT FUNCTIONAL DATA (NOISY AND OTHERWISE)
    
    nu[i] <- sum(X_s[i, ] * Beta[i, ], na.rm = 1) / (M[i]) # NOT NOISY
  }
  
  
  if (family == "gaussian") {
    var_e <- (1 / Rsq - 1) * stats::var(nu)
    y <- nu + rnorm(N, sd = sqrt(var_e))
  } else if (family == "poisson") {
    y <- rpois(N, exp(nu))
  }
  
  data <- data.frame(y = y)
  data[["X_se"]] <- X_se
  data[["X_s"]] <- X_s
  data[["Beta"]] <- Beta
  
  ############## HERE ENDS THE GENERATION OF THE DATA
  
  ############ HERE BEGINS THE ESTIMATION OF THE MODEL
  
  # ASSIGNING THE TRAIN AND TEST SETS
  
  groups <- kfold(data, k = k) # creating the K-fold groups
  
  error_group_VDFR <- rep(0, k)
  error_group_FF_VDFR <- rep(0, k)
  error_group_Carmen <- rep(0, k)
  
  for (group in 1:k) {
    test <- which(groups == group)
    train <- which(groups != group)
    
    if (M[tail(train, 1)] != M_max) {
      aux <- test[length(test)]
      test[length(test)] <- train[length(train)]
      train[length(train)] <- aux
      test <- sort(test)
    }
    
    if (M[head(train, 1)] != min(M)) {
      aux <- test[1]
      test[1] <- train[1]
      train[1] <- aux
      test <- sort(test)
    }
    
    print(c("group = ", group))
    
    X_test <- X_se[test, ]
    X_train <- X_se[train, ]
    
    X_reg_train <- X_reg[train, ]
    X_reg_test <- X_reg[test, ]
    
    M_test <- M[test]
    M_train <- M[train]
    
    y_test <- y[test]
    y_train <- y[train]
    
    # MODEL ESTIMATION BY OUR APPROACH
    
    data_train <- data.frame(y_train)
    data_train[["X_train"]] <- X_train
    
    data_test <- data.frame(y_test)
    data_test[["X_test"]] <- X_test
    
    formula <- y_train ~ ffvd(X_train, t, nbasis = c(c1, c2, c3))
    BB <- ffvd(X_train, t[1:max(M_train)], nbasis = c(c1, c2, c3))
    
    
    if (family == "gaussian") {
      res <- vd_fit(formula = formula, data = data_train, family = stats::gaussian())
      fit_Gellar <- pfr(y_train ~ lf.vd(X_train, k = 89), family = stats::gaussian())
      fit_Gellar_te <- pfr(
        y_train ~ lf(X_reg_train, bs = "ps", k = 15, presmooth = "bspline", presmooth.opts = list(nbasis = 15)),
        family = stats::gaussian()
      )
      
    } else if (family == "poisson") {
      res <- vd_fit(formula = formula, data = data_train, family = stats::poisson())
      fit_Gellar <- pfr(y_train ~ lf.vd(X_train, k = 89), family = stats::poisson())
      fit_Gellar_te <- pfr(
        y_train ~ lf(X_reg_train, bs = "ps", k = 15, presmooth = "bspline", presmooth.opts = list(nbasis = 15)),
        family = stats::poisson()
      )
    }
    
    # MODEL ESTIMATION USING GELLAR APPROACH
    
    
    #
    
    # MODEL ESTIMATION USING GOLDSMITH APPROACH
    
    #################### HERE ENDS THE ESTIMATION OF THE MODELS
    
    #################### HERE WE CALCULATE THE ESTIMATION ERRORS
    
    error_2_ellos <- error_2_te_ellos <- matrix(0, nrow = (k - 1) * fold, ncol = 1)
    error_2_sop <- matrix(0, nrow = (k - 1) * fold, ncol = 1)
    
    error_2_ellos_test <- error_2_te_ellos_test <- matrix(0, nrow = fold, ncol = 1)
    error_2_sop_test <- matrix(0, nrow = fold, ncol = 1)
    
    ERROR_2_sop <- 0
    ERROR_2_ellos <- 0
    
    ############# HERE WE CALCULATE THE ERROR OF THE ESTIMATED RESPONSE VARIABLE
    
    BB_test <- ffvd(X_test, t[1:max(M_test)], nbasis = c(c1, c2, c3))
    
    
    # ESTIMATED REPSONSE VARIABLE USING HE GELLAR AND GOLDSMITH APPROACHES
    
    Beta_G <- coef(fit_Gellar, n = c(length(t), length(unique(M))))$value
    Beta_G_te <- coef(fit_Gellar_te, n = M_max)$value
    
    y_h_ellos <- rep(0, fold)
    y_h_ellos_te <- rep(0, fold)
    y_h_sop <- rep(0, fold)
    
    for (j in 1:fold) {
      ind <- which(unique(M) == M_test[j])
      ind_t <- max(M) # OJO AQUI ESTAMOS SUPONIENDO QUE LAS OBSERVACIONES SON CONSECUTIVAS DESDE 1:MAX(M)
      
      Beta_refund <- Beta_G[(ind_t * (ind - 1) + 1):(ind * ind_t)]
      Beta_refund <- Beta_refund[1:M_test[j]]
      
      Beta_refund_te <- Beta_G_te[1:M_test[j]]
      
      y_h_ellos[j] <- sum(X_test[j, 1:M_test[j]] * Beta_refund, na.rm = 1) / M_test[j]
      
      y_h_ellos_te[j] <- sum(X_test[j, 1:M_test[j]] * Beta_refund_te, na.rm = 1) / M_test[j]
      
      B_T_fold <- splines::spline.des(BB$B_T$knots, M_test[j], 3 + 1, 0 * M_test[j])$design
      L_Phi_fold <- splines::spline.des(BB$L_Phi$knots, t[1:M_test[j]], 3 + 1, 0 * t[1:M_test[j]])$design
      
      Beta_sop <- as.matrix(kronecker(L_Phi_fold, B_T_fold)) %*% res$theta_ffvd
      
      y_h_sop[j] <- sum(BB_test$X_hat[j, 1:M_test[j]] * Beta_sop, na.rm = 1) / M_test[j]
    }
    
    
    if (family == "gaussian") {
      # ESTIMATION ERRORS FOR THE TEST DATA SET IN THE CASE OF NORMAL RESPONSE
      
      error_group_VDFR[group] <- sqrt(sum((y_test - y_h_ellos)^2) / fold)
      error_group_FF_VDFR[group] <- sqrt(sum((y_test - y_h_sop)^2) / fold)
      error_group_Carmen[group] <- sqrt(sum((y_test - y_h_ellos_te)^2) / fold)
    } else if (family == "poisson") {
      # ESTIMATION ERRORS FOR THE TEST DATA SET IN THE CASE OF POISSON RESPONSE
      
      error_group_VDFR[iter,group,iter_out]=sqrt(sum((y_test-(exp(y_h_ellos[iter,,iter_out])))^2)/fold)
      error_group_Carmen[iter,group,iter_out]=sqrt(sum((y_test-(exp(y_h_ellos_te[iter,,iter_out])))^2)/fold)
      error_group_FF_VDFR[iter,group,iter_out]=sqrt(sum((y_test-(exp(y_h_sop[iter,,iter_out])))^2)/fold)
    }
    
    
    
    
  } # HERE END THE FOR OF group = 1:fold
  
  Y_ERROR_2_ellos <- mean(error_group_VDFR)
  Y_ERROR_2_te_ellos <- mean(error_group_Carmen)
  Y_ERROR_2_sop <- mean(error_group_FF_VDFR)
  
  # HERE WE ESTIMATE THE FUNCTIONAL COEFFICIENT
  
  formula <- y ~ ffvd(X_se, t, nbasis = c(c1, c2, c3))
  BB_all <- ffvd(X_se, t, nbasis = c(c1, c2, c3))
  
  if (family == "gaussian") {
    res_Beta <- vd_fit(formula = formula, data = data ,family = stats::gaussian())
    fit_Gellar_all <- pfr(y ~ lf.vd(X_se, k = 89), family = stats::gaussian())
  } else if (family == "poisson") {
    res_Beta <- vd_fit(formula = formula, data = data ,family = stats::poisson())
    fit_Gellar_all <- pfr(y ~ lf.vd(X_se, k = 89), family = stats::poison())
  }
  
  Beta_FFVD[, 1:M_max] <- res_Beta$Beta_ffvd[[1]]
  
  error_Beta_FFVD <- sum(((Beta - Beta_FFVD[, 1:M_max])^2) / (J * (J + 1)), na.rm = TRUE)
  
  #####
  Beta_G <- coef(fit_Gellar_all, n = c(length(t), length(unique(M))))$value
  
  for (j in 1:M_max) {
    ind <- which(unique(M) == M[j])
    ind_t <- max(M) # OJO AQUI ESTAMOS SUPONIENDO QUE LAS OBSERVACIONES SON CONSECUTIVAS DESDE 1:MAX(M)
    
    Beta_refund <- Beta_G[(ind_t * (ind - 1) + 1):(ind * ind_t)]
    Beta_VD[j, 1:M[j]] <- Beta_refund[1:M[j]]
  }
  
  error_Beta_VD <- sum(((Beta - Beta_VD[, 1:M_max])^2) / (J * (J + 1)), na.rm = TRUE)
  
  list(
    Y_ERROR_2_ellos = Y_ERROR_2_ellos,
    Y_ERROR_2_te_ellos = Y_ERROR_2_te_ellos,
    Y_ERROR_2_sop = Y_ERROR_2_sop,
    error_Beta_VD = error_Beta_VD,
    error_Beta_FFVD = error_Beta_FFVD,
    res_Beta = res_Beta
  )
} # HERE ENDS THE MIDDLE FOR ITERATION (RUNS ON R)