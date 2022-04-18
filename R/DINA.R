#' DINA
#'
#' This function implements DINA ---  an estimator of heterogeneous treatment effects for general responses. DINA makes use of existing machine learning tools and enjoys desirable statistical properties.
#'
#' @param data
#' @param family "gaussian"
#' @param z splines of the time; default indicator of treatment; TODO
#' @param prop "logistic regression"
#' @param baseline "gradient boosting"
#' @param k 2
#' @param crossFitting TRUE
#' @param bootstrap to bagging
#' @param B 100
#' @param inference FALSE
#' @param ... pass on to inference
#'
#' @export
DINA = function(data, family = "gaussian", prop = NULL, baseline = "gradient boosting", k = 2, crossFitting = TRUE, bootstrap = F, B = 1, full = F, weight.mechanism = "random"){
  n = length(data$Y); d = dim(data$X)[2]
  params = list(); paramsW = list()
  params$max_depth = 5; params$eta = 0.99
  if(family == "gaussian"){
    params$objective = "reg:squarederror"
    link.function = function(x){x}; var.function = function(x){x * 0 + 1}
  }
  if(family == "binomial"){
    params$objective = "binary:logistic"
    params$eval_metric = "logloss"
    link.function = function(x){log(x/(1-x))}; var.function = function(x){x*(1-x)}
  }
  if(family == "poisson"){
    params$objective = "count:poisson"
    link.function = function(x){log(x)}; var.function = function(x){x}
  }
  if(family == "cox"){
    params$objective = "survival:cox"
    # TODO!!!
  }

  if(!crossFitting){
    nuisance.function = DINAEstimator.step1(data$X, data$Y, data$W, family = family, prop = prop, params = params, paramsW = paramsW, link.function = link.function, var.function = var.function, baseline = baseline)
    estimator = DINAEstimator.step2(X = data$X, Y = data$Y, W = data$W, family = family, nuisance.function = nuisance.function, full = F)
    return (estimator)
  }

  foldSize = floor(n/k); record = list()
  bootstrap.membership = matrix(0, nrow = B, ncol = n)
  for(l in 1:B){
    foldIndex = list(); estimator = list()
    if(bootstrap){
      weight = rmultinom(n = 1, size = n, prob = rep(1/n,n))[,1]
      bootstrap.membership[l,] = weight
    } else {weight = rep(1, n)}
    index = sample(n, n, replace = F); index.resample = rep(index, weight[index])
    if(weight.mechanism == "random"){ #!!! to write as a helper function...
      start = 0; end = 0
      for(i in 1:(k-1)){
        end = end + sum(weight[index[seq(1+foldSize*(i-1), foldSize*i)]])
        foldIndex[[i]] = seq(start+1, end); start = end}
      foldIndex[[k]] = seq(start+1, n)}

    if(weight.mechanism == "cumsum"){
      start = 0; end = 0; temp = cumsum(weight[index])
      for(i in 1:k){
        end = temp[which.min(abs(temp - (n/k)*i))]
        foldIndex[[i]] = seq(start+1, end); start = end}}
    if(!((!bootstrap) * full)){
      for(i in 1:k){
        nuisance.function = DINAEstimator.step1(data$X[index.resample[-foldIndex[[i]]],], data$Y[index.resample[-foldIndex[[i]]]], data$W[index.resample[-foldIndex[[i]]]], family = family, prop = prop, params = params, paramsW = paramsW, link.function = link.function, var.function = var.function, baseline = baseline)
        estimator[[i]] = DINAEstimator.step2(X = data$X[index.resample[foldIndex[[i]]],], Y = data$Y[index.resample[foldIndex[[i]]]], W = data$W[index.resample[foldIndex[[i]]]], family = family, nuisance.function = nuisance.function)
      }
    } else {
      Sigma1 = Sigma2 = list()
      for(i in 1:k){
        nuisance.function = DINAEstimator.step1(data$X[index.resample[-foldIndex[[i]]],], data$Y[index.resample[-foldIndex[[i]]]], data$W[index.resample[-foldIndex[[i]]]], family = family, prop = prop, params = params, paramsW = paramsW, link.function = link.function, var.function = var.function, baseline = baseline)
        temp = DINAEstimator.step2(X = data$X[index.resample[foldIndex[[i]]],], Y = data$Y[index.resample[foldIndex[[i]]]], W = data$W[index.resample[foldIndex[[i]]]], family = family, nuisance.function = nuisance.function, full = T) # TODO: can merge the two? pass on full or call inference T
        estimator[[i]] = temp$estimator
        Sigma1[[i]] = temp$Sigma1
        Sigma2[[i]] = temp$Sigma2
      }
    }
    estimator.ave = apply(matrix(unlist(estimator), byrow = T, nrow = k), 2, mean, na.rm = T)
    names(estimator.ave) = names(estimator[[1]]);
    record[[l]] = estimator.ave
  }
  if(B==1){
    if(!full){return(record[[1]])}
    # inference using the asymptotic variance matrix
    return(result) # result$Sigma, Sigma1, Sigma2, estimator
  }
  if(!full){return(record)}
  result = list(); result$estimator = record; result$bootstrap.membership = bootstrap.membership; return (result)
}

DINAEstimator.step1 = function(X, Y, W, family, prop, baseline = "gradient boosting", params, paramsW = NA, link.function, var.function){
  if(is.character(prop)){
    if(prop == "logistic regression"){
      data.prop = data.frame(W, X); colnames(data.prop) = c("W", paste("X", seq(1,d), sep = ""))
      prop.model = glm(W~., data = data.prop, family = "binomial")
    } else if(prop == "gradient boosting"){#????
      prop.model = xgboost::xgboost(data = X, label = W, nrounds = 100, params = paramsW, verbose = 0)
    }
  }else if(is.function(prop)){prop.model = prop
  }else {stop("please input a method or a function for propensity score")}
  if(family == "gaussian"){
    if(is.character(baseline)){
      if(baseline == "gradient boosting"){
        m.model = xgboost::xgboost(data = X, label = Y, nrounds = 100, params = params, verbose = 0)
      }
      if(baseline == "regression"){
        data.baseline = data.frame(Y, X); colnames(data.baseline) = c("Y", paste("X", seq(1,d), sep = ""))
        m.model = glm(Y~., data = data.baseline, family = family)
      }
    }
    if(is.function(baseline)){m.model = baseline}
    nuisance.function = function(X.test, prop.model.test = prop.model, m.model.test = m.model, prop.test = prop, baseline.test = baseline){
      data.test = data.frame(X.test); colnames(data.test) = c(paste("X", seq(1,d), sep = ""))
      if(is.character(prop.test)){
        if(prop.test == "logistic regression"){
          a = predict(prop.model.test, newdata = data.test, type = "response")
        } else if(prop.test == "gradient boosting"){#!!!!
          a = predict(prop.model.test, newdata = X.test)}
      }
      if(is.function(prop.test)){a = prop.test(X.test)}
      if(is.character(baseline.test)){
        if(baseline.test == "gradient boosting"){nu = predict(m.model.test, newdata = X.test)}
        if(baseline.test == "regression"){nu = predict(m.model.test, newdata = data.test, type = "response")}
      }
      if(is.function(baseline.test)){nu = m.model.test(X.test)$nu}
      result = list(); result$a = a; result$nu = nu; return(result)
    }
  } else {
    if(is.character(baseline)){
      if(baseline == "gradient boosting"){
        trt.model = xgboost::xgboost(data = X[W==1, , drop = F], label = Y[W==1], nrounds = 100, params = params, verbose = 0)
        cnt.model = xgboost::xgboost(data = X[W==0, , drop = F], label = Y[W==0], nrounds = 100, params = params, verbose = 0)
      }
      if(baseline == "regression"){
        data.baseline = data.frame(Y, X); colnames(data.baseline) = c("Y", paste("X", seq(1,d), sep = ""))
        trt.model = glm(Y~., data = data.baseline[W==1,], family = family)
        cnt.model = glm(Y~., data = data.baseline[W==0,], family = family)}}
    if(is.function(baseline)){trt.model = NULL; cnt.model = NULL}

    nuisance.function = function(X.test, prop.model.test = prop.model, trt.model.test = trt.model, cnt.model.test = cnt.model, prop.test = prop, baseline.test = baseline){
      data.test = data.frame(X.test); colnames(data.test) = c(paste("X", seq(1,d), sep = ""))
      if(is.character(prop.test)){
        if(prop.test == "logistic regression"){a = predict(prop.model.test, newdata = data.test, type = "response")
        } else if(prop.test == "gradient boosting"){
          a = predict(prop.model.test, newdata = X.test) #!!!!
        }
      }
      if(is.function(prop.test)){a = prop.model.test(X.test)}
      if(is.character(baseline.test)){
        if(baseline.test == "gradient boosting"){
          mu1 = predict(trt.model.test, newdata = X.test)
          mu0 = predict(cnt.model.test, newdata = X.test)
        }
        if(baseline.test == "regression"){
          mu1 = predict(trt.model.test, newdata = data.test, type = "response")
          mu0 = predict(cnt.model.test, newdata = data.test, type = "response")
        }
        a = (a * var.function(mu1))/(1e-6 + (1-a) * var.function(mu0)); a = a/(1+a)
        nu = a * link.function(mu1) + (1-a) * link.function(mu0)
      }
      if(is.function(baseline.test)){
        temp = baseline.model.test(X.test)
        a = (a * var.function(temp$mu1))/(1e-6 + (1-a) * var.function(temp$mu0)); a = a/(1+a)
        nu = temp$nu
      }
      result = list(); result$a = a; result$nu = nu; return(result)
    }
  }
  return (nuisance.function)
}


# debug
# test3 = DINAEstimator.step2(X = data$X, Y = data$Y, W = data$W, family = family, nuisance.function = test)
DINAEstimator.step2 = function(X, Y, W, family, nuisance.function, full = F){
  temp = nuisance.function(X)
  data.HTE = data.frame(Y, temp$nu, W - temp$a, I(X)); colnames(data.HTE) = c("Y", "nu", "resW", "X")
  glm.model = tryCatch(glm(Y ~ X:resW + resW + offset(nu) - 1, data = data.HTE, family = family), error=function(e) e, warning=function(w) w)
  if(is(glm.model,"warning")){estimator = rep(NA, d+1)
  } else {estimator = coef(glm.model)}
  if(!full){return (estimator)}
  if(family != "gaussian"){stop("current sd estimates only support family = gaussian")}
  result = list(); result$estimator = estimator
  return(result) # $Sigma1, $Sigma2
}

