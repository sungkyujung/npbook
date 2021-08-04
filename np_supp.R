# nonparametric book by S. Jung
# Supplementary functions
# 
# ----------------


# One-way ANOVA F statistic -------------------------------------------------------------

F.stat = function(x, grp, data = NULL){
  pars <- as.list(match.call()[-1])
  if(!is.null(data)){
    x = eval(pars$data)[,as.character(pars$x)]
    grp = eval(pars$data)[,as.character(pars$grp)]
  }
  groupmeans_ss = aggregate(x = x, 
                            by = list(grp), 
                            FUN = function(x) {c(mean(x),length(x))})
  n = length(x)
  k = nrow(groupmeans_ss)
  x.tt = mean(x)
  SST = sum( (x - x.tt)^2 )
  SSTrt = sum(groupmeans_ss$x[,2]*(groupmeans_ss$x[,1] - x.tt)^2) 
  SSE = SST - SSTrt
  return(  SSTrt / (k-1) / SSE * (n-k) )
} 


# One-way ANOVA Kruskal-Wallis statistic ------------------------------------------------

KW.stat = function(x, grp, data = NULL){
  pars <- as.list(match.call()[-1])
  if(!is.null(data)){
    x = eval(pars$data)[,as.character(pars$x)]
    grp = eval(pars$data)[,as.character(pars$grp)]
  }
  x = rank(x)
  groupmeans_ss = aggregate(x = x, 
                            by = list(grp), 
                            FUN = function(x) {c(mean(x),length(x))})
  n = length(x)
  k = nrow(groupmeans_ss)
  SSTrt = sum(groupmeans_ss$x[,2]*(groupmeans_ss$x[,1] - (n+1)/2)^2) 
  return( SSTrt / var(x) )
} 


# Rank tranformation within blocks ----------------------------------------

rankinblock = function(x,block){
  require(dplyr)
  out = data.frame(x,block) %>% group_by(block) %>% mutate(ranks = rank(x))
  return(out$ranks)
}


# Complete Block Two Way ANOVA F statistic --------------------------------

F.stat.CBD = function(x, trt, blk, data = NULL){
  pars <- as.list(match.call()[-1])
  if(!is.null(data)){
    x = eval(pars$data)[,as.character(pars$x)]
    trt = eval(pars$data)[,as.character(pars$trt)]
    blk = eval(pars$data)[,as.character(pars$blk)]
  }
  
  x.tt = mean(x)
  x.trt = aggregate(x = x, by = list(trt), FUN = mean)
  x.blk = aggregate(x = x, by = list(blk), FUN = mean)
  n = length(x)
  k = nrow(x.trt)
  b = nrow(x.blk)
  
  SStrt = b * sum(   (x.trt$x - x.tt)^2)
  SSblk = k * sum(   (x.blk$x - x.tt)^2)
  SST = sum( (x - x.tt)^2 )
  SSE = SST - SStrt - SSblk
  return( SStrt / (k-1) / SSE * ( (k-1)*(b-1)))
}


# Complete Block Two Way ANOVA exact permutation F statistic --------------------------------

expand.row = function(pb, each){
  m = NULL
  for ( j in 1:nrow(pb)){
    m = rbind(m,matrix(rep(pb[j,],each), nrow = each, byrow = T) )
  }
  return(m)
}
repeat.mat = function(pb, copy){
  m = pb
  if (copy > 1){
    for (i in 2:copy){
      m = rbind(m,pb)
    }}
  return(m)
}

F.test2.exact = function(x,trt, blk, data = NULL){
  # exact permutations (no need to present this in class)
  # works only for THREE blocks (b = 3)
  pars <- as.list(match.call()[-1])
  if(!is.null(data)){
    x = eval(pars$data)[,as.character(pars$x)]
    treatment = eval(pars$data)[,as.character(pars$trt)]
    block = eval(pars$data)[,as.character(pars$blk)]
  }
  
  darr = data.frame(x, treatment, block) %>% arrange(block)
  x = darr$x
  treatment = darr$treatment
  block = darr$block
  trt.id = unique(treatment)
  k = length(trt.id)
  pb = permutations(k,k,trt.id)  
  allperms = cbind(
    expand.row(pb, each = (factorial(k))^2) ,
    repeat.mat( expand.row(pb, each = factorial(k)), factorial(k)),
    repeat.mat( pb, (factorial(k))^2))
  Fperms = apply(allperms, 1, function(perm.trt){
    F.stat.CBD(x,perm.trt,block)
  })
  pval = mean(Fperms >=Fobs)
  nrow(allperms)
  list(Fobs = F.stat.CBD(x, treatment, block),
       Fperms = Fperms,
       nperms = nrow(allperms),
       pvalue = pval)
}


# Complete Block Two Way ANOVA Friedman statistic --------------------------------

FM.stat.CBD = function(x, trt, blk, data = NULL){
  pars <- as.list(match.call()[-1])
  if(!is.null(data)){
    x = eval(pars$data)[,as.character(pars$x)]
    trt = eval(pars$data)[,as.character(pars$trt)]
    blk = eval(pars$data)[,as.character(pars$blk)]
  }
  x = rankinblock(x,blk)
  x.tt = mean(x)
  x.trt = aggregate(x = x, by = list(trt), FUN = mean)
  var.blk = aggregate(x = x, by = list(blk), FUN = var)
  n = length(x)
  k = nrow(x.trt)
  b = nrow(var.blk)
  
  SStrt = b * sum(   (x.trt$x - x.tt)^2)
  avg.var = mean(var.blk$x)
  return( SStrt / avg.var )
}

# Complete Block Two Way ANOVA Friedman exact test --------------------------------

Friedman.test.exact = function(x,trt, blk, data = NULL){
  # exact permutations (no need to present this in class)
  # works only for THREE blocks (b = 3)
  pars <- as.list(match.call()[-1])
  if(!is.null(data)){
    x = eval(pars$data)[,as.character(pars$x)]
    treatment = eval(pars$data)[,as.character(pars$trt)]
    block = eval(pars$data)[,as.character(pars$blk)]
  }
  
  darr = data.frame(x, treatment, block) %>% arrange(block)
  x = darr$x
  treatment = darr$treatment
  block = darr$block
  trt.id = unique(treatment)
  k = length(trt.id)
  pb = permutations(k,k,trt.id)  
  allperms = cbind(
    expand.row(pb, each = (factorial(k))^2) ,
    repeat.mat( expand.row(pb, each = factorial(k)), factorial(k)),
    repeat.mat( pb, (factorial(k))^2))
  Fperms = apply(allperms, 1, function(perm.trt){
    FM.stat.CBD(x,perm.trt,block)
  })
  pval = mean(Fperms >=Fobs)
  nrow(allperms)
  list(Fobs = FM.stat.CBD(x, treatment, block),
       Fperms = Fperms,
       nperms = nrow(allperms),
       pvalue = pval)
}

# random permutation within block -----------------------------------------
randperm.inblock = function(x, block){
  block.id = unique(block)
  b = length(block.id)
  xperm = x
  for (j in 1:b){
    xperm[block == block.id[j]] = sample( xperm[block == block.id[j]], replace = FALSE)
  }
  return(xperm)
}


# Complete Block Two Way ANOVA Friedman randperm test --------------------------------

Friedman.test.randperm = function(x,treatment, block, M = 1999, data = NULL){
  # exact permutations (no need to present this in class) 
  pars <- as.list(match.call()[-1])
  if(!is.null(data)){
    x = eval(pars$data)[,as.character(pars$x)]
    treatment = eval(pars$data)[,as.character(pars$treatment)]
    block = eval(pars$data)[,as.character(pars$block)]
  }
  Fobs = FM.stat.CBD(x, treatment, block)
  FMperm = replicate(M, {
                  xperm = randperm.inblock(x,block)
                  FM.stat.CBD(xperm, treatment, block)})
  pval = mean(c(FMobs,FMperm) >= FMobs)
  list(FMobs = Fobs,
       FMperms = FMperm,
       nperms = M,
       pvalue = pval)
} 


# Page L statistic --------------------------------------------------------

Page.L.stat = function(x,treatment.order,block,data = NULL){
  # treatment.order must consist of 1:k 
  pars <- as.list(match.call()[-1])
  if(!is.null(data)){
    x = eval(pars$data)[,as.character(pars$x)]
    treatment.order = eval(pars$data)[,as.character(pars$treatment.order)]
    block = eval(pars$data)[,as.character(pars$block)]
  }
  x = rankinblock(x,block)
  x.trt = aggregate(x = x, by = list(treatment.order), FUN = sum)
  return( sum( x.trt[,1] * x.trt[,2]) ) 
}


# Page trend rand-perm test -----------------------------------------------


Page.trend.randperm = function(x,treatment.order, block, M = 1999, data = NULL){
  # treatment.order must consist of 1:k 
  pars <- as.list(match.call()[-1])
  if(!is.null(data)){
    x = eval(pars$data)[,as.character(pars$x)]
    treatment.order = eval(pars$data)[,as.character(pars$treatment.order)]
    block = eval(pars$data)[,as.character(pars$block)]
  } 
  Lobs = Page.L.stat(x, treatment, block)
  Lperm = replicate(M, {
    xperm = randperm.inblock(x,block)
    Page.L.stat(xperm, treatment, block)})
  pval = mean(c(Lobs,Lperm) >= Lobs)
  list(Lobs = Lobs,
       Lperms = Lperm,
       nperms = M,
       pvalue = pval)
} 



# MackSkillings test statistic ------------------------------------------------------
SM.stat = function(x,treatment,block,data = NULL){
  pars <- as.list(match.call()[-1])
  if(!is.null(data)){
    x = eval(pars$data)[,as.character(pars$x)]
    treatment = eval(pars$data)[,as.character(pars$treatment)]
    block = eval(pars$data)[,as.character(pars$block)]
  }
  
  x = rankinblock(x,block)
  x.tt = mean(x)
  x.trt = aggregate(x = x, by = list(treatment), FUN = mean)
  var.blk = aggregate(x = x, by = list(block), FUN = var)
  n = length(x)
  k = nrow(x.trt)
  b = nrow(var.blk)
  c = n / (k*b)
  SStrt = b * sum(   (x.trt$x - x.tt)^2)
  var.blk
  avg.var = mean(var.blk$x)
  SStrt / avg.var
  return( SStrt / avg.var )
}

# MackSkillings test (random permutation) ------------------------------------------------------
SM.perm.test = function(x,treatment,block,M = 1999, data = NULL){
  pars <- as.list(match.call()[-1])
  if(!is.null(data)){
    x = eval(pars$data)[,as.character(pars$x)]
    treatment = eval(pars$data)[,as.character(pars$treatment)]
    block = eval(pars$data)[,as.character(pars$block)]
  }
  
  SMobs = SM.stat(x,treatment,block)
  SMperm = replicate(M, {xperm = randperm.inblock(x,block)
    SM.stat(xperm,treatment,block)})
  pval = mean(c(SMobs,SMperm) >= SMobs)
  list(SMobs = SMobs,
       SMperms = SMperm,
       nperms = M,
       pvalue = pval)
}

 
pooled.se = function(x,y){
  m = length(x)
  n = length(y)
  sqrt( sum( c( (x - mean(x))^2, (y - mean(y))^2)) / (m+n-2) * (1/m + 1/n) ) 
}




# Functions used in smoothing ---------------------------------------------

# Gaussian Kernel evaluation with bandwidth h
Kh = function(x,h){ exp( - x^2 /h^2)/ sqrt(2*pi)/h}

# Nadaraya-Watson (local constant) estimator
# for regressing y on x
# returns mhat (estimated regression function) at xgrid, and weights. 
lcreg = function(xgrid, x, y, h){
  kx = sapply(x,function(xi){ Kh(xgrid-xi,h)})
  weight <- kx / rowSums(kx)
  mhat = drop(weight %*% y) # xgrid의 각 점에서 sum w * y 계산
  return(list(mhat = mhat, weight = weight))
} 


# local linear estimator
# for regressing y on x
# returns mhat (estimated regression function) at xgrid, and weights. 
llreg = function(xgrid,x,y,h){
  nx = length(xgrid)
  mhat = xgrid
  X = model.matrix(y ~ x)
  
  weight = matrix(NA,nrow = nx, ncol = length(x)) 
  # weight w_i^1(x) matrix 
  # each row corresponds to (w_1^1(x),...w_n^1(x)),
  
  for (i in 1:nx){
    xp = xgrid[i]
    kx = Kh(xp-x,h)
    W = diag(kx)   
    t = MASS::ginv(t(X) %*% W %*% X) %*% t(X) %*% W 
    # used ginv() in place of solve() to avoid numerical instability
    beta.xp = t %*% y
    mhat[i] = beta.xp[1] + beta.xp[2] * xp
    weight[i,] = c(1, xp) %*% t
  }
  # mhat == weight %*% y (sanity check)
  return(list(mhat = mhat, weight = weight))
}


