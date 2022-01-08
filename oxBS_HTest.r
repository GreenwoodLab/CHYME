LikeoxBS_null <- function (theta, betaBS, betaOxBS, signalBS, signalOxBS) 
{
  theta <- pmin(100, pmax(-100, theta))
  p <- exp(c(0, theta))
  p <- p/sum(p)
  a1 <- signalBS * (p[2])
  b1 <- signalBS * (p[1])
  a2 <- signalOxBS * (p[2])
  b2 <- signalOxBS * (p[1])
  -(dbeta(betaBS, a1, b1, log = TRUE) + dbeta(betaOxBS, a2, 
                                              b2, log = TRUE))
}

scoreOxBS_null <- function (theta, betaBS, betaOxBS, signalBS, signalOxBS) 
{
  theta <- pmin(100, pmax(-100, theta))
  p <- exp(c(0, theta))
  p <- p/sum(p)
  a1 <- signalBS * (p[2])
  b1 <- signalBS * (p[1])
  a2 <- signalOxBS * (p[2])
  b2 <- signalOxBS * (p[1])
  dp <- matrix(c(-p[1] * p[2], p[2] * (p[1]), 0, -p[1], -p[2], 0), 
               3, 2)
  
  
  ua1 <- (diffBeta1(betaBS, a1, b1) * signalBS * (dp[2, 1] + 
                                                    dp[3, 1]) + diffBeta1(betaOxBS, a2, b2) * signalOxBS * 
            (dp[2, 1]))
  ub1 <- (diffBeta2(betaBS, a1, b1) * signalBS * (dp[1, 1]) + 
            diffBeta2(betaOxBS, a2, b2) * signalOxBS * (dp[1, 1] + 
                                                          dp[3, 1]))
  ua2 <- (diffBeta1(betaBS, a1, b1) * signalBS * (dp[2, 2] + 
                                                    dp[3, 2]) + diffBeta1(betaOxBS, a2, b2) * signalOxBS * 
            (dp[2, 2]))
  ub2 <- (diffBeta2(betaBS, a1, b1) * signalBS * (dp[1, 2]) + 
            diffBeta2(betaOxBS, a2, b2) * signalOxBS * (dp[1, 2] + 
                                                          dp[3, 2]))
  #c(ua1 + ub1, ua2 + ub2)
  ua1 + ub1
  }


fitOneOxBS_null <- function (betaBS, betaOxBS, signalBS, signalOxBS, eps = 1e-05) 
{
  est5mC <- max(min(betaOxBS, 1 - eps), eps)
  estTotMeth <- max(min(betaBS, 1 - eps), eps)
  #est5hmC <- max(estTotMeth - est5mC, eps)
  #theta <- log(c(est5mC, est5hmC)) - log(1 - est5mC - est5hmC)
  theta <- log(est5mC) - log(1 - est5mC)
  
  #gradient is different
    opt <- try(optim(theta, LikeoxBS_null, method = "BFGS", gr = scoreOxBS_null, 
                   betaBS = betaBS, betaOxBS = betaOxBS, signalBS = signalBS, 
                   signalOxBS = signalOxBS))
    
  if (inherits(opt, "try-error")) {
    out <- c(1 - estTotMeth, est5mC)
  }
  else out <- exp(c(0, opt$par))
  out/sum(out)
}


hypothesis_oxBS <- function( betaBS, betaOxBS, signalBS, signalOxBS, sig.level = 0.05 ){
  

Res_OxyBS_null <- fitOneOxBS_null( betaBS, betaOxBS ,
                         signalBS , signalOxBS )


Res_OxyBS <- fitOneOxBS( betaBS , betaOxBS ,
                         signalBS , signalOxBS )


Null_ML <- LikeoxBS_null( log(Res_OxyBS_null[2]/Res_OxyBS_null[1]), 
                          betaBS , betaOxBS ,
                          signalBS , signalOxBS )
  

All_ML <- likeOxBS( log(Res_OxyBS[-1]/Res_OxyBS[1]),
                    betaBS , betaOxBS ,
                    signalBS , signalOxBS )
  
  
LL_ratio <-  2*(Null_ML-All_ML)

rej.crit <- qchisq(1-sig.level, df=1)
  
p.value <- 1-pchisq(LL_ratio, df=1)

rej.null <-  p.value < sig.level

list( LL_ratio=LL_ratio, rej.crit=rej.crit, p.value=p.value,  rej.null=rej.null )
}



  