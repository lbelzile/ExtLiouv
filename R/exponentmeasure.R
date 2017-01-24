# 
# #Exponent measure
# V.dir <- function(u, alpha, rho, B=1000){
#   sum(sapply(1:length(u),
#              function(ind){
#                temp <- exp((lgamma(alpha[-ind]+rho)-lgamma(alpha[-ind])+log(u[-ind])-
#                               (lgamma(alpha[ind]+rho)-lgamma(alpha[ind])+log(u[ind])))/rho)
#                mean(replicate(n=1000, mean(sapply(rgamma(B, shape=alpha[ind]+rho, 1) ,
#                 function(y){exp(sum(pgamma(y*temp, shape=alpha[ind], rate = 1,
#                                            lower.tail=TRUE, log.p=TRUE)))}))))/u[ind]
#                 #exp(-sum(lgamma(alpha[-ind])))
#              }))
# 
# }
# 
# 
# 
# V.negdir <- function(u, alpha, rho, B=1000){
#   sum(sapply(1:length(u),
#              function(ind){
#                temp <- exp((lgamma(alpha[ind]-rho)-lgamma(alpha[ind])+log(u[ind])-
#                               (lgamma(alpha[-ind]-rho)-lgamma(alpha[-ind])+log(u[-ind])))/rho)
#                mean(replicate(n=1000, mean(sapply(rgamma(B,shape=alpha[ind]-rho,1) ,
#                 function(y){exp(sum(pgamma(y*temp, shape=alpha[ind], rate = 1, lower.tail=FALSE, log.p=TRUE)))}))))/u[ind]
#                  #exp(-sum(lgamma(alpha[-ind])))
#              }))
# }
# 

V.dir2 <- function(x, alpha, rho, B=1000){
  k <- exp(lgamma(alpha+rho)-lgamma(alpha))
  gamma(sum(alpha)+rho)/gamma(sum(alpha))*mean(replicate(
    mean(apply(t(exp(rho*log(mev::rdir(n=B,alpha=alpha,normalize=TRUE))))/(x*k),2,max)),n=100))
}


V.dir_is <- function(x, alpha, rho, B=1000){
  k <- exp(lgamma(alpha+rho)-lgamma(alpha))
  d <- length(alpha)
  samp_num <- table(sample.int(d, size=B*100, TRUE))
  sum(sapply(1:d,  function(i){
    mixt <- rep(0, d); mixt[i] <- rho
    samp <- exp(rho*log(mev::rdir(n=samp_num[i],alpha=alpha+mixt,normalize=TRUE)))
  sum(d*apply(t(samp)/(x*k),2,max)/colSums(t(samp)/k))
  }))/(100*B)
}

V.negdir2 <- function(x, alpha, rho, B=1000){
  k <- exp(lgamma(alpha-rho)-lgamma(alpha))
  gamma(sum(alpha)-rho)/gamma(sum(alpha))*mean(replicate(
    mean(apply(t(exp(-rho*log(mev::rdir(n=B,alpha=alpha,normalize=TRUE))))/(x*k),2,max)),n=100))
}
