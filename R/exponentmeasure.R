
#Exponent measure
V.dir <- function(u, alpha, rho){
  sum(sapply(1:length(u),
             function(ind){
               temp <- exp((lgamma(alpha[-ind]+rho)-lgamma(alpha[-ind])+u[-ind]-
                              (lgamma(alpha[ind]+rho)-lgamma(alpha[ind])+u[ind]))/rho)
               mean(replicate(n=1000, mean(sapply(rgamma(10000,shape=alpha[ind]+rho,1) , 
                                                  function(y){exp(sum(pgamma(y*temp, shape=alpha[ind], rate = 1, lower.tail=TRUE, log.p=TRUE)))}))))/u[ind]*
                 exp(-sum(lgamma(alpha[-ind])))
             }))
  
}

V.negdir <- function(u, alpha, rho){
  sum(sapply(1:length(u),
             function(ind){
               temp <- exp((lgamma(alpha[ind]-rho)-lgamma(alpha[ind])+u[ind]-
                              (lgamma(alpha[-ind]-rho)-lgamma(alpha[-ind])+u[-ind]))/rho)
               mean(replicate(n=1000, mean(sapply(rgamma(10000,shape=alpha[ind]-rho,1) , 
                                                  function(y){exp(sum(pgamma(y*temp, shape=alpha[ind], rate = 1, lower.tail=FALSE, log.p=TRUE)))}))))/u[ind]*
                 exp(-sum(lgamma(alpha[-ind])))
             }))
}
