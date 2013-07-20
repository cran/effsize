
## Data from:
##     M.Hess, J.Kromrey. 
##     Robust Confidence Intervals for Effect Sizes: 
##     A Comparative Study of Cohen’s d and Cliff’s Delta Under Non-normality and Heterogeneous Variances
##     American Educational Research Association, San Diego, April 12 – 16, 2004
##
## Available at: http://www.coedu.usf.edu/main/departments/me/documents/cohen.pdf
##
## Kristine Y. Hogarty and Jeffrey D. Kromrey
## Using SAS to Calculate Tests of Cliff’s Delta
## Proceedings of the Twenty-Fourth Annual SAS® Users Group ￼International Conferenc, Miami Beach, Florida, 1999
## Availabel at: http://www2.sas.com/proceedings/sugi24/Posters/p238-24.pdf

cliff.delta <- function(treatment, ... ) UseMethod("cliff.delta")

cliff.delta.default <- function( treatment, control, conf.level=.95, 
                         use.unbiased=TRUE, use.normal=FALSE, return.dm=FALSE, ...){

  if(conf.level>.9999 | conf.level<.5) stop("conf.level must be within 50% and 99.99%")
  if(class(control)=="factor" | class(control)=="character"){
     f = control
     if(length(f)!=length(treatment)){
         stop("Factor must have the same length as the data")
     }
     if(length(unique(f))!=2){
        stop("The factor must have exactly two levels")
     }
     vals = split(treatment,f)
     treatment = vals[[1]]
     control = vals[[2]]
  }
  
  treatment = sort(treatment)
  control = sort(control)
  dominance = sign(outer(treatment, control, FUN="-"))
  row.names(dominance) = treatment
  colnames(dominance) = control

d = mean(dominance)

n1 = length(treatment)
n2 = length(control)


d_i. = apply(dominance,1,mean)
d_.j = apply(dominance,2,mean)

if(use.unbiased){
# method 1: unbiased estimate:
S_d = ( n2^2 * sum( (d_i. - d)^2) + n1^2*sum( (d_.j-d)^2) - sum( (dominance-d)^2 ) ) / (n1*n2*(n1-1)*(n2-1))
}else{
# method 2: consistent estimate
S_i. = sum( (apply(dominance,1,mean)-d)^2 ) / (n1-1);
S_.j = sum( (apply(dominance,2,mean)-d)^2 ) / (n2-1);
S_ij = sum( (dominance-d)^2 ) / ( ( n1-1)*(n2-1 ) ) ### Cliff 1996
#S_ij = sum( (dominance-delta)^2 ) / ( n1*n2-1 )  ## Long et al 2003
S_d = ( (n2-1)*S_i. + (n1-1)*S_.j + S_ij) / ( n1 * n2)
}

if(use.normal){
  ## assume a normal distribution
Z = -qnorm((1-conf.level)/2)
}else{
  ## assume a Student t distribution See (Feng & Cliff, 2004) 
Z = -qt((1-conf.level)/2,n1+n2-2)
}
conf.int = c(
  ( d - d^3 - Z * sqrt(S_d) * sqrt((1-d^2)^2+Z^2*S_d )) / ( 1 - d^2+Z^2*S_d),
  ( d - d^3 + Z * sqrt(S_d) * sqrt((1-d^2)^2+Z^2*S_d )) / ( 1 - d^2+Z^2*S_d)
)
names(conf.int) = c("inf","sup")
  
  levels = c(0.147,0.33,0.474)
  magnitude = c("negligible","small","medium","large")
res= 
 list(
  estimate = d,
  conf.int = conf.int,
  var = S_d,
  conf.level = conf.level,
  magnitude = magnitude[findInterval(abs(d),levels)+1],
  method = "Cliff's Delta",
  variance.estimation = if(use.unbiased){ "Unbiased"}else{"Consistent"},
  CI.distribution = if(use.normal){ "Normal"}else{"Student-t"}
  )
  if(return.dm){
    res$dm = dominance;
  }
  res$name = "delta"
  
  class(res) <- "effsize"
  return(res)
}


cliff.delta.formula <-function(formula, data=list(),conf.level=.95, 
                                use.unbiased=TRUE, use.normal=FALSE, 
                               return.dm=FALSE, ...){
#cliff.delta.formula <-function(f, data=list(),...){
  mf <- model.frame(formula=formula, data=data)
  if(dim(mf)[2]!=2){
     stop("Fomula must be a variable vs a factor")
  }
  x <- mf[[2]]
  y <- mf[[1]]
  res = cliff.delta.default(y,x,conf.level, use.unbiased, use.normal, return.dm)
  return(res)
}

print.effsize <- function(x, ...){
  cat("\n")
  cat(x$method)
  cat("\n\n")
  cat(x$name)
  cat(" estimate: ")
  cat(x$estimate)
  cat(" (")
  cat(x$magnitude)
  cat(")\n")
  if("conf.level" %in% names(x)){
    conf = x$conf.level*100
    cat(conf)
    cat(" percent confidence interval:\n")
    print(x$conf.int)
  }
}

# test_cliff.delta <- function(){
#   treatment = c(10,10,20,20,20,30,30,30,40,50)
#   control = c(10,20,30,40,40,50)
#   
#   res = cliff.delta(treatment,control,use.unbiased=F)
#   
#   print(res)
#   
# }


