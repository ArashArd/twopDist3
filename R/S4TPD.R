#TPDIST FAMILY modiffied 04 2020 by arash

source("D:\\TwoPDist2\\TPLapD.R")
source("D:\\TwoPDist2\\TPLogisD.R")
source("D:\\TwoPDist2\\TPNLapD.R")
#source("E:\\Paper under review\\TwopDist_Codes\\TPNLD.R")
#source("E:\\Paper under review\\TPLogistic\\Rcode\\twoplogis.R")
 
###############################  Validity   #############################

###############################  Validity   #############################
check_inputs <- function(object) {
  errors <- character()
  length_data <- length(object@data)
  if (length_data < 15) {
    msg <- paste("data should be a vector with at least 15 elements", sep = "")
    errors <- c(errors, msg)
  }

  Var<- object@var
  if (!all(Var > 0)) {
    msg <- paste("The Variance must be positive!", sep = "")
    errors <- c(errors, msg)
  }
      
  if (length(errors) == 0) TRUE else errors
}

###############################      Main Class   #############################

###############################      Main Class   #############################


twopDist = setClass("twopDist",
		slots = list(Est = "vector", logLike = "numeric",
			    var = "vector", data = "numeric", Fitted = "numeric", KST = "vector"),
				validity = check_inputs)

twoplaplace = setClass("twoplaplace",   contains ="twopDist")

twopnormlap = setClass("twopnormlap",   contains ="twopDist")
twoplogistic = setClass("twoplogistic",   contains ="twopDist")
		#slots = c(Est = "vector", logLike = "numeric",
		#	    var = "vector", data = "numeric"),
			#	validity = check_inputs)

#new("twoplaplace", data = 1:10)

#new("twoplaplace", var = c(1, 2, -3))
###############################      Show method   #############################

###############################      Show method   #############################

 setMethod("show",
 signature = "twopDist",
 definition = function(object) {
 cat('The MLEs of the' , "\n", sep = "")
 cat("location, scale and skewpar of ", class(object), " distribution are:", "\n", sep = "")
 #cat(c("location", "scale", "skewpar"), '\n')
 cat(object@Est, "\n", sep = " ")
cat("The Log-likelihood is : ", object@logLike, "\n", sep = "")
 
cat("The approximate variance of the ", class(object), " distribution", "\n", sep = "")
cat(object@var, "\n", sep = " ")

}
)




########################## for the plot #######################
#################### SET GENERIC   #########################
setGeneric(name = "dplot"
                       , def=function(x, colhist = "gray90", ...)
                       { standardGeneric("dplot") }
                       )
##########################  plot method   #####################
#
setMethod("dplot", #f =
 signature = c(x = "twopDist"),
definition = function (x, colhist = "gray90", ...){
  dat     = x@data
  den.dat = density(dat)
  xx      = den.dat$x
  yy      = x@Fitted
  max.y   = c(max(yy, den.dat$y))
  My      = max(max.y[!(is.na(max.y))])

  par(mar = c(4.2, 4, 2, 1))

  hist(dat, freq = FALSE, col = colhist , ylim = c(0, My), ... )
  lines(den.dat, lty = 4)
  lines(xx, yy, type = "l", col = 3, lty = 1, lwd = 2)

  ################   KS-TEST   ###########################
###############################################################
mykstest = x@KST


  txt <- paste("K--S test", "\n",
             "D = ", round(mykstest[1], 4), "\n",
             "P-value = ", round(mykstest[2], 3), sep="")

 text(ifelse(median(dat) > x@Est[1], quantile(dat, 0.99), quantile(dat, 0.01)),
                quantile(den.dat$y, 0.75), txt, adj = c(0,1), pos = 3)

 legend(ifelse(median(dat) > x@Est[1], "topright", "topleft"),
          c(class(x), "Non-parametric."), col = c("green", 1),
          text.col = 1, lty = c(1, 4), pch = c(-1, -1), lwd = c(2, 1),
          bg = 'gray90',   merge = TRUE)

 rug(dat) }
 )



########################## ConfINT method   #####################
#
########################## for the Conf int #####################

#################### SET GENERIC   #########################
setGeneric(name = "ConfInt",
                       def=function(object, level= 0.95)
                       {
                               standardGeneric("ConfInt")
                       }
                       )

setMethod(f = "ConfInt",
 signature = c(object = "twopDist"),
definition = function(object, level = 0.95){
dat   = object@data
est   = object@Est
sdtpl = sqrt(object@var)
tt = paste(c("LB ", "UB "), level * 100, "% conf int", sep ="")
Mat = matrix(NA, nc = 3, nr = 3)
dimnames(Mat) = list(c("location", "scale", "skewpar"), c("MLE", tt))

Low = est + sdtpl * qnorm((1 - level)/2)
Up  = est - sdtpl * qnorm((1 - level)/2)
Mat[,1] = est;  Mat[,2] = Low; Mat[,3] = Up
print(Mat)
}
)

########################## summary method   #####################
#
########################## for the Conf int #####################

setMethod(f = "summary",
 signature = c(object = "twopDist"),
definition = function(object, level = 0.95){
dat   = object@data
est   = object@Est
sdtpl = sqrt(object@var)
tt = paste(c("LB ", "UB "), level * 100, "% conf.int", sep ="")
Mat = matrix(NA, nc = 4, nr = 3)
dimnames(Mat) = list(c("location", "scale", "skewpar"), c("MLE", "S.dev", tt))

Low = est + sdtpl * qnorm((1 - level)/2)
Up  = est - sdtpl * qnorm((1 - level)/2)
Mat[,1] = est; Mat[,2] = sdtpl; Mat[,3] = Low; Mat[,4] = Up
print(Mat)
}
)

########################## Main  MLE   #####################
#
########################## for the TPL #####################
mle_twoplaplace =
function(x, k = 1000, e.rate = 1e-10,
                         el.rate = length(x), max.iter=1000)
{
  n  = length(x)
  m  = max.iter
  pp = seq(.01, .99, length=k)
  M  = matrix (0.0, nr = k, nc = 4)
  for (j in 2 : length(pp) - 1){
    p   = pp[j]
    scale = sd(x)
    mu  = mean(x)

    t0 = c(mu = mu, scale = scale, p = p)
    t1 = t0;
    i = 1; e = 1; el = 1;
    while(i <= m && (e >= e.rate || el >= el.rate))
    {
    t1["p"]   = p
    t1["mu"]  = quantile(x, t1["p"])
    t1["scale"] = scale.tpl(x, t1["mu"], t1["p"])
    l0 = ell2.tpl(t0, x); l1 = ell2.tpl(t1, x)
    el = abs(l0 - l1); e = sum(abs(t1 - t0))
    t0 = t1;
    lt1 = ell2.tpl(t1, x)
    i = i + 1
    }
    M[j,] = c(t1, ell2.tpl(t1, x))
  }
  myname <- c("location", "scale", "skewpar", "-Log-liklihood")
  colnames(M) <- myname
  wh = which.min(M[-1000, 4])
  Ans = M[wh, ]
  #Ans
tt1 = Ans[3] * (1 - Ans[3])
varMLE = c(8* tt1 * Ans[2]^2, Ans[2]^2, tt1) / n
dat     = x
  den.dat = density(dat)
  xx      = den.dat$x
  yy      = dtpl(xx, Ans[1], Ans[2], Ans[3])

mykstest = (ks.test(dat, "ptpl", Ans[1], Ans[2], Ans[3] ))
mKST = c(mykstest$statistic, mykstest$p.value)

################   Final Answer  ###########################
###############################################################
Answer = new("twoplaplace", 
		Est = Ans[1:3], logLike = Ans[4], var = varMLE, data = x,
          Fitted = yy, KST = mKST)
Answer
}


########################## Main  MLE for the TPNL #####################
#
########################## for the TPNL #####################
mle_twopnormlap = function(x, e = 1e-10, max_iter = 100)
{
x=sort(x);
m=max_iter;
n=length(x);
ne=length(e);
if(ne==2){epr=e[1];elr=e[2];}
if(ne==1){epr=e;elr=e;}
minl=Inf;
for(k in 2:(n-1))
{
t0=c(x[k],Scale(x[k],k/n,x),k/n);
th=alg(x,t0,epr,elr,m);
l1=th[4];
if(l1<=minl)
{
theta=th;
minl=l1;
}
}

Ans = theta
varMLE = diag(solve(FI(Ans[2], Ans[3]))) / n

  dat     = x
  den.dat = density(dat)
  xx      = den.dat$x
  yy      = dtwopnormlap(xx, Ans[1], Ans[2], Ans[3])
mykstest = (ks.test(dat, "ptwopnormlap", Ans[1], Ans[2], Ans[3] ))
mKST = c(mykstest$statistic, mykstest$p.value)

################   Final Answer  ###########################
###############################################################
Answer = new("twopnormlap",
		Est = Ans[1:3], logLike = Ans[4], var = varMLE, data = x,
          Fitted = yy, KST = mKST)
Answer

}

########################## Main  MLE for the TPLogis #####################
#
##########################    for the TPLogis       #####################
mle_twoplogistic = function(x, maxit = 1e3, intialVal = NULL){
if(!is.null(intialVal)){skewpar = intialVal[3]; scale = intialVal[2]}
if(!is.null(intialVal)){if (any(skewpar <= 0 |
          skewpar >= 1 |
          scale   <= 0 ,
           na.rm = TRUE))
  stop("Intial valus are out of bound")
  }

n = length(x)
Theta = matrix(NA, ncol = 3, nrow = maxit)
tempsk = mean(x < mean(x))
tempIV = c(mean(x), sd(x), tempsk)
if(is.null(intialVal)){ Theta[1, ] = tempIV} else{
                        Theta[1, ] = intialVal}

for (i in 1:(maxit-1)){
	temp1 = Utplog(x, Theta[i,1], Theta[i,2],  Theta[i,3])#
	temp2 = IEIMtplog(Theta[i,2], Theta[i,3])
	Theta[i+1,] = Theta[i, ] + (temp2 %*%   as.matrix(temp1)/n^2)

}


Ans = Theta[maxit,]
logLike = sum(dtwoplogis(x, Ans[1], Ans[2], Ans[3], log.arg = TRUE))
varMLE = diag(IEIMtplog(Ans[2], Ans[3])) / n

  dat     = x
  den.dat = density(dat)
  xx      = den.dat$x
  yy      = dtwoplogis(xx, Ans[1], Ans[2], Ans[3])
mykstest  = (ks.test(dat, "ptwoplogis", Ans[1], Ans[2], Ans[3] ))
mKST      = c(mykstest$statistic, mykstest$p.value)

################   Final Answer  ###########################
###############################################################
Answer = new("twoplogistic",
		Est = Ans[1:3], logLike = logLike, var = varMLE, data = x,
          Fitted = yy, KST = mKST)
Answer
}

