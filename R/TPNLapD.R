library(stats4)


# NORMAL LAPLACE FAMILY
# Reference: Ardalan et. al, The Two-Piece Normal-Laplace Distribution paper
# (Comm. in Statist.: Theory & Methods. preprint)
#=============================================================
#           dtwopnormlap
#=============================================================


dtwopnormlap <- function(x, location = 0, scale = 1, skewpar = 0.5,
                 log.arg = FALSE) {
  if (any(skewpar <= 0 |
          skewpar >= 1 |
          scale   <= 0 ,
           na.rm = TRUE))
  stop("some parameters out of bound")

# Recycle the vectors to equal lengths
  LLL = max(length(x), length(location),
            length(scale), length(skewpar))
  if (length(x)        != LLL)        x = rep(x,        length.out = LLL)
  if (length(location) != LLL) location = rep(location, length.out = LLL)
  if (length(scale)    != LLL)    scale = rep(scale,    length.out = LLL)
  if (length(skewpar)  != LLL)  skewpar = rep(skewpar,  length.out = LLL)

  zedd <- (x - location) / scale
  y = -zedd^2 / (8 * skewpar^2)
  y[zedd >  0] = (-zedd / (sqrt(2 * pi) * (1 - skewpar)))[zedd > 0]
  if(log.arg) {
    y - 0.5 * log(2 * pi) - log(scale)
  } else {
    1 / (sqrt(2 * pi) * scale) * exp(y)
  }
}

#=============================================================
#         ptwopnormlap
#=============================================================


ptwopnormlap <- function(q, location = 0, scale = 1, skewpar = 0.5) {
  if (any(skewpar <= 0 |
          skewpar >= 1 |
          scale   <= 0 ,
          na.rm = TRUE))
    stop("some parameters out of bound")

# Recycle the vectors to equal lengths
  LLL = max(length(q), length(location),
            length(scale), length(skewpar))
  if (length(q)        != LLL)        q = rep(q,        length.out = LLL)
  if (length(location) != LLL) location = rep(location, length.out = LLL)
  if (length(scale)    != LLL)    scale = rep(scale,    length.out = LLL)
  if (length(skewpar)  != LLL)  skewpar = rep(skewpar,  length.out = LLL)

  zedd <- (q - location) / scale
  pR <- (1 - skewpar)

  s1 <- 2 * skewpar * pnorm(zedd, sd = 2 * skewpar) #/ scale
  s2 <- skewpar +  pR *(1 - exp(-zedd / (sqrt(2 * pi) * pR)))

  ans <- rep(0.0, length(zedd))
  ans[zedd <= 0] <- s1[zedd <= 0]
  ans[zedd >  0] <- s2[zedd >  0]
  ans
}



#============================================================
#         qtwopnormlap
#============================================================

# Pos.  <- function(x) ifelse(x >  0, x, 0.0)
# Lone. <- function(x) ifelse(x <= 1, x, 0.0)

qtwopnormlap <- function(p, location = 0, scale = 1, skewpar = 0.5) {

# posfun <- function(x) pmax(x, 0)
# Lone.  <- function(x) ifelse(x <= 1, x, 0.0)

  Pos.  <- function(x) ifelse(x >  0, x, 0.0)

  pp = p
  if (any(pp      <= 0 |
          pp      >= 1 |
          skewpar <= 0 |
          skewpar >= 1 |
          scale   <= 0 ,
             na.rm = TRUE))
  stop("some parameters out of bound")

# Recycle the vectors to equal lengths
  LLL = max(length(pp), length(location),
            length(scale), length(skewpar))
  if (length(pp)       != LLL)       pp = rep(pp,       length.out = LLL)
  if (length(location) != LLL) location = rep(location, length.out = LLL)
  if (length(scale)    != LLL)    scale = rep(scale,    length.out = LLL)
  if (length(skewpar)  != LLL)  skewpar = rep(skewpar,  length.out = LLL)

  pL <- skewpar
  pR <- (1 - skewpar)
  QL = pp / (2 * pL);
  QL = ifelse(QL >= 1 | QL <= 0, 0.1, QL)

  qtpL <- qnorm(QL, mean = 0, sd = 2 * skewpar) * (pp <= pL)
# qtpR <- (-sqrt(2*pi) * pR * log(1 - (Pos.(pp) - pL) / pR)) * (pp > pL)
  qtpR <- (-sqrt(2*pi) * pR * log1p(- (Pos.(pp) - pL) / pR)) * (pp > pL)

  qtp <- qtpL + qtpR
  qtp * scale + location
}


#=============================================================
#         rtwopnormlap
#=============================================================

rtwopnormlap <- function(n, location = 0, scale = 1, skewpar = 0.5) {
  qtwopnormlap(p = runif(n), location = location, scale = scale,
               skewpar = skewpar)
}





#############################################
dnl = function(x = dat, location = 0, scale = 1, skewpar = 0.5, log = FALSE)
{
n = length(x)
if(scale <= 0||skewpar <= 0||skewpar >= 1)
if(log == TRUE) return(rep(-10000,n))
else return(rep(0,n))
else
{
z = (x-location)/scale
y = z
z1 = z[z <= 0]
z2 = z[z > 0]
y[z<= 0]=-z1^2/(8*skewpar^2)
y[z>0]=-z2/(sqrt(2*pi)*(1-skewpar))
if(log==TRUE) y=y-1/2*log(2*pi)-log(scale)
else y=1/(sqrt(2*pi)*scale)*exp(y)
}
y

}
################################################
qnl<-function(skewpar, location = 0,scale=1,skewpara=.5){
gskewpar<-c(1/(8*skewpar^2),1/(sqrt(2*pi)*(1-skewpar)))
skewpar1=skewpar[skewpar<=skewpara]
skewpar2=skewpar[skewpar>skewpara]
q=skewpar
q[skewpar<=skewpara]=2*skewpara*qnorm(skewpar1/(2*skewpara))
q[skewpar>skewpara]=-sqrt(2*pi)*(1-skewpara)*log(1-(skewpar2-skewpara)/(1-skewpara))
q*scale+location
}
###################################

rnl<-function(n,location=0,scale=1,skewpar=.5){
u<-runif(n)
qnl(u,location,scale,skewpar)
}
###################################################
L = function(location, scale, skewpar, x = dat) -sum(dnl(x, location, scale, skewpar, TRUE))

L2 = function(theta, x = dat)
L(theta[1], theta[2], theta[3], x)
###################################################
I.location=function(location,x=dat)
c(sum((x-location)^2*(x<=location)),sum((x-location)*(x>location)))
###############################################
Scale = function(location,skewpar,x=dat)
{
n=length(x)
i=I.location(location,x)
i1=i[1];i2=i[2]
a=(4*n*skewpar^2)*(n*(1-skewpar)*sqrt(2*pi))
b=-i2*(4*n*skewpar^2)
c=-i1*(n*(1-skewpar)*sqrt(2*pi))
D=b^2-4*a*c
if(D<0)
{cat('Warning: Delta for scale is negative\n')
return(.0001)
}
if(a==0) return(0.0001)
return((-b + sqrt(D))/(2 * a))
}
###############################################
SKEWPAR = function(location, scale, x=dat)
{
n=length(x)
i.location=I.location(location,x)
i1=i.location[1];i2=i.location[2]
a=4*scale^2*i2
b=-scale*sqrt(2*pi)*i1
skewpar=root.skewpar(function(x) a*x^3+b*(1-x)^2)
skewpar
}

##############################  SKEWPARu  #################

SKEWPARu = function(location, scale, x = dat)
{
n = length(x)
i.location = I.location(location,x)
i1 = i.location[1]; i2 = i.location[2]
a = 4 * scale^2 * i2
b = -scale * sqrt(2 * pi) * i1

skewpar=uniroot(function(x) a*x^3+b*(1-x)^2, c(0,1), tol = 1e-10 )[[1]]
skewpar
}


###############################################
Location=function(scale,skewpar,x=dat)
{
if(skewpar <= 0) return(min(x))
if(skewpar >= 1) return(max(x))
n=length(x);lx=rep(0,0);x=sort(x)
a=8*skewpar^2*scale^2
b=1/(sqrt(2*pi)*(1-skewpar)*scale)
for(i in 1:n)lx[i]=L(x[i],scale,skewpar,x)
llocation=min(lx)
location.hat=x[which.min(lx)]
for(i in 1:(n-1))
{
m.hat=(n-i)*b/(2*i)*a+mean(x[1:i])
if(m.hat>x[i]&m.hat<x[i+1])
{
l=L(m.hat,scale,skewpar,x)
#if(llocation>l) {
#cat(m.hat,"location.hat is not a point","\n")
#llocation=l;location.hat=m.hat;j=i}
lm=l
}
}
location.hat
}
################################################
mle.alg=function(x=dat,e.rate=1e-10,el.rate=length(x),max.iter=1000)
{
n=length(x)
m=max.iter
skewpar=mean(x<=mean(x))
c1=(4*skewpar^2/sqrt(2*pi)-(1-skewpar)^2*sqrt(2*pi))
c2=1/((4*skewpar^3+4*(1-skewpar)^3)*pi-((1-skewpar)^2*sqrt(2*pi)-(4*skewpar^2)/(sqrt(2*pi)))^2)
scale=sqrt(c2)*sd(x)

location=mean(x)+c1*sd(x)
#t0=c(location=mean(x),scale=sd(x),skewpar=sum(x<=mean(x))/n)
t0=c(location=location,scale=scale,skewpar=skewpar)
t1=t0;
i=1;e=1;el=1;
while(i<=m&&(e>=e.rate||el>=el.rate))
#while(i!=100)
{
t1['skewpar']=SKEWPAR(t1['location'],t0['scale'],x)
t1['scale']=Scale(t1['location'],t1['skewpar'],x)
t1['location']=Location(t1['scale'],t1['skewpar'],x)
l0=L2(t0);l1=L2(t1)
cat(i,t1,l1,"\n")
el=abs(l0-l1);e=sum(abs(t1-t0))
t0=t1;
lt1=L2(t1)
i=i+1
}
if(i>=m) cat('Warning: Max iteration reached\n')
t1
}
###############################################
diff.L=function(location,scale,skewpar,x=dat)
{
y=rep(0,0)
i=I.location(location,x);i1=i[1];i2=i[2]
y11=1/(4*skewpar^2*scale^2)*sum((x-location)*(x<=location))
y12=sum(x>location)/(scale*(1-skewpar)*sqrt(2*pi))
y[1]=y11+y12
y[2]=scale^2-scale*i1/(4*n*skewpar^2)-i2/(n*(1-skewpar)*sqrt(2*pi))
y[3]=4*scale^2*skewpar^3*i2-scale*(1-skewpar)^2*sqrt(2*pi)*i1
y
}
################################################
root.skewpar=function(f,bound=c(0,1),max.iter=40)
{
a=bound[1]
b=bound[2]
m=max.iter
ab=mean(c(a,b))
if(f(a)==0) return(a)
if(f(b)==0) return(b)
if(f(a)*f(b)>0) stop('No roots were found')
for(i in 1:m)
{
fab=f(ab)
if(fab==0) return(ab)
fa=f(a)
if(fa*fab>0) a=ab
else b=ab
ab=mean(c(a,b))
}
ab
}
############################################## The psi function


gmsp=function(location,scale,skewpar,x=dat){
n=length(x)
Ilocation=sum(c((x-location)/(4*skewpar^2*scale^2)*(x<=location)  , 1/(scale*(1-skewpar)*sqrt(2*pi))*(x>location)))
I1=sum((x-location)^2*(x<=location))
I2=sum((x-location)*(x>location))

gscale=scale^2-(scale*I2)/(n*(1-skewpar)*sqrt(2*pi))-I1/(4*n*skewpar^2)

gskewpar=4*scale^2*skewpar^3*I2-scale*(1-skewpar)^2*sqrt(2*pi)*I1
eq=list("Ilocation"=Ilocation,"gscale"=gscale,"gskewpar"=gskewpar)
eq
}
###################################################
mle.alg.skewpars=function(x=dat,location=0,e.rate=1e-10,el.rate=length(x),max.iter=1000)
{
n=length(x)
m=max.iter
#skewpar=sum(x<=quantile(x,.65))/n
#c1=(4*skewpar^2/sqrt(2*pi)-(1-skewpar)^2*sqrt(2*pi))
#c2=1/((4*skewpar^3+4*(1-skewpar)^3)*pi-((1-skewpar)^2*sqrt(2*pi)-(4*skewpar^2)/(sqrt(2*pi)))^2)
#scale=sqrt(c2)*sd(x)

#location=mean(x)+c1*sd(x)
t0=c(location=location,scale=sd(x),skewpar=sum(x<=mean(x))/n)
#t0=c(location=location,scale=scale,skewpar=skewpar)
t1=t0;
i=1;e=1;el=1;
while(i<=m&&(e>=e.rate||el>=el.rate))
#while(i!=100)
{
t1['skewpar']=SKEWPAR(t1['location'],t0['scale'],x)
t1['scale']=Scale(t1['location'],t1['skewpar'],x)
t1['location']=location#Location(t1['scale'],t1['skewpar'],x)
l0=L2(t0);l1=L2(t1)
#cat(i,t1,l1,"\n")
el=abs(l0-l1);e=sum(abs(t1-t0))
t0=t1;
lt1=L2(t1)
i=i+1
}
if(i>=m) cat('Warning: Max iteration reached\n')
c(t1,lt1)
}
#####################################   new alg ######################
mle.algn=function(x=dat,n1=20)
{
Like=rep(0,0)
n=length(x)
x=sort(x)
skewpar=sum(x<=quantile(x,.65))/n
c1=(4*skewpar^2/sqrt(2*pi)-(1-skewpar)^2*sqrt(2*pi))
c2=1/((4*skewpar^3+4*(1-skewpar)^3)*pi-((1-skewpar)^2*sqrt(2*pi)-(4*skewpar^2)/(sqrt(2*pi)))^2)
scale=sqrt(c2)*sd(x)

location=mean(x)+c1*sd(x)
t1=matrix(ncol=3,nrow=n1)
t1[1,]=c(location,scale,2/n)
t1[1,]=c(.907251261,1.0174235380,.8782498)


for( i in 2 : n1)
{

t1[i,1]=Location(t1[i-1,2],t1[i-1,3],x)  #Location(t1[i-1,2],t1[i-1,3],x)

t1[i,3]= i/n #SKEWPAR(t1[i,1],t1[i-1,2],x)
t1[i,2]=Scale(t1[i,1],t1[i,3],x)
Like[i]=l1=L2(t1[i,])
cat(i,t1[i,],l1,"\n")
}
m=which.min(Like[-c(1,n1)])
skewpars=t1[m+1,]
return(skewpars)
}

######################    mle.alg with start value new alg  #########

mle.algnn=function(x=dat,e.rate=1e-10,el.rate=length(x),max.iter=1000)
{
n=length(x)
m=max.iter
skewpar=sum(x<=quantile(x,.65))/n
c1=(4*skewpar^2/sqrt(2*pi)-(1-skewpar)^2*sqrt(2*pi))
c2=1/((4*skewpar^3+4*(1-skewpar)^3)*pi-((1-skewpar)^2*sqrt(2*pi)-(4*skewpar^2)/(sqrt(2*pi)))^2)
scale=sqrt(c2)*sd(x)

location=mean(x)+c1*sd(x)
t0=mle.algn(x,n)
#t0=c(location=location,scale=scale,skewpar=skewpar)
t1=t0;
i=1;e=1;el=1;
#while(i<=m&&(e>=e.rate||el>=el.rate))
while(i!=100)
{
t1['skewpar']=SKEWPAR(t1['location'],t0['scale'],x)
t1['scale']=Scale(t1['location'],t1['skewpar'],x)
t1['location']=Location(t1['scale'],t1['skewpar'],x)
l0=L2(t0);l1=L2(t1)
cat(i,t1,l1,"\n")
#el=abs(l0-l1);e=sum(abs(t1-t0))
t0=t1;
#lt1=L2(t1)
i=i+1
}
if(i>=m) cat('Warning: Max iteration reached\n')
t1
}
#################################################
mle.algs=function(x=dat,scale=scale,e.rate=1e-10,el.rate=length(x),max.iter=1000)
{
n=length(x)
m=max.iter
skewpar=sum(x<=quantile(x,.65))/n
c1=(4*skewpar^2/sqrt(2*pi)-(1-skewpar)^2*sqrt(2*pi))
c2=1/((4*skewpar^3+4*(1-skewpar)^3)*pi-((1-skewpar)^2*sqrt(2*pi)-(4*skewpar^2)/(sqrt(2*pi)))^2)
#scale=sqrt(c2)*sd(x)
scale=scale
location=mean(x)+c1*sd(x)
#t0=c(location=mean(x),scale=sd(x),skewpar=sum(x<=mean(x))/n)
t0=c(location=location,scale=scale,skewpar=skewpar)
t1=t0;
i=1;e=1;el=1;
while(i<=m&&(e>=e.rate||el>=el.rate))
#while(i!=100)
{
t1['skewpar']=SKEWPAR(t1['location'],t0['scale'],x)
t1['scale']=scale
#Scale(t1['location'],t1['skewpar'],x)
t1['location']=Location(t1['scale'],t1['skewpar'],x)
l0=L2(t0);l1=L2(t1)
#cat(i,t1,l1,"\n")
el=abs(l0-l1);e=sum(abs(t1-t0))
t0=t1;
lt1=L2(t1)
i=i+1
}
if(i>=m) cat('Warning: Max iteration reached\n')
t1
}
################################## 27/12/2010     ##########################
mle.alg.skewpar=function(x=dat,k=1000,e.rate=1e-10,el.rate=length(x),max.iter=1000)
{
n = length(x)
m = max.iter
skewparskewpar = seq(0,1,length=k)
M = matrix (0.00000,nr=k, nc=4)
for (j in 2: length(skewparskewpar)-1){
skewpar=skewparskewpar[j]
c1=(4*skewpar^2/sqrt(2*pi)-(1-skewpar)^2*sqrt(2*pi))
c2=1/((4*skewpar^3+4*(1-skewpar)^3)*pi-((1-skewpar)^2*sqrt(2*pi)-(4*skewpar^2)/(sqrt(2*pi)))^2)
scale=sqrt(c2)*sd(x)
location = mean(x)+c1*sd(x)

t0=c(location=location,scale=scale,skewpar=skewpar)
t1=t0;
i=1;e=1;el=1;
while(i<=m&&(e>=e.rate||el>=el.rate))
{
t1['skewpar']= skewpar
t1['scale']=Scale(t1['location'],t1['skewpar'],x)
t1['location']=Location(t1['scale'],t1['skewpar'],x)
l0=L2(t0);l1=L2(t1)
el=abs(l0-l1);e=sum(abs(t1-t0))
t0=t1;
lt1=L2(t1)
i=i+1
}
M[j,]= c(t1,L2(t1))
}
wh = which.min(M[,4])
Ans = M[wh,]
Ans
#M
}
##############################   for SMSA alg  #######################
mle.algmskewpar=function(x=dat,location,skewpar,e.rate=1e-10,el.rate=length(x),max.iter=1000)
{
n=length(x)
m=max.iter
skewpar=skewpar
c1=(4*skewpar^2/sqrt(2*pi)-(1-skewpar)^2*sqrt(2*pi))
c2=1/((4*skewpar^3+4*(1-p)^3)*pi-((1-skewpar)^2*sqrt(2*pi)-(4*skewpar^2)/(sqrt(2*pi)))^2)
scale=sqrt(c2)*sd(x)

location=location
t0=c(location=location,scale=scale,skewpar=skewpar)
t1=t0;
i=1;e=1;el=1;
while(i<=m&&(e>=e.rate||el>=el.rate))
#while(i!=100)
{
t1['skewpar']=SKEWPAR(t1['location'],t0['scale'],x)
t1['scale']=Scale(t1['location'],t1['skewpar'],x)
t1['location']=Location(t1['scale'],t1['skewpar'],x)
l0=L2(t0);l1=L2(t1)
#cat(i,t1,l1,"\n")
el=abs(l0-l1);e=sum(abs(t1-t0))
t0=t1;
lt1=L2(t1)
i=i+1
}
if(i>=m) cat('Warning: Max iteration reached\n')
c(t1,L2(t1))
}
###################3 final for smsa alg ####################33
mle.alg.smsa = function(dat){
skewpars = matrix(ncol=4,nrow=length(dat))
dat=sort(dat)
for(i in 2:length(dat))
skewpars[i,] = mle.algmskewpar(dat,location=dat[i],skewpar=i/length(dat))

skewparsmain=skewpars[which.min(skewpars[,4]),]
skewparsmain}
################################   Algorithm     #############################

alg=function(x, t0, epr, elr, m)
{
i = 1; e = 1;el = 1;t1 = t0;
while(i<=m&&(e>=epr||el>=elr))
{
t1[1]=Location(t0[2],t0[3],x);
t1[2]=Scale(t1[1],t0[3],x);
t1[3]=SKEWPAR(t1[1],t1[2],x);
l0=L2(t0,x);
l1=L2(t1,x);
el=abs(l0-l1);
e=sum(abs(t1-t0));
t0=t1;
i=i+1;
}
if(i>=m) cat('Warninig: Reached max iteration');
c(t1,l1);
}



#############################################################################

############################### MLE ALG NEW       ##############################


mle.alg.n =function(x,e=1e-10,max_iter=100,l=0)
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
if(l==0) theta=theta[1:3];
theta
}
###############################################
################################################
root.skewpar = function(f, bound = c(0,1), e.x = 1e-12, e.f = 1e-10)
{
a=bound[1];b=bound[2]
ab=mean(c(a,b))
if(f(a)==0) return(a)
if(f(b)==0) return(b)
if(f(a)*f(b)>0) stop('No roots were found')
e.xx=1;e.ff=1
while(e.xx>=e.x||e.ff>=e.f)
{
fab=f(ab)
if(fab==0) return(ab)
fa=f(a)
if(fa*fab>0) a=ab
else b=ab
ab=mean(c(a,b))
e.xx=b-a
e.ff=abs(f(ab))
}
ab
}


################################  fisher information ####################

FI = function(scalema,skewpar){
IL = (1/(4*skewpar)+1/(((1-skewpar)*2)*pi))/scalema^2
ILS = -1/(scalema^2*sqrt(2*pi))
ILSK = (-2+skewpar)/(scalema*skewpar*(1-skewpar)*sqrt(2*pi))
IS = (skewpar+1)/scalema^2
ISSK = 1/scalema
ISK = 3/skewpar+2/(1-skewpar)
F=matrix(c(c(IL, ILS, ILSK),
c(ILS, IS, ISSK),
c(ILSK, ISSK, ISK)),ncol=3)
F
}

