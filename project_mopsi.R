library(grDevices)
library(RColorBrewer)
require(plotrix)

trajectoire <- function(x,h,Te,Ti,r){
  iter = 0;
  #print(x);
  X=x;
  Y=x;
  #print(abs(X[1]))
  
  o1=c(-0.5,0.5);
  o2=c(0.5,0.5);
  o3=c(0.5,-0.5);
  o4=c(-0.5,-0.5);
  while (abs(X[1])<=1 && abs(X[2])<=1 && sqrt(sum((X-o1)*(X-o1)))>r && sqrt(sum((X-o2)*(X-o2)))>r && sqrt(sum((X-o3)*(X-o3)))>r && sqrt(sum((X-o4)*(X-o4)))>r ){
    X=X+sqrt(h)*rnorm(2,0,1);
    Y=c(Y,X);
    iter = iter + 1
  }
  #print(iter)
  if( (abs(X[1])>1) || (abs(X[2])>1)){
    return (c(Te,Y));
  }else {
    return (c(Ti,Y));
  }
}

h=0.01;
Te=0;
Ti=0.5;
x <- c(0,0);
r=0.1;
#tracer une trajectoire
affich_traj <- function(x,h,Te,Ti,r){
  y=trajectoire(x,h,Te,Ti,r);
  print(y);
  ind_imp = seq(from = 2,to = length(y), by = 2);
  x = y[ind_imp];
  ind_p = seq(from = 3,to = length(y), by = 2);
  z = y[ind_p];
  plot(x,z,type='l',axes = FALSE,col = "RED")
  axis(1,  pos = 0, cex.axis = 0.8)
  axis(2,  pos = 0, cex.axis = 0.8, las = 2)
}
#affich_traj(x,h,Te,Ti,r)
#temperature
temperature <- function(x,h,Te,Ti,r,N){
  X <- list(x);
  #print(X)
  vec <- rep(X,N);
  l <- sapply(vec,function(a) trajectoire(a,h,Te,Ti,r)[1]);
  #print(l)
  return (mean(l));
}
temp <- temperature(x,0.01,Te,Ti,r,500);
print(temp);
#faire la grille abline pour des lignes
#comment faire une couleur avec un entier
#ggplot2 ?
seq1 <- seq(from = -1, to = 1,by=0.002);
seq2<- seq(from = -1, to = 1,by=0.002);
#grille <- data.matrix(merge(data.frame(seq1),data.frame(seq2)));
grille <- data.matrix(merge(seq(from = -1, to = 1,by=0.01),seq(from = -1, to = 1,by=0.01)));
#list <-
M <-apply(grille,1,function(a) temperature(a,0.01,Te,Ti,r,50));

#concatener le vecteur des palettes
#> pie(rep(1,16),col=heat.colors(15))
#> pie(rep(1,50),col=heat.colors(50))
#> palette(heat.colors(12))
#> palette()
#ou alors :
#display.brewer.pal(n = 11, name = 'RdYlBu')
require(akima)
data <- data.frame(x = grille[,1],y = grille[,2],distance = M)
resolution <- 0.1 # you can increase the resolution by decreasing this number (warning: the resulting dataframe size increase very quickly)
a <- interp(x=data$x, y=data$y, z=data$distance, 
            xo=seq(min(data$x),max(data$x),by=resolution), 
            yo=seq(min(data$y),max(data$y),by=resolution), duplicate="mean")
#image(a) #you can of course modify the color palette and the color categories. See ?image for more explanation
filled.contour(a, col=palette())

r=c(0.1,0.1,0.1,0.1)
ax=c(-0.5,-0.5,0.5,0.5)
by=c(-0.5,0.5,-0.5,0.5)
plot(c(-1,1), c(-1,1), type= "n", xlab = "", ylab = "")
#draw.circle(0,0,0.02,col='blue',border='blue')
rect(-1,-1,1,1, border="blue")
draw.circle(ax,by,r,border='red')
y <- trajectoire(x,0.001,Te,Ti,r);
ind_imp = seq(from = 2,to = length(y), by = 2);
z1 <- y[ind_imp];
ind_p = seq(from = 3,to = length(y), by = 2);
z2 <- y[ind_p];
lines(z1,z2);

temperature2 <- function(x,h,Te,Ti,r,N){
  X <- list(x);
  #print(X)
  vec <- rep(X,N);
  l <- sapply(vec,function(a) trajectoire(a,h,Te,Ti,r)[1]);
  #print(mean(l))
  #print(sd(l))
  return (c(mean(l),sd(l)));
}
#intervalle de confiance pour en fonction de Ti
n = 500
h=0.001
u_0 <- temperature2(x,h,Te,Ti,r,n)
print(u_0[2]/sqrt(n)* qnorm(0.995)+ u_0[1])
print(-u_0[2]/sqrt(n)* qnorm(0.995)+ u_0[1])
print(u_0[1])
Temp = seq(from = 0,to = 1, by = 0.01);
l <- sapply(Temp,function(a) temperature2(x,h,Te,a,r,n));
plot(Temp,l[1,],type='l')
lines(Temp,l[2,]/sqrt(n)* qnorm(0.995)+ l[1,])
lines(Temp,-l[2,]/sqrt(n)* qnorm(0.995)+ l[1,])

#evolution de u(x) selon le rayon des disques
ray = seq(from = 0.01,to = 0.49,by=0.01)
l <- sapply(ray,function(r) temperature2(x,h,Te,Ti,r,n));
plot(ray,l[1,],type='l')

trajectoire_2 <- function(x,h,Ti,alpha){
  X=x;
  Y=x;
  while (sum(X*X)<1 ){
    X=X+sqrt(h)*rnorm(2,0,1);
    Y=c(Y,X);
  }
    return (c(X[1]+h*alpha*sum(Y*Y),X[1]-Ti));
}
alpha = -1
temperature_2 <- function(x,h,Ti,N,alpha){
#donne u(x) par trois estimateurs différents et les variances associées
    X <- list(x);
    vec <- rep(X,N);
    l <- sapply(vec,function(a) trajectoire_2(a,h,Ti,alpha));
    lambda = cov(l[1,],l[1,]-l[2,])/var(l[2,]); 
    return (c(mean(l[1,]),mean(l[2,])-alpha*sum(x*x)*sum(x*x)/8+alpha/8+Ti,mean((1-lambda)*l[1,]+lambda*l[2,])+lambda*(-alpha*sum(x*x)*sum(x*x)/8+alpha/8+Ti),var(l[1,]),var(l[2,]),var((1-lambda)*l[1,]+lambda*l[2,])));
  
    }

trajectoire_3 <- function(x,h,Ti,r,alpha){
  X=x;
  Y=x;
  while (sum(X*X)<1 ){
    X=X+sqrt(h)*rnorm(2,0,1);
    Y=c(Y,X);
  }
  return (X[1]-Ti);
}
# temperature_3 <- function(x,h,Te,Ti,r,N,alpha){
#   X <- list(x);
#   vec <- rep(X,N);
#   l <- sapply(vec,function(a) trajectoire_3(a,h,Te,Ti,r,alpha));
#   return (mean(l)-alpha*sum(x*x)*sum(x*x)/8+alpha/8+Ti);
# }

u0 <- temperature_2(x,h,Ti,n,1)
#ureduc <- temperature_3(x,h,Te,Ti,r,n,1)
# trajectoire_1 <- function(x,h,Te,Ti,r){
#   iter = 0;
#   #print(x);
#   X=x;
#   Y=x;
#   #print(abs(X[1]))
#   
#   o1=c(-0.5,0.5);
#   o2=c(0.5,0.5);
#   o3=c(0.5,-0.5);
#   o4=c(-0.5,-0.5);
#   while (abs(X[1])<=1 && abs(X[2])<=1 && sqrt(sum((X-o1)*(X-o1)))>r && sqrt(sum((X-o2)*(X-o2)))>r && sqrt(sum((X-o3)*(X-o3)))>r && sqrt(sum((X-o4)*(X-o4)))>r ){
#     X=X+sqrt(h)*rnorm(2,0,1);
#     Y=c(Y,X);
#   }
#   if( (abs(X[1])>1) || (abs(X[2])>1)){
#     return (c(Te,Y));
#   }else {
#     return (c(Ti,Y));
#   }
# }
alpha = -1
grille2 <- data.matrix(merge(seq(from = -1, to = 1,by=0.01),seq(from = -1, to = 1,by=0.01)));
#list <-
M2 <-apply(grille2,1,function(a) temperature_2(a,0.01,Ti,50,alpha)[1]);

require(akima)
data <- data.frame(x = grille2[,1],y = grille2[,2],distance = M2)
resolution <- 0.1 # you can increase the resolution by decreasing this number (warning: the resulting dataframe size increase very quickly)
a <- interp(x=data$x, y=data$y, z=data$distance, xo=seq(min(data$x),max(data$x),by=resolution), yo=seq(min(data$y),max(data$y),by=resolution), duplicate="mean")
filled.contour(a, col=palette())
draw.circle(-0.45,0.0,0.65)
