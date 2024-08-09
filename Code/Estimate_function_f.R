##########
### Estimate regression function f
# Load library ggplot2
library(ggplot2)
library(ICSNP)
library(tidyverse)
set.seed(45)  # Thiet lap seed de tai tao ket qua
T.test <- function(X, mu=0){
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  df2 <- n - p
  if(df2 < 1L) stop("Need nrow(X) > ncol(X).")
  if(length(mu) != p) mu <- rep(mu[1], p)
  xbar <- colMeans(X)
  S <- cov(X)
  T2 <- n * t(xbar - mu) %*% solve(S) %*% (xbar - mu)
  Fstat <- T2 / (p * (n-1) / df2)
  pval <- 1 - pf(Fstat, df1=p, df2=df2)
  data.frame(T2=as.numeric(T2), Fstat=as.numeric(Fstat),
             df1=p, df2=df2, p.value=as.numeric(pval), row.names="")
}


T.ci <- function(mu, Sigma, n, avec=rep(1,length(mu)), level=0.95){
  p <- length(mu)
  if(nrow(Sigma)!=p) stop("Need length(mu) == nrow(Sigma).")
  if(ncol(Sigma)!=p) stop("Need length(mu) == ncol(Sigma).")
  if(length(avec)!=p) stop("Need length(mu) == length(avec).")
  if(level <=0 | level >= 1) stop("Need 0 < level < 1.")
  cval <- qf(level, p, n-p) * p * (n-1) / (n-p)
  zhat <- crossprod(avec, mu)
  zvar <- crossprod(avec, Sigma %*% avec) / n
  const <- sqrt(cval * zvar)
  c(lower = zhat - const, upper = zhat + const)
}

rma <- function(x, p, r = TRUE){
  # Extend matrix follow row or column
  # matrix                  x: the input matrix
  # line of extend          p: the number of line that extend procedure do
  # direction of extension  r: TRUE if extend by row, FALSE if extend by column
  # state: checked
  c_or_r = x
  if(r == TRUE){
    for(i in 2:p){
      x = rbind(x,c_or_r)
    }
  }else{
    for(i in 2:p){
      x = cbind(x,c_or_r)
    }
  }
  return(x)
}

sum_broadcast <- function(a,b){
  # a and b is corresponding matrix dim (p,1) and (1,p)
  return(rma(a,dim(b)[2],r=FALSE)+rma(b,dim(a)[1]))
}

D <- function(X,t,p=5){
  # state: checked
  # input X: a random variable, EX: X[1,1]
  g = 1
  t = matrix(t,nrow=length(t), ncol=1)
  s = (1/g)*c(sin(2*pi*(X-t[1,])))
  for(i in 2:p){
    s = c(s,(1/g)*sin(2*pi*(X-t[i,])))
  }
  return(diag(s))
}

my_f <- function(x) {    # Create own function
  out <- cos(2*pi*(x))*(cos(2*pi*1*(x)) + cos(2*pi*2*(x)) + cos(2*pi*3*(x)) + cos(2*pi*4*(x)) + cos(2*pi*5*(x)))
  return(out)
}
matrix_pika <- function(x_bef, x_aft){
  # state: checked
  
  x_bef <- matrix(x_bef,ncol = 1)
  x_aft <- matrix(x_aft,ncol = 1)
  x_bef <- (abs(matrix(pi_K_vectorized(x_aft), ncol = 1)) < 0.25)*x_aft + (abs(matrix(pi_K_vectorized(x_aft), ncol = 1)) >= 0.25)*x_bef
  return(x_bef)
}
f1 <- integrate(my_f,          # Apply integrate in R
                lower = -1/2,
                upper = 1/2)

f <- function(x){
  out2 <- cos(2*pi*1*(x))+cos(2*pi*2*(x))+cos(2*pi*3*(x))+cos(2*pi*4*(x))+cos(2*pi*5*(x))
  return(out2)
}

Vectorize_f <- Vectorize(f)


p =  5
n = 2000
g = 1
rowx <- runif(n, min=-0.5,max= 0.5)  # Tao gia tri x
x = matrix(rowx,nrow = 1,ncol = length(rowx))
x = rma(x,p=p)

# Chuan bi cac vecto v, theta, a and esp
v = matrix(c(0,1/3,-1,2,-9/10),nrow=p,ncol=1)
theta = matrix(c(0,1/5,-1/20,-1/7,1/6),nrow=p,ncol=1)
a = matrix(c(1,-4,3,-5/2,2),nrow=p,ncol=1)


##
f_vec <- Vectorize(f)
gamma_j <-  1/(2*pi*a*(1/2))
##

esp = matrix(rnorm(n*p, mean = 0, sd = 1), nrow=p)
Y = (rma(a, n, r= FALSE))*f_vec(x-rma(theta,n,r=FALSE)) + rma(v,n, r=FALSE) + esp

# Khi ma bind xong thi Y la matrix

# dataframe
df <- data.frame(
  x = c(t(x)),
  y = c(t(Y)),
  wave = factor(c(rep('a_1,theta_1,v_1', n),
                  rep('a_2,theta_2,v_2', n),
                  rep('a_3,theta_3,v_3', n),
                  rep('a_4,theta_4,v_4', n),
                  rep('a_5,theta_5,v_5', n))))

# Create a scatter plot
ggplot(df, aes(x = x, y = y, color = wave)) +  # Initialize ggplot object with data and aesthetic mappings
  geom_line() +  # Add a layer for scatter plot points
  xlab('x') +
  ylab('y') +
  theme_minimal() +
  scale_color_manual(values = c('blue', 'red', 'green', 'purple', 'yellow'))





##################################################
# Estimate transition v parameters
# Calculating v_n for each type
v_n = matrix(diag(diag(p))-0.75, nrow=p, ncol=1)
cv = 0
for(i in 2:n){
  cv <- 0
  for(j in 1:i){
    cv <- cv + (1/i)*Y[,j]/g  # This calculates average per row
  }
  v_n <- cbind(v_n, (cv))  
}






##################################
# Estimate theta parameters
f1_hat <- function(x,Y,n, g = 1){
  # x is vector EX: x[1,]
  # y is vector EX: Y[1,]
  x = matrix(x,nrow = 1)
  Y = matrix(Y,ncol = 1)
  return((1/n)*((cos(2*pi*x)/g)%*%Y))
}

# f1_hat = (1/n)*matrix(cos(2*pi*x[1,]),nrow = 1)%*%matrix(Y[1,],nrow = length(Y[1,]))  

C <- function(X,t,p=5){
  g = 1
  t = matrix(t,nrow=length(t), ncol=1)
  s = (1/g)*c(cos(2*pi*(X-t[1,])))
  for(i in 2:p){
    s = c(s,(1/g)*cos(2*pi*(X-t[i,])))
  }
  return(diag(s))
}





# Estimate theta
pi_K <- function(x) {
  #state: checked
  ifelse(abs(x) <= 1/4, x, ifelse(x > 1/4, 1/4, -1/4))
}
pi_K_vectorized <- Vectorize(pi_K) # Tra ve mot day so co tri tuyet doi nho hon hoac bang 0.25

N = n
gamma <- 1/ (1:n)  # Step size, gamma_n = 1/n
theta_hat = c(0,0,0,0,0)
matrix_theta =  matrix(theta_hat, ncol=1)

# Estimate theta
for (nn in 1:(N-1)) { #(N-1)
  T_np <- D(x[1,nn+1],theta_hat)%*%matrix(Y[,nn+1],nrow=p,ncol=1) # Assuming Y_n = 1 and g(X_n) = 1 (uniform density)
  index = nn + 1
  #theta_hat <- matrix(theta_hat,nrow=p,ncol=1) + sign(a*f1$value)*(gamma_j*gamma[index]) * T_np
  theta_hat <- pi_K_vectorized(matrix(theta_hat,nrow=p,ncol=1) + (gamma_j*gamma[index]) * T_np)
  matrix_theta = cbind(matrix_theta, theta_hat)
}















#########################################################################################
# Estimate asim
C <- function(X,t,p=5){
  g = 1
  t = matrix(t,nrow=length(t), ncol=1)
  s = (1/g)*c(cos(2*pi*(X[1,]-t[1,])))
  if(p >= 2){
    for(i in 2:p){
      s = c(s,(1/g)*cos(2*pi*(X[i,]-t[i,])))
    }  
  }
  return(diag(s))
}

vector_asim = (1/(2*f1$value))*diag(C(matrix(x[1,2:2],ncol=1), matrix_theta[1,1:(2-1)],p=2-1))%*%matrix(Y[1,2:2],ncol = 1)
matrix_asim = matrix(0, nrow = 1, ncol = n - 2)
for(j in 1:p){
  for(i in 3:n){
    asim_i = (1/(i*f1$value)*diag(C(matrix(x[1,2:i],ncol=1), matrix_theta[j,1:(i-1)],p=i-1))%*%matrix(Y[j,2:i],ncol = 1))
    vector_asim = c(vector_asim, asim_i)
  }
  matrix_asim <- rbind(matrix_asim, matrix(vector_asim,nrow = 1))
  vector_asim = (1/(2*f1$value))*diag(C(matrix(x[1,2:2],ncol=1), matrix_theta[j,1:(2-1)],p=2-1))%*%matrix(Y[j,2:2],ncol = 1)
}
matrix_asim = matrix_asim[2:dim(matrix_asim)[1], ]










last_col <- matrix(matrix_asim[,n-(n - dim(matrix_asim)[2] + 1)], ncol = 1)
for(i in 1:(n - dim(matrix_asim)[2])){
  matrix_asim <- cbind(matrix_asim, last_col)
}




###########################################################################
# Da chuan bi xong v_n, matrix_theta, matrix_asim
fn_hat <- function(x, kernel, n = 1500){
  weight = diag(1/p*diag(1,p,p))
  alpha = 9/10
  x = matrix(x,nrow = 1, ncol = n)
  
  
  
  # uniform_kernel_vectorized return sequence c()
  h = matrix(c(1:n), nrow = 1)**alpha
  xn = matrix(x[,1:n], nrow = 1)
  for(j in 1:p){
    thetan_j = matrix(matrix_theta[j,1:n], nrow = 1)
    xn = xn - thetan_j - x
    inpK = (1/h)*xn
  }  
}

uniform_kernel <- function(u) {
  # Applying the uniform kernel function
  # K(u) = 1/2 if |u| <= 1, otherwise K(u) = 0
  ifelse(abs(u) <= 2, 1/4, 0)
}

uniform_kernel_vectorized <- Vectorize(uniform_kernel)

gaussian_kernel <- function(u) {
  # Gaussian kernel function
  1 / sqrt(2 * pi) * exp(-u^2 / 2)
}

gaussian_kernel_vectorized <- Vectorize(gaussian_kernel)

inpK <- function(xn, theta, x, h){
  # Kiem soat input dau vao
  xn = matrix(xn, nrow = 1)
  theta = matrix(theta, nrow = 1)
  x = matrix(x, nrow = 1)
  h = matrix(h, nrow = 1)
  
  xn = xn - theta - x
  inpK = (1/h)*xn
  return(inpK)
}
# uniform_kernel_vectorized return sequence c()
weight = matrix(diag(1/p*diag(1,p,p)), nrow = 1)
alpha = 8/10
matrix_fn_hat = seq()
nopoint = 2000
for(i in 3:nopoint){
  h = (1/matrix(c(1:n), nrow = 1, ncol = n))**alpha
  xn = matrix(x[1,1:n], nrow = 1)
  x1 = matrix(seq(-0.5,0.5, length.out = nopoint)[i],nrow = 1, ncol = n) 
  matrix_fn_hat_j = seq()
  for(j in 1:p){
    thetan_j = matrix(c(matrix_theta[j,1],matrix_theta[j,1:n-1]), nrow = 1, ncol = n)
    inpKp = inpK(xn, thetan_j,  x1, h)
    inpKm = inpK(xn, thetan_j, -x1, h)
    matrixKp = matrix(uniform_kernel_vectorized(inpKp),nrow = 1)
    matrixKm = matrix(uniform_kernel_vectorized(inpKm),nrow = 1)
    Swp = (1/h)*matrixKp
    Swm = (1/h)*matrixKm
    Sw = (Swp + Swm)
    
    fsim_j = matrix(Y[j,1:n],ncol = 1) - matrix(c(v_n[j,1],v_n[j,1:n-1]), ncol = 1)
    
    fn_hat_j = (1/matrix_asim[j,n])*(Sw%*%fsim_j)/sum(Sw)
    matrix_fn_hat_j = c(matrix_fn_hat_j, fn_hat_j[1])
  }
  fn_hat = weight%*%matrix(matrix_fn_hat_j[2:length(matrix_fn_hat_j)], ncol = 1, nrow = length(matrix_fn_hat_j[2:length(matrix_fn_hat_j)]))
  matrix_fn_hat = c(matrix_fn_hat, fn_hat)
}
  
######################################################################################
# dataframe
df <- tibble(
  x = c(seq(-0.5,0.5, length.out=nopoint-1)),
  f_hat_n = c(matrix(matrix_fn_hat[2:n],ncol = 1, nrow = n-1)),
  f = Vectorize_f(c(seq(-0.5,0.5, length.out=nopoint-1))))

df %>%
  pivot_longer(f_hat_n:f, names_to = "Function") %>%
# Create a scatter plot
ggplot(aes(x = x, y = value, color = `Function`)) +  # Initialize ggplot object with data and aesthetic mappings
  geom_line(size=1)+
  #geom_line(data = df, aes(x = x, y = lowerconfinterval, color = 'red')) +
  #geom_line(data = df, aes(x = x, y = upperconfinterval, color = 'red')) +
  #geom_hline(yintercept=a[1]) +
  #scale_y_continuous(expand = c(0,0), limits = c(-0, 2))+
  #scale_x_continuous(expand = c(0,0), limits = c(0, 2200))+
  xlab('x') +
  ylab('value') +
  theme_minimal() +
  scale_color_manual(values = c('red','blue'), labels = c("f","f_hat_n"))
