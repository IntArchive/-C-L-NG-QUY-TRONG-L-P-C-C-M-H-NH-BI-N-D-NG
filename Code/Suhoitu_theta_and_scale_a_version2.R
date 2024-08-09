# Load library ggplot2
library(ggplot2)
library(ICSNP)
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
  color = factor(c(rep('a_1,theta_1,v_1', n),
                  rep('a_2,theta_2,v_2', n),
                  rep('a_3,theta_3,v_3', n),
                  rep('a_4,theta_4,v_4', n),
                  rep('a_5,theta_5,v_5', n))))

# Create a scatter plot
ggplot(df, aes(x = x, y = y, color = color)) +  # Initialize ggplot object with data and aesthetic mappings
  geom_line() +  # Add a layer for scatter plot points
  xlab('X') +
  ylab('Y') +
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
###########################################################################################
# Tinh covariance cho v_n
t = n - 2
covariance_v_j = matrix(0,p,p)

matrix_covariance_v = list()
matrix_covariance_v[[1]] <- diag(diag(matrix(1,p,p)))
for(i in 2:t){
  for(j in 1:i){
    #covariance_a_j = covariance_a_j +  (matrix_asim[,j] - a)%*%t((matrix_asim[,j] - a))
    covariance_v_j = covariance_v_j +  (v_n[,j] - v)%*%t((v_n[,j] - v))
    #covariance_v_j = covariance_v_j +  (v_n[,j] - v_n[,i])%*%t((v_n[,j] - v_n[,i]))
  }
  matrix_covariance_v[[i]] <- covariance_v_j
  covariance_v_j = matrix(0,p,p)
}
#########################################################################################
# Da tinh Covariance for v
X <- t(v_n)
n <- nrow(X)
p <- ncol(X)

xbar <- X[2,]
S <- matrix_covariance_v[[1]]
cf = T.ci(mu=xbar, Sigma=S, n=n, avec=c(1,1,1,1,1))
matrix_confint_v = matrix(c(cf[1],cf[2]),ncol=1)


for(i in 2){
  for(j in 2:t){
    xbar <- X[j,]
    S <- (1/j)*matrix_covariance_v[[j]]
    cf = T.ci(mu=xbar, Sigma=S, n=n, avec=c(1,1,0,0,0))
    matrix_confint_v = cbind(matrix_confint_v,matrix(c(cf[1],cf[2]),ncol=1))
  }
  
}


######################################################################################
# dataframe
df <- data.frame(
  x = c(2:t),
  y = c(t(v_n[2,2:t])),
  lowerconfinterval = matrix_confint_v[1,2:t],
  upperconfinterval = matrix_confint_v[2,2:t])

# Create a scatter plot
ggplot() +  # Initialize ggplot object with data and aesthetic mappings
  geom_line(data = df, aes(x = x, y = y), color = 'blue') +  # Add a layer for scatter plot points
  geom_line(data = df, aes(x = x, y = lowerconfinterval),color = 'red') +
  geom_line(data = df, aes(x = x, y = upperconfinterval), color = 'red') +
  geom_hline(yintercept=v[2]) +
  scale_y_continuous(expand = c(0,0), limits = c(-2.5, 2.5))+
  scale_x_continuous(expand = c(0,0), limits = c(0, 2000))+
  xlab('n') +
  ylab('v_2') +
  theme_minimal() +
  scale_color_manual(values = c('blue', 'red', 'green', 'purple', 'green'))





















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






###########################################################################################
# Tinh covariance cho theta
V_t <- function(x,t,Y,a=matrix(c(1,-4,3,-5/2,2),nrow=p,ncol=1),fone=1/2){
  return(diag(c(sign(a*fone)))%*%D(x,t)%*%Y)
} 

hsocovariance_matrix = 4*pi*pi*rma(matrix(abs(a),nrow=1),p, r = TRUE)*rma(abs(a),p, r = FALSE)*abs(f1$value)**2
phi_theta_hat_n = V_t(x[1,1],matrix_theta[,1],Y[,1])%*%t(V_t(x[1,1],matrix_theta[,1],Y[,1]))
matrix_phi_theta_hat = list()
matrix_phi_theta_hat[[1]] <- phi_theta_hat_n
for(i in 2:n){
  phi_theta_hat_n = (1/i)*V_t(x[1,i],matrix_theta[,i],Y[,i])%*%t(V_t(x[1,i],matrix_theta[,i],Y[,i]))
  for(j in 1:i){
    phi_theta_hat_n =phi_theta_hat_n+ (1/i)*V_t(x[1,j],matrix_theta[,i],Y[,j])%*%t(V_t(x[1,j],matrix_theta[,i],Y[,j]))
  }
  matrix_phi_theta_hat[[i]] <- phi_theta_hat_n
}
matrix_covariance_theta <- list()
for(i in 1:n){
  matrix_covariance_theta[[i]] <- hsocovariance_matrix*matrix_phi_theta_hat[[i]]
}

#########################################################################################
# Da tinh lai Covariance for theta
X <- t(matrix_theta)
n <- nrow(X)
p <- ncol(X)

xbar <- X[1,]
S <- matrix_covariance_theta[[1]]
cf = T.ci(mu=xbar, Sigma=S, n=n, avec=c(1,0,0,0,0))
matrix_confint_theta = matrix(c(cf[1],cf[2]),ncol=1)

for(i in 1){
  for(j in 2:n){
    xbar <- X[j,]
    S <- (1/j)*matrix_covariance_theta[[i]]
    cf = T.ci(mu=xbar, Sigma=S, n=n, avec=c(1,0,0,0,0))
    matrix_confint_theta = cbind(matrix_confint_theta,matrix(c(cf[1],cf[2]),ncol=1))
  }
  
}




# dataframe
df <- data.frame(
  x = c(1:n),
  y = c(t(matrix_theta[1,])),
  lowerconfinterval = matrix_confint_theta[1,],
  upperconfinterval = matrix_confint_theta[2,],
  wave = factor(c(rep('theta_1', n),
                  rep('theta_2', n),
                  rep('theta_3', n),
                  rep('theta_4', n),
                  rep('theta_5', n))))

# Create a scatter plot
ggplot() +  # Initialize ggplot object with data and aesthetic mappings
  geom_line(data = df, aes(x = x, y = y),color = 'blue') +  # Add a layer for scatter plot points
  geom_line(data = df, aes(x = x, y = lowerconfinterval), color = 'red') +
  geom_line(data = df, aes(x = x, y = upperconfinterval), color = 'red') +
  #geom_hline(yintercept=theta[1]) +
  scale_y_continuous(expand = c(0,0), limits = c(-2, 2))+
  scale_x_continuous(expand = c(0,0), limits = c(0, 2010))+
  xlab('n') +
  ylab('theta_1') +
  theme_minimal() +
  scale_color_manual(values = c('blue', 'red', 'green', 'purple', 'green'))

























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

# ######################################################################################
# dataframe
df <- data.frame(
  x = c(1:(n-2)),
  y = c(t(matrix_asim[1,])))
#lowerconfinterval = matrix_confint_theta[1,],
#upperconfinterval = matrix_confint_theta[2,],
#  wave = factor(c(rep('a_1', n-2),
#                  rep('a_2', n-2),
#                  rep('a_3', n-2),
#                  rep('a_4', n-2),
#                  rep('a_5', n-2))))

# Create a scatter plot
ggplot() +  # Initialize ggplot object with data and aesthetic mappings
  geom_line(data = df, aes(x = x, y = y, color = 'blue')) +  # Add a layer for scatter plot points
  #geom_line(data = df, aes(x = x, y = lowerconfinterval, color = 'red')) +
  #geom_line(data = df, aes(x = x, y = upperconfinterval, color = 'red')) +
  geom_hline(yintercept=a[1]) +
  scale_y_continuous(expand = c(0,0), limits = c(-25, 5))+
  scale_x_continuous(expand = c(0,0), limits = c(0, 2200))+
  xlab('x') +
  ylab('a_1') +
  theme_minimal() +
  scale_color_manual(values = c('blue', 'red', 'green', 'purple', 'green'))

#########################################################################################

###########################################################################################
#Tinh ma tran covariance cho a

V_t <- function(x,t,Y,a=matrix(c(1,-4,3,-5/2,2),nrow=p,ncol=1),f1=1/2){
  return(diag(c(sign(a*f1)))%*%D(x,t)%*%Y)
} 
# hsocovariance_matrix 
#hsocovariance_matrix = 2*pi*rma(matrix(abs(a),nrow=1),p, r = TRUE) + rma(abs(a),p, r = FALSE)*abs(f1$value) - 1
e_1 =  matrix(0,p,1)
e_1[1] = e_1[1] + 1 # Da chuan bi xong eT_1
M_p = diag(diag(matrix(1,p,p))) - a%*%t(e_1)

# covariance_a_j = (1/(f1$value**2))*cov(C(matrix(x[,1],ncol=1),matrix(matrix_theta[,1],ncol=1)))
# matrix_covariance_a = list()
# matrix_covariance_a[[1]] <- covariance_a_j
# for(i in 1:(n-2)){
#   covariance_a_j = (1/(f1$value**2))*cov(C(matrix(x[,i],ncol=1),matrix(matrix_theta[,i],ncol=1)))
#   matrix_covariance_a[[i]] <- covariance_a_j
# }

###########################################################################################
# Cach thu hai tinh ma tran covariance cho a 
t = n - 2
covariance_a_j = matrix(0,p,p)

matrix_covariance_a = list()
matrix_covariance_a[[1]] <- diag(diag(matrix(1,p,p)))
for(i in 2:t){
  for(j in 1:i){
    #covariance_a_j = covariance_a_j +  (matrix_asim[,j] - a)%*%t((matrix_asim[,j] - a))
    covariance_a_j = covariance_a_j +  (matrix_asim[,j] - matrix_asim[,i])%*%t((matrix_asim[,j] - matrix_asim[,i]))
  }
  matrix_covariance_a[[i]] <- covariance_a_j
  covariance_a_j = matrix(0,p,p)
}


#########################################################################################
# Da tinh Covariance for a
set.seed(123456)
X <- t(matrix_asim)
n <- nrow(X)
p <- ncol(X)

xbar <- X[1,]
S <- matrix_covariance_a[[1]]
cf = T.ci(mu=xbar, Sigma=S, n=n, avec=c(1,0,0,0,0))
matrix_confint_a = matrix(c(cf[1],cf[2]),ncol=1)


for(i in 1){
  for(j in 2:(n)){
    xbar <- X[j,]
    S <- matrix_covariance_a[[i]]
    cf = T.ci(mu=xbar, Sigma=S, n=n, avec=c(1,0,0,0,0))
    matrix_confint_a = cbind(matrix_confint_a,matrix(c(cf[1],cf[2]),ncol=1))
  }
  
}


######################################################################################
# dataframe
df <- data.frame(
  x = c(2:(n)),
  y = c(t(matrix_asim[1,2:n])),
  lowerconfinterval = matrix_confint_a[1,2:n],
  upperconfinterval = matrix_confint_a[2,2:n])
#  wave = factor(c(rep('a_1', n-2),
#                  rep('a_2', n-2),
#                  rep('a_3', n-2),
#                  rep('a_4', n-2),
#                  rep('a_5', n-2))))

# Create a scatter plot
ggplot() +  # Initialize ggplot object with data and aesthetic mappings
  geom_line(data = df, aes(x = x, y = y), color = 'blue') +  # Add a layer for scatter plot points
  geom_line(data = df, aes(x = x, y = lowerconfinterval),color = 'red') +
  geom_line(data = df, aes(x = x, y = upperconfinterval), color = 'red') +
  geom_hline(yintercept=a[1]) +
  scale_y_continuous(expand = c(0,0), limits = c(-2.5, 2.5))+
  scale_x_continuous(expand = c(0,0), limits = c(0, 2000))+
  xlab('n') +
  ylab('a_1') +
  theme_minimal() +
  scale_color_manual(values = c('blue', 'red', 'green', 'purple', 'green'))



















#########################################################################################
# Kiem tra su hoi tu cua a
kiemtra <- function(x, accuracy_x, x_name = 'a', y_name = 'a_hat', name_file_save = "something.pdf"){
  #param  x: matrix of a paramter
  
  dataV <- data.frame(
    x = c(x),
    y = c(accuracy_x),
    x1 = c((diag(matrix(-max(abs(x)),p-1,p-1))),max(abs(x))), 
    y1 = c((diag(matrix(-max(abs(x)),p-1,p-1))),max(abs(x)))
  )
  
  # Create a scatter plot
  image = ggplot() +  # Initialize ggplot object with data and aesthetic mappings
    geom_point(data = dataV, aes(x = x, y = y), color = "red") +  # Add a layer for scatter plot points
    geom_line(data = dataV, aes(x = x1, y = y1),color = "blue") +
    xlab(x_name) +
    ylab(y_name) +
    theme_minimal()
  ggsave(file=paste("D:/Minh_Hai_D/10_HCMUE_Semester7/1_KLTN/image1/",name_file_save), plot=image, width=11.33, height=7)
}

kiemtra(theta,matrix_theta[,n], x_name = 'theta', y_name = 'theta_hat', name_file_save = "Approx_theta_with_thetahat2000.pdf")
kiemtra(a, matrix_asim[,n-2])
kiemtra(v, v_n[,n], x_name = 'v', y_name = 'v_hat')


