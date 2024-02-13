# Funkcja do rozwiązywania układów równań mieszanych
mme = function(y, X, Z, A, sigma_a, sigma_e) {
    alpha = sigma_e / sigma_a
    invA = ginv(A)
    C = rbind(cbind(t(X)%*%X, t(X)%*%Z),
              cbind(t(Z)%*%X, t(Z)%*%Z+invA*c(alpha)))
    rhs = rbind(t(X)%*%y, t(Z)%*%y)
    invC = ginv(C)
    estimators = invC%*%rhs
    list(mC = C, est = estimators)
}

# Estymacja paramterów wariancji

sigma_a = 1.01  #starting value for random effect
sigma_e = 10.01  #starting value for error variance

EM = function(y, X, Z, A, sigma_a, sigma_e) {
    n = nrow(X)
    p = ncol(X) 
    q = nrow(A) 
    
    t = 1 #iteration number 1
    tmp = 0.1 #test for convergance
    
    while (tmp > 0.0001) {
        mme_new = mme(y, X, Z, A, sigma_a, sigma_e)
        C_new = ginv(mme_new$mC)
        Ck = C_new[(p+1):(p+q), (p+1):(p+q)]
        mme2 = mme_new$est
        
        a = as.matrix(mme2[(p+1):(p+q)])
        sigma_a_new = (t(a)%*%ginv(A)%*%a + sum(diag(ginv(A)%*%Ck))*c(sigma_e))/q
        
        res = as.matrix(y-X%*%as.matrix(mme2[1:p]) - Z%*%as.matrix(mme2[(p+1):(p+q)]))
        X.tmp1 = cbind(X,Z) %*% C_new
        X.tmp2 = t(cbind(X,Z))
        sigma_e_new = (t(res)%*%res + sum(diag(X.tmp1%*%X.tmp2))*c(sigma_e))/n
        
        tmp = max(abs(sigma_a - sigma_a_new), abs(sigma_e - sigma_e_new))
        sigma_a = sigma_a_new
        sigma_e = sigma_e_new
        
        t = t + 1
    }
    test = list(t = t, sigma_a = sigma_a, sigma_e = sigma_e)
    return(test)
}

# Obliczamy wariancję
start <- Sys.time()
wyniki <- EM(y, X, Z, A, sigma_a, sigma_e)
stop <- Sys.time()
stop - start

# Ocena (to samo jest w pliku ocena_modelu)#

# Sprawdzamy dokładność wartości hodowlanych
C <- as.matrix(wynik$mC)
invC <- ginv(C)
k <- ncol(X)
n <- nrow(C)

invC22 <- invC[(k+1):n, (k+1):n]
alpha <- wyniki$sigma_e / wyniki$sigma_a
r2 <- diag(1 - invC22*alpha)
r2[r2 < 0] <- 0
r <- sqrt(r2)
summary(r)

# Odziedziczalność
h2 <- sigma_a / (wyniki$sigma_a + wyniki$sigma_e)
h2

# Test Walda
G <- A*wyniki$sigma_a
R <- diag(nrow(IFL))*wyniki$sigma_e
V <- Z%*%G%*%t(Z) + R

varB <- ginv(t(X)%*%ginv(V)%*%X)
seB <- sqrt(diag(varB))

testWalda <- wynik$est[1:k] / seB
p_value <- 2*pnorm(abs(testWalda), lower.tail = FALSE)

which(p_value > 0.05)