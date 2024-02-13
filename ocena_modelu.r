#### OCENA MODELU JEDNOCECHOWEGO ####
library(matlib)
install.packages('matlib')
#### IFL ####

sigma_a <- 235.55
sigma_e <- 3632.90

# Sprawdzamy dokładność wartości hodowlanych
C <- as.matrix(res.IFL$mC)
invC <- ginv(C)
k <- ncol(X)
n <- nrow(C)

invC22 <- invC[(k+1):n, (k+1):n]
alpha <- sigma_e / sigma_a
r2 <- diag(1 - invC22*alpha)
r2[r2 < 0] <- 0
r <- sqrt(r2)

summary(r)
hist(r)

# Odziedziczalność
h2 <- sigma_a / (sigma_a + sigma_e)
h2

# Test Walda
G <- A*sigma_a
R <- diag(nrow(IFL))*sigma_e
V <- Z%*%G%*%t(Z) + R

varB <- ginv(t(X)%*%ginv(V)%*%X)
seB <- sqrt(diag(varB))

testWalda <- res.IFL$est[1:k] / seB
p_value <- 2*pnorm(abs(testWalda), lower.tail = FALSE)

which(p_value > 0.05)

#### ICF ####

sigma_a <- 412.24
sigma_e <- 8134.30

# Sprawdzamy dokładność wartości hodowlanych
C <- as.matrix(res.ICF$mC)
invC <- ginv(C)
k <- ncol(X)
n <- nrow(C)

invC22 <- invC[(k+1):n, (k+1):n]
alpha <- sigma_e / sigma_a
r2 <- diag(1 - invC22*alpha)
r2[r2 < 0] <- 0
r <- sqrt(r2)

summary(r)
hist(r)

# Odziedziczalność
h2 <- sigma_a / (sigma_a + sigma_e)
h2

# Test Walda
G <- A*sigma_a
R <- diag(nrow(ICF))*sigma_e
V <- Z%*%G%*%t(Z) + R

varB <- ginv(t(X)%*%ginv(V)%*%X)
seB <- sqrt(diag(varB))

testWalda <- res.ICF$est[1:k] / seB
p_value <- 2*pnorm(abs(testWalda), lower.tail = FALSE)

which(p_value > 0.05)


#### OCENA MODELU DWUCECHOWEGO ####

