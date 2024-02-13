genotypes <- read.csv('http://theta.edu.pl/wp-content/uploads/2024/01/Genotypes.csv', 
                      sep = ';', header = FALSE)

p <- numeric(ncol(genotypes)-1)

tmp <- data.frame(matrix(0, 3, 1))
colnames(tmp) <- c("Var1")
rownames(tmp) <- c("AA", "AB", "BB")
tmp[, 1] <- c(0, 1, 2)

head(genotypes)
head(tmp)

for (i in 2:length(p)) {
    tmp2 <- merge(tmp, as.data.frame(table(genotypes[, i])), by = "Var1", all.x = TRUE)
    tmp2[is.na(tmp2)] <- 0
    p[i-1] <- (tmp2$Freq[1]*2 + tmp2$Freq[2]/2) / (2*sum(tmp2$Freq))
}

k <- which(p < 0.05)
k2 <- which(p > 0.95)

if (length(k) != 0 & length(k2) != 0) {
    p <- p[-c(k, k2)]
    genotypes2 <- genotypes[, -c(k, k2)]
} else if (length(k) != 0 & length(k2) == 0) {
    p <- p[-k]
    genotypes2 <- genotypes[, -k]
} else if (length(k) == 0 & length(k2) != 0) {
    p <- p[-k2]
    genotypes2 <- genotypes[, -k2]
} else {
    p <- p
    genotypes2 <- genotypes
}

genotypes2 <- genotypes2[,-1]

M <- matrix(0, nrow(genotypes), length(p))

for (i in 1:nrow(genotypes)) {
    for (j in 1:length(p)) {
        if (genotypes2[i, j] == 0) {
            M[i, j] <- 2 - 2*p[j]
        } else if (genotypes2[i, j] == 1) {
            M[i, j] <- 1 - 2*p[j]
        } else {
            M[i, j] <- - 2*p[j]
        }
    }
}

G <- M%*%t(M) / (2*sum(p*(1-p)))


# Funkcja do rozwiązywania układów równań mieszanych z genotypami
mme_genotype = function(y, X, Z1, Z2, A, G, sigma_a, sigma_g, sigma_e) {
    alpha1 = sigma_e / sigma_a
    alpha2 = sigma_e / sigma_g
    invA = ginv(A)
    invG = ginv(G)
    C = rbind(cbind(t(X)%*%X, t(X)%*%Z1, t(X)%*%Z2),
              cbind(t(Z1)%*%X, t(Z1)%*%Z1+invA*c(alpha1), t(Z1)%*%Z2),
              cbind(t(Z2)%*%X, t(Z2)%*%Z1, t(Z2)%*%Z2 + invG*c(alpha2)))
    rhs = rbind(t(X)%*%y, t(Z1)%*%y, t(Z2)%*%y)
    invC = ginv(C)
    estimators = invC%*%rhs
    list(C = C, est = estimators)
}

sigma_a <- 235.55
sigma_e <- 3632.90

# niby taki wzór ale się nie pokrywa z tym co Suchy ma
sigma_g <- sigma_a / length(p)

# Sprawdzanie czy efekty markerów genetycznych SNP są statystycznie istotne

k = 0.2
snp = ginv((t(M)%*%ginv(A)[321:499, 321:499]%*%M)/k + 1 / ((1-k)/length(p2)) ) / k 
snp = snp%*%t(M)%*%ginv(A)[321:499, 321:499]
snp = snp %*% est[8:12]

(s = sd(snp))
(m = mean(snp))

snp2 = snp / s
sd(snp2)

W = snp2
p.value = 2*pnorm(abs(W), lower.tail = FALSE)

length(which(p.value < 0.05))


