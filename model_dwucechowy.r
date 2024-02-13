library(AGHmatrix)
library(MASS)
library(data.table)

# Model 2-cechowy

#### MACIERZ G i R ####

G <- matrix(c(235.55,-238.38, -283.38, 412.24), 2, 2)
R <- matrix(c(3632.90, -1107.40, -1107.40, 8134.30), 2, 2)
dim(G)
G

#### WCZYTANIE I OBRÓBKA DANYCH Z POKREWIEŃSTWA I EFEKTÓW STAŁYCH ####

ped <- read.csv('http://theta.edu.pl/wp-content/uploads/2023/01/pedigree.csv',
                sep = ';', header=FALSE)

colnames(ped) <- c('id', 'sire', 'dam', 'birthY', 'breed')

pednew <- c(rbind(ped$dam, ped$sire, ped$id))
pednew <- data.frame(pednew)
pednew <- unique(pednew)
colnames(pednew) <- 'id'
pednew <- merge(pednew, ped, by='id', all.x=TRUE)

pednew[is.na(pednew)] = 0
pednew <- pednew[order(pednew$birthY),]

pednew$idNew <- 1:nrow(pednew)

#### WCZYTANIE I OBRÓBKA DANYCH Z INFO O LAKTACJI ####

# Przygotowanie danych z fenotypem
lak <- read.table('http://theta.edu.pl/wp-content/uploads/2023/01/laktacje.csv',
                  sep=';', header=TRUE)

lak$caldate <- as.Date(lak$calvingDate, format='%d.%m.%Y')
lak$insdate <- as.Date(lak$inseminationDate, format='%d.%m.%Y')

lak <- as.data.table(lak)
date_icf <- lak[, .(icfDate = unique(caldate)), by = cowId ]
minDate_ins <- lak[, .(minDate = min(insdate)), by = cowId ]
maxDate_ins <- lak[, .(maxDate = max(insdate)), by = cowId ]

ICF <- merge(date_icf, minDate_ins, by='cowId', all.x=TRUE)
ICF$ICF <- as.integer(ICF$icfDate - ICF$minDate)
IFL <- merge(minDate_ins, maxDate_ins, by='cowId', all.x=TRUE)
IFL$IFL <- as.integer(IFL$maxDate - IFL$minDate)

cechy <- merge(IFL[, c(1,4)], ICF[,c(1,4)], by='cowId', all.x=TRUE)

cechy <- merge(cechy, pednew[, c(1,4:6)], by.x='cowId', by.y='id', all.x=TRUE)
cechy <- data.frame(cechy)
cechy <- cechy[order(cechy$idNew),]

rm(date_icf, minDate_ins, maxDate_ins, ICF, IFL)

#### MACIERZ SPOKREWNIEŃ ####

A <- Amatrix(pednew[,1:3])
dim(A)

#### MACIERZ CECH ####

y <- as.matrix(c(cechy$IFL, cechy$ICF))
(dim(y))

#### MACIERZ EFEKTÓW STAŁYCH ####

X1 <- model.matrix(~ factor(cechy$birthY) - 1)
X2 <- model.matrix(~ factor(cechy$breed) - 1)
x1 <- cbind(X1, X2)
x2 <- x1

w <- nrow(y)
k <- ncol(x1)*2

X <- matrix(0, w, k)
X[1:179, 1:19] <- x1
X[180:w, 20:k] <- x2
dim(X)

rm(X1, x1, X2, x2)

#### MACIERZ Z ####

w <- nrow(cechy)
k <- nrow(pednew)

Z1 <- matrix(0, w, k)
I <- diag(w)
Z1[1:w, 321:499] = I
Z2 <- Z1

Z <- matrix(0, w*2, k*2)
Z[1:w, 1:k] <- Z1
Z[(w+1):(w*2), (k+1):(k*2)] <- Z2

rm(Z1, Z2)
dim(Z)

#### FUNKCJA MME2 ####

mme2 = function(y, X, Z, A, G, R) {
    invA = ginv(A)
    invG = ginv(G)
    R = kronecker(R, diag(179))
    invR = ginv(R)
    C = rbind(cbind(t(X)%*%invR%*%X, t(X)%*%invR%*%Z),
              cbind(t(Z)%*%invR%*%X, t(Z)%*%invR%*%Z+kronecker(invG, invA)))
    rhs = rbind(t(X)%*%invR%*%y, t(Z)%*%invR%*%y)
    invC = ginv(C)
    estimators = invC%*%rhs
    list(C = C, est = round(estimators,4))
}

wyniki <- mme2(y, X, Z, A, G, R)

# IFL
year.IFL2 <- wyniki$est[1:13]
names(year.IFL2) <- sort(unique(ped$birthY))
breed.IFL2 <- wyniki$est[14:19]
names(breed.IFL2) <- paste('breed', 1:6)
est.IFL2 <- wyniki$est[39:537]
names(est.IFL2) <- pednew$id

# ICF
year.ICF2 <- wyniki$est[20:32]
names(year.ICF2) <- sort(unique(ped$birthY))
breed.ICF2 <- wyniki$est[33:38]
names(breed.ICF2) <- paste('breed', 1:6)
est.ICF2 <- wyniki$est[538:length(wyniki$est)]
names(est.ICF2) <- pednew$id

























