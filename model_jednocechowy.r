library(AGHmatrix)
library(MASS)
library(data.table)
library(pheatmap)
install.packages('AGHmatrix')
install.packages('MASS')
install.packages('data.table')
install.packages('pheatmap')


#### Przygotowanie danych pedigree ####

ped <- read.table('http://theta.edu.pl/wp-content/uploads/2023/01/pedigree.csv',
                sep=';', header=FALSE)
colnames(ped) <- c('id', 'sire', 'dam', 'birthY', 'breed')

pednew <- c(rbind(ped$dam, ped$sire, ped$id))
pednew <- data.frame(pednew)
pednew <- unique(pednew)
colnames(pednew) <- 'id'
pednew <- merge(pednew, ped, by='id', all.x=TRUE)

pednew[is.na(pednew)] = 0
pednew <- pednew[order(pednew$birthY),]

pednew$idNew <- 1:nrow(pednew)

#### Macierz spokrewnień ####

A <- Amatrix(pednew[,1:3])
dim(A)

library(pheatmap)

# Utworzenie heatmap dla macierzy A
pheatmap(A, main="Heatmap of Matrix A", cluster_rows = TRUE, cluster_cols = TRUE)

# Przygotowanie danych z fenotypem
lak <- read.table('http://theta.edu.pl/wp-content/uploads/2023/01/laktacje.csv',
                sep=';', header=TRUE)

lak$caldate <- as.Date(lak$calvingDate, format='%d.%m.%Y')
lak$insdate <- as.Date(lak$inseminationDate, format='%d.%m.%Y')

lak <- as.data.table(lak)
minDate <- lak[, .(minDate = min(insdate)), by = cowId ]
maxDate <- lak[, .(maxDate = max(insdate)), by = cowId ]
#IFL
IFL <- merge(minDate, maxDate, by='cowId', all.x=TRUE)
IFL$IFL <- as.integer(IFL$maxDate - IFL$minDate)
IFL <- IFL[, c(1, 4)]
IFL <- merge(IFL, pednew[, c(1,4:6)], by.x='cowId', by.y='id', all.x=TRUE)
IFL <- data.frame(IFL)
IFL <- IFL[order(IFL$idNew),]

summary(IFL$IFL)

#ICF
date_icf <- lak[, .(icfDate = unique(caldate)), by = cowId ]
ICF <- merge(date_icf, minDate, by='cowId', all.x=TRUE)
ICF$ICF <- as.integer(ICF$icfDate - ICF$minDate)
ICF <- ICF[, c(1,4)]
ICF <- merge(ICF, pednew[, c(1,4:6)], by.x='cowId', by.y='id', all.x=TRUE)
ICF <- data.frame(ICF)
ICF <- ICF[order(ICF$idNew),]
summary(ICF$ICF) #( 270.0, 768.0)

rm(maxDate, minDate, date_icf)

#### Macierz efektów stałych ####

X1 <- model.matrix(~ factor(IFL$birthY) - 1)
X2 <- model.matrix(~ factor(IFL$breed) - 1)
X <- cbind(X1, X2)
dim(X)

rm(X1, X2)


#### Macierz Z ####
Z <- matrix(0, nrow(IFL), nrow(pednew))
I <- diag(179)
Z[1:179, 321:499] = I
dim(Z)
rm(I)

#### zmienna y ####
y <- as.matrix(IFL$IFL)

#### Funkcja do rozwiązywania układów równań mieszanych ####
mme = function(y, X, Z, A, sigma_a, sigma_e) {
    alpha = sigma_e / sigma_a
    invA = ginv(A)
    C = rbind(cbind(t(X)%*%X, t(X)%*%Z),
              cbind(t(Z)%*%X, t(Z)%*%Z+invA*c(alpha)))
    rhs = rbind(t(X)%*%y, t(Z)%*%y)
    invC = ginv(C)
    estimators = invC%*%rhs
    list(mC = C, est = round(estimators,4))
}


# IFL
sigma_a <- 235.55
sigma_e <- 3632.90

res.IFL <- mme(y, X, Z, A, sigma_a, sigma_e)
year.IFL <- res.IFL$est[1:13]
names(year.IFL) <- sort(unique(ped$birthY))
breed.IFL <- res.IFL$est[14:19]
names(breed.IFL) <- paste('breed', 1:6)
est.IFL <- res.IFL$est[20:length(res.IFL$est)]
names(est.IFL) <- pednew$id

summary(year.IFL)
year.IFL

# ICF
sigma_a <- 412.24
sigma_e <- 8134.30

res.ICF <- mme(ICF$ICF, X, Z, A, sigma_a, sigma_e)
year.ICF <- res.ICF$est[1:13]
names(year.ICF) <- sort(unique(ped$birthY))
breed.ICF <- res.ICF$est[14:19]
names(breed.ICF) <- paste('breed', 1:6)
est.ICF <- res.ICF$est[20:length(res.ICF$est)]
names(est.ICF) <- pednew$id



