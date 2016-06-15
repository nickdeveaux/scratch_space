source('mi_and_clr.R')

X <- read.table('X.csv', sep = ',', header = 1, row.names = 1)
Y <- read.table('Y.csv', sep = ',', header = 1, row.names = 1)

PARS <- list()
PARS$mi.bins <- 10
PARS$cores <- 10

# fill mutual information matrices
Ms <- mi(t(Y), t(X), nbins=PARS$mi.bins, cpu.n=PARS$cores)
diag(Ms) <- 0
Ms_bg <- mi(t(X), t(X), nbins=PARS$mi.bins, cpu.n=PARS$cores)
diag(Ms_bg) <- 0

# get CLR matrix
clr.mat = mixedCLR(Ms_bg,Ms)
dimnames(clr.mat) <- list(rownames(Y), rownames(X))

write.table(clr.mat, 'clr_matrix.tsv', sep = '\t')