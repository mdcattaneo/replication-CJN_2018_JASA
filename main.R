#########################################################################################
## R CODE FOR CATTANEO-JANSSON-NEWEY (2016)
## Inference in Linear Regression Models with Many Covariates and Heteroskedasticity
## DATE: 20-Dec-2016
## MAIN SIMULATION CODE
#########################################################################################
# setwd("/mainfiles/cattaneo/Dropbox/research/HCSEManyCov/simulations")
# cd /mainfiles/cattaneo/Dropbox/research/HCSEManyCov/simulations; R
# setwd("Z:/research/HCSEManyCov/simulations")
#########################################################################################
# RUN LOCALLY: R CMD BATCH --vanilla main.R main.Rout &
#########################################################################################
# rm(list=ls(all=TRUE))

source("main.fun.R")

#########################################################################################
## SEED SETUP FOR PARALLELIZATION
#########################################################################################
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) { args <- list("1") }
message("The argument is ", args[[1]])

# Set seed
# set.seed(666)
set.seed(as.integer(args[[1]])*666)

#########################################################################################
## MONTECARLO SETUP
#########################################################################################
## Number of simulations (S)
S = 30
## Sample size (n)
n = 700
## Bootstrap replications (B)
B = 500
## d = dim(beta)
d = 1; beta = 1
## Nuisance covariates
K.grid = 1:5;
K.grid.LM = floor(n*c(0,.1,.2,.3,.4));
K.grid.PL = c(0,1.5,2.5,3.5,4.5);
K.grid.FE.T = c(n,10,5,4,3);
K.grid.FE.N = floor(n/K.grid.FE.T);
## Load population scaling factors
popval = read.csv("main.popval.csv", row.names=1)
## Models to run
#models = c(1,4,7,10,13)
models = 1:15

#########################################################################################
# Run Monte Carlo Experiment
#########################################################################################
# m = 1
for (m in models){
## Output Tables
## HO0 = HOM/n | HO1 = HOM/(n-K-d)
## HC0 = HET/n | HC1 = HET/(n-K-d) | HC2 = HET/Mii | HC3 = HET/Mii²
## HCK = HET w/ W=ginv(M*M)
col.names1 = c("n","d","beta","m","vartheta","s","K","rank")
col.names2 = c("beta.hat","V.HO0","V.HO1","V.HC0","V.HC1","V.HC2","V.HC3","V.HC4","V.HCK",
                          "T.HO0","T.HO1","T.HC0","T.HC1","T.HC2","T.HC3","T.HC4","T.HCK")
col.names3 = c("V.hat.B","T.HO0.B.q025","T.HO1.B.q025","T.HC0.B.q025","T.HC1.B.q025","T.HC2.B.q025","T.HC3.B.q025","T.HC4.B.q025","T.HCK.B.q025",
                         "T.HO0.B.q975","T.HO1.B.q975","T.HC0.B.q975","T.HC1.B.q975","T.HC2.B.q975","T.HC3.B.q975","T.HC4.B.q975","T.HCK.B.q975")
col.names=c(col.names1,col.names2,col.names3)
out.HO = matrix(NA, nrow=S*length(K.grid), ncol=length(col.names), dimnames=list(NULL,col.names))
out.HE = matrix(NA, nrow=S*length(K.grid), ncol=length(col.names), dimnames=list(NULL,col.names))
col12 = length(c(col.names1,col.names2));

message("Simulations began for Model ", m,". Time: ", Sys.time()); showevery=.50
row=1;

# s=1; K.i=3;
for (s in 1:S) {
    if (max(s==seq(0,S,S*showevery))==1) {message("Simulations Completed: ",s," of ",S," (",round(s/S*100),"%) - ", Sys.time())}

    for (K.i in K.grid) {
        if (1  <= m & m <= 6 ) {dgp = dgp.fn(m=m,n=n,K.i=K.i,K=K.grid.LM[K.i]); W = dgp$W}
        if (7  <= m & m <= 9 ) {N = K.grid.FE.N[K.i]; T = K.grid.FE.T[K.i];
                                dgp = dgp.fn(m=m,n=N*T,K.i=K.i,K=50);           W = fixed.effects(N,T,G=0)}
        if (10 <= m & m <= 12) {N = K.grid.FE.N[K.i]; T = K.grid.FE.T[K.i];
                                dgp = dgp.fn(m=m,n=N*T,K.i=K.i,K=50);           W = fixed.effects(N,T,G=floor(N/3))}
        if (13 <= m & m <= 15) {dgp = dgp.fn(m=m,n=n,K.i=K.i,K=6);              W = gen.P(dgp$W[,-1],K.grid.PL[K.i])}

        K = ncol(W); qr.W = qr(W); if (qr.W$rank==K){
            M = -tcrossprod(qr.Q(qr.W)); diag(M) = 1 + diag(M);
            kappa = try(chol2inv(chol(M*M)), TRUE);
            out.HO[row,1:col12] = c(n,d,beta,m,0,s,K,qr.W$rank,lsfit(dgp$Y.HO,dgp$X.HO,M,kappa,beta,K));
            out.HE[row,1:col12] = c(n,d,beta,m,1,s,K,qr.W$rank,lsfit(dgp$Y.HE,dgp$X.HE,M,kappa,beta,K));
        }

        ## BOOTSTRAP INFERENCE -- BEGIN
        if(B>0 & s<floor(S/1)) {
        col.names.B = c("beta.hat.B","T.HO0.B","T.HO1.B","T.HC0.B","T.HC1.B","T.HC2.B","T.HC3.B","T.HC4.B","T.HCK.B")
        out.HO.B = matrix(NA, nrow=B, ncol=length(col.names.B), dimnames=list(NULL,col.names.B))
        out.HE.B = matrix(NA, nrow=B, ncol=length(col.names.B), dimnames=list(NULL,col.names.B))
        for (b in 1:B){
            id.B = sample.int(nrow(W), size = nrow(W), replace = TRUE)
            if (1  <= m & m <= 12) {W.B = W[id.B,]}
            if (13 <= m & m <= 15) {W.B = gen.P(dgp$W[id.B,-1],K.grid.PL[K.i])}
            qr.W.B = qr(W.B); if (qr.W.B$rank==K){
            M = -tcrossprod(qr.Q(qr.W.B)); diag(M) = 1 + diag(M);
            kappa = try(chol2inv(chol(M*M)), TRUE);
            out.HO.B[b,] = lsfit(dgp$Y.HO[id.B],dgp$X.HO[id.B],M,kappa,beta=out.HO[row,"beta.hat"],K)[c(1,10:17)]
            out.HE.B[b,] = lsfit(dgp$Y.HE[id.B],dgp$X.HE[id.B],M,kappa,beta=out.HE[row,"beta.hat"],K)[c(1,10:17)]
            }
        }

        out.HO[row,(col12+1):ncol(out.HO)] = c(var(out.HO.B[,1]), c(t(apply(out.HO.B[,2:9],2,function(x) quantile(x, probs=c(.025,.975), na.rm=TRUE, type=1)))))
        out.HE[row,(col12+1):ncol(out.HE)] = c(var(out.HE.B[,1]), c(t(apply(out.HE.B[,2:9],2,function(x) quantile(x, probs=c(.025,.975), na.rm=TRUE, type=1)))))
        }
        ## BOOTSTRAP INFERENCE -- END

        row=row+1;
    }
}

## Save final table
filename = paste0("output/parts/output_m",m,"_hom_",as.integer(args[[1]]),".csv");
if (file.exists(filename)) write.table(out.HO, file=filename, append=TRUE, col.names=FALSE, row.names=FALSE, sep=",") else write.table(out.HO, file=filename, row.names=FALSE, sep=",")
filename = paste0("output/parts/output_m",m,"_het_",as.integer(args[[1]]),".csv");
if (file.exists(filename)) write.table(out.HE, file=filename, append=TRUE, col.names=FALSE, row.names=FALSE, sep=",") else write.table(out.HE, file=filename, row.names=FALSE, sep=",")
}

merge.output(cpus=200,models=models)

## TO TEST:
## rm(list=ls(all=TRUE)); source("main.R"); source("tables.R"); rownames(out.EC)=rowname; colnames(out.EC)=c(c("HO0","HO1","HC0","HC1","HC2","HC3","HC4","HCK"),rep(".",8)); print(out.EC)


