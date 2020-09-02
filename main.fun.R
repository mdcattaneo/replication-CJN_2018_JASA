#########################################################################################
## R CODE FOR CATTANEO-JANSSON-NEWEY (2016)
## Inference in Linear Regression Models with Many Covariates and Heteroskedasticity
## DATE: 20-Dec-2016
## FUNCTIONS
#########################################################################################
path = "/mainfiles/cattaneo/Dropbox/research/HCSEManyCov/simulations/"
#path = "Z:/research/HCSEManyCov/simulations/"


#########################################################################################
## ERRORS DISTRIBUTIONS
#########################################################################################
rmixnorm = function(n,m1,m2,s1,s2,alpha) {I = runif(n)<alpha; rnorm(n,mean=ifelse(I,m1,m2),sd=ifelse(I,s1,s2));}
dmixnorm = function(x,m1,m2,s1,s2,alpha) return(alpha * dnorm(x,m1,s1) + (1-alpha)*dnorm(x,m2,s2))

if(FALSE) {
    postscript(file=paste0(path,"Figure1.ps"), horizontal=F)
    par(mfrow=c(3,1))
    curve(dnorm, from=-3, to=3, main="Normal Distribution")
    fun = function(x) dmixnorm(x,m1=-1/2,m2=1/2,s1=1/2,s2=1,alpha=1/2); curve(fun, from=-3, to=3, main="Asymmetric Distribution")
    fun = function(x) dmixnorm(x,m1=-3/2,m2=3/2,s1=1/2,s2=1/2,alpha=1/2); curve(fun, from=-3, to=3, main="Bimodal Distribution")
    dev.off()
}

#########################################################################################
## DATA GENERATING PROCESSES
#########################################################################################
## TO TEST:
## popval = read.csv("main.popval.csv", row.names=1); source("main.fun.R"); tmp = dgp.fn(m=7,n=1000,K.i=4,K=6); cbind(tmp$Y.HE,tmp$X.HE,tmp$v.HE,tmp$u.HE)[1:10,]
## var(cbind(tmp$Y.HE,tmp$X.HE,tmp$v.HE,tmp$u.HE))
dgp.fn = function(m,n,K.i,K,varkappa.v=NULL,varkappa.u=NULL) {
    if(is.null(varkappa.v)) varkappa.v = popval[is.na(popval[,"m"])==FALSE & popval[,"m"]==m & popval[,"K.i"]==K.i & popval[,"vartheta"]==1, "varkappa.v"];
    if(is.null(varkappa.u)) varkappa.u = popval[is.na(popval[,"m"])==FALSE & popval[,"m"]==m & popval[,"K.i"]==K.i & popval[,"vartheta"]==1, "varkappa.u"];
    trunc = function(x){c=1.5; return(-c*(x < -c) + x*(abs(x) <= c) + c*(x > c))};
    sigma.v.fn = function(Wg,vartheta=.75) return((1+Wg^2)^vartheta);
    sigma.u.fn = function(Xb,Wg,vartheta=.75) return((1+(trunc(Xb)+Wg)^2)^vartheta);

    #########################################################################################
    ## MODELS 1--12: Linear Regression Models
    #########################################################################################
    if(1 <= m & m <= 12){
        if(m==1 | m==4 | m==7 | m==10) {v = matrix(rnorm(n),n,1);                             u = matrix(rnorm(n),n,1);}
        if(m==2 | m==5 | m==8 | m==11) {v = rmixnorm(n,m1=-1/2,m2=1/2,s1=1/2,s2=1,alpha=1/2); u = rmixnorm(n,m1=-1/2,m2=1/2,s1=1/2,s2=1,alpha=1/2);}
        if(m==3 | m==6 | m==9 | m==12) {v = rmixnorm(n,m1=-3/2,m2=3/2,s1=1/2,s2=1,alpha=1/2); u = rmixnorm(n,m1=-3/2,m2=3/2,s1=1/2,s2=1,alpha=1/2);}

        W = matrix(1,n,1+K); if (K>0) {
            if(1  <= m & m <= 3 ) W[,2:(K+1)] = abs(matrix(rnorm(n*K,mean=0,sd=1),n)) >= qnorm(.99);
            if(4  <= m & m <= 6 ) W[,2:(K+1)] = matrix(runif(n*K,min=-1,max=1),n);
            if(7  <= m & m <= 9 ) W[,2:(K+1)] = matrix(runif(n*K,min=-1,max=1),n);
            if(10 <= m & m <= 12) W[,2:(K+1)] = matrix(runif(n*K,min=-1,max=1),n);
        }
        Wg = W%*%rep.int(1,1+K);

        ## DGP: HOMOSKEDASTIC ERRORS
        X.HO = 0        + v;
        Y.HO = X.HO + 0 + u;

        ## DGP: HETEROSKEDASTIC ERRORS
        v.HE = varkappa.v * sigma.v.fn(Wg) * v     ; X.HE = 0        + v.HE;
        u.HE = varkappa.u * sigma.u.fn(X.HE,Wg) * u; Y.HE = X.HE + 0 + u.HE;
    }
    #########################################################################################
    ## MODELS 13-14-15: Partially Linear
    #########################################################################################
    g = function(Z) {apply(Z,1,function(x){tmp=exp(-crossprod(x)^(1/2))});} #;return(tmp/(1+tmp))
    h = function(Z) {apply(Z,1,function(x){tmp=exp(crossprod(x)^(1/2))});} #;return(tmp/(1+tmp))

    if(m==13 | m==14 | m==15){
        if(m==13) {v = matrix(rnorm(n),n,1);                             u = matrix(rnorm(n),n,1);}
        if(m==14) {v = rmixnorm(n,m1=-1/2,m2=1/2,s1=1/2,s2=1,alpha=1/2); u = rmixnorm(n,m1=-1/2,m2=1/2,s1=1/2,s2=1,alpha=1/2);}
        if(m==15) {v = rmixnorm(n,m1=-3/2,m2=3/2,s1=1/2,s2=1,alpha=1/2); u = rmixnorm(n,m1=-3/2,m2=3/2,s1=1/2,s2=1,alpha=1/2);}

        W = matrix(1,n,1+K); if (K>0) W[,2:(K+1)] = matrix(runif(n*K,min=-1,max=1),n);
        Wg = W%*%rep.int(1,1+K);

        ## DGP: HOMOSKEDASTIC ERRORS
        X.HO =        h(W) + v;
        Y.HO = X.HO + g(W) + u;

        ## DGP: HETEROSKEDASTIC ERRORS
        v.HE = varkappa.v * sigma.v.fn(Wg) * v     ; X.HE =        h(W) + v.HE;
        u.HE = varkappa.u * sigma.u.fn(X.HE,Wg) * u; Y.HE = X.HE + g(W) + u.HE;
    }

    return(list(Y.HO=Y.HO, X.HO=X.HO,
                Y.HE=Y.HE, X.HE=X.HE, W=W,
                v.HE=v.HE, u.HE=u.HE))
}

#########################################################################################
## DGPs: Population Scaling
#########################################################################################
## setwd("/mainfiles/cattaneo/Dropbox/research/HCSEManyCov/simulations"); source("main.fun.R"); gen.table(n=700,models=1:15,K.grid=1:5)

gen.table = function(n,models,K.grid){
    ## NOTE: this function uses environment variables n (sample size) and K.grid (grid for K)
    vartheta=1
    dimnames = list(NULL,c("m","K.i","vartheta","varkappa.v","varkappa.u"))
    popval = matrix(NA,nrow=length(models)*length(K.grid),ncol=length(dimnames[[2]]),dimnames=dimnames)

    I = 500000; row=1;
    for (m in models){
        message("\nComputing Scaling Constants for Model ",m,".")
        for (K.i in K.grid){
            K.grid.LM = floor(n*c(0,.1,.2,.3,.4));
            K.grid.PL = c(0,1.5,2.5,3.5,4.5);
            K.grid.FE.T = c(n,10,5,4,3)
            K.grid.FE.N = floor(n/K.grid.FE.T)
            if (1  <= m & m <= 6)  K = K.grid.LM[K.i];
            if (7  <= m & m <= 9)  K = 6;
            if (10 <= m & m <= 15) K = 50;

            dgp = dgp.fn(m=m,n=I,K.i=K.i,K=K,varkappa.v=1,varkappa.u=1)
            sigma2.v = as.numeric(crossprod(dgp$v.HE)/I);

            dgp = dgp.fn(m=m,n=I,K.i=K.i,K=K,varkappa.v=1/sqrt(sigma2.v),varkappa.u=1)
            sigma2.u = as.numeric(crossprod(dgp$u.HE)/I);
            popval[row,] = c(m,K.i,vartheta,1/sqrt(sigma2.v),1/sqrt(sigma2.u))

            row=row+1
        }
        write.csv(popval, file=paste0(path,"main.popval.csv"))
    }
}

#########################################################################################
## QR-based (X'X)^(-1)
#########################################################################################
#XXinv = function(x) tcrossprod(solve(qr.R(qr(x))))
XXinv = function(x) chol2inv(chol(crossprod(x)))

################################################################################
## BASIS OF APPROXIMATION
################################################################################
gen.P = function(Z,K) {
    dim.basis = factorial(floor(K)+ncol(Z))/(factorial(floor(K))*factorial(ncol(Z)))
    dim.basis.adj = dim.basis + ncol(Z)*(K-floor(K) == 0.5)

    out = matrix(1,nrow(Z),dim.basis.adj)
    if (dim.basis>1) out[,2:dim.basis] = poly(Z,degree=floor(K),raw=TRUE);
    if (dim.basis.adj>dim.basis) for (j in 1:ncol(Z)) out[,dim.basis+j] = Z[,j]^ceiling(K);

    return(out)
}
## Z=matrix(runif(3*6),3); for (K in seq(0,6,by=.5)) {tmp=gen.P(Z,K); message(K, " | ",ncol(tmp), " | ", round(ncol(tmp)/700,3))}

################################################################################
## FIXED EFFECTS
################################################################################
fixed.effects = function(N,T,G) {
    ## Mapping: n = N*T | K = (N+G)/n | G <= N-1 | T >= 3
    GS = floor(N/G)
    out = matrix(0,N*T,N+G)
    ## Unit fixed effects
    out[,1] = rep.int(1,N*T)
    if (N>1) for (j in 2:N) out[((j-1)*T+1):((j-1)*T+T),j] = 1

    ## Group fixed effects
    if (G>0) for (j in 1:G) out[seq.int(1+(j-1)*GS*T,(j<G)*j*GS*T+(j==G)*N*T,by=T),N+j] = 1
    return(out)
}
## source("main.fun.R"); W = fixed.effects(N=20,T=3,G=50); W; solve(crossprod(W))


#########################################################################################
## LSFIT
#########################################################################################
lsfit = function(Y,X,M,kappa,beta,K) {
    MX = M%*%X; XMX_inv = XXinv(MX); XMY = crossprod(MX,Y); n = length(Y)

    beta.hat = XMX_inv%*%XMY

    u2.hat = (M%*%Y - MX%*%beta.hat)^2

    VHO0 = sum(u2.hat)/n * XMX_inv
    THO0 = (beta.hat - beta) / sqrt(VHO0)

    VHO1 = VHO0 * n/(n-1-K)
    THO1 = (beta.hat - beta) / sqrt(VHO1)

    VHC0 = XMX_inv * crossprod(MX*u2.hat,MX) * XMX_inv
    THC0 = (beta.hat - beta) / sqrt(VHC0)

    VHC1 = VHC0 * n/(n-1-K)
    THC1 = (beta.hat - beta) / sqrt(VHC1)

    M.diag = diag(M); M.exp = apply(matrix(M.diag,ncol=1),1,function(x) min(c(4,n*x/K)))
    VHC2 = XMX_inv * crossprod(MX*(u2.hat/M.diag),MX) * XMX_inv
    THC2 = (beta.hat - beta) / sqrt(VHC2)

    VHC3 = XMX_inv * crossprod(MX*(u2.hat/(M.diag^2)),MX) * XMX_inv
    THC3 = (beta.hat - beta) / sqrt(VHC3)

    VHC4 = XMX_inv * crossprod(MX*(u2.hat/(M.diag^(M.exp))),MX) * XMX_inv
    THC4 = (beta.hat - beta) / sqrt(VHC4)

    if (is.matrix(kappa)){
        VHK = XMX_inv * crossprod(MX*(kappa%*%u2.hat),MX) * XMX_inv
        THK = (beta.hat - beta) / sqrt(VHK)
    } else {VHK = NA; THK=NA;}

    return(c(beta.hat,VHO0,VHO1,VHC0,VHC1,VHC2,VHC3,VHC4,VHK,
             THO0,THO1,THC0,THC1,THC2,THC3,THC4,THK))
}

#########################################################################################
## MERGE OUTPUT FILES
#########################################################################################
## source('main.fun.R'); merge.output()
merge.output = function(cpus=200,models=1:15){
    output_hom = list(NA)
    output_het = list(NA)
    
    for (m in models) {
        files_hom = NULL;
        files_het = NULL
        for (s in 1:cpus) {
            filename_hom = paste0("output/parts/output_m",m,"_hom_",s,".csv")
            if (file.exists(filename_hom)) files_hom = c(files_hom, filename_hom)
            filename_het = paste0("output/parts/output_m",m,"_het_",s,".csv")
            if (file.exists(filename_het)) files_het = c(files_het, filename_het)
        }
        if (length(filename_hom)>0 & length(filename_het)>0){
            write.table(do.call("rbind", lapply(files_hom, function(x) read.csv(x, row.names=NULL))), file = paste0("output/output_m",m,"_hom.csv"), row.names=FALSE, sep=",")
            write.table(do.call("rbind", lapply(files_het, function(x) read.csv(x, row.names=NULL))), file = paste0("output/output_m",m,"_het.csv"), row.names=FALSE, sep=",")
            message("Model ",m," done -- ",length(files_hom)," | ",length(files_het)," Files read.")
        } else message("Model ",m," not available.")
    }
}

