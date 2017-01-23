# Experimental code, subject to random changes
# 
# Author: mg14
###############################################################################

CoxMFX <- function(Z, X=NULL, surv, groups = rep(1, ncol(Z)), which.mu = unique(groups), tol=1e-3, max.iter=50, sigma0 = 0.1, nu = 0,  penalize.mu = FALSE, sigma.hat=c("df","MLE","REML","BLUP"), verbose=FALSE){
	if(class(Z)=="data.frame"){
		Z = as.matrix(Z)
		Z.df <- TRUE
	}else
		Z.df <- FALSE
	if(is.null(colnames(Z)))
		colnames(Z) <- make.names(1:ncol(Z))
	sigma.hat = match.arg(sigma.hat)
	o <- order(groups)
	Z <- Z[,o]
	groups <- factor(groups[o])
	uniqueGroups <- levels(groups)
	ZZ <- lapply(uniqueGroups, function(i) Z[,groups==i, drop=FALSE])
	names(ZZ) <- uniqueGroups
	sumZ <- sapply(which.mu, function(i) rowSums(ZZ[[i]]))
	nGroups = length(uniqueGroups)
	sigma2 <- sigma0ld <- rep(ifelse(sigma0>0, sigma0,1), nGroups)
	iter = 1
	mu <- mu0ld <- rep(0, nGroups)
	names(mu) <- uniqueGroups
	if(!is.null(X)){
		f <- ncol(X)
		if(is.null(colnames(X))) colnames(X) <- paste0("X.", 1:ncol(X))
	} 
	else f <- 0
	beta = rep(1,ncol(Z)+length(which.mu) + f)
	beta0ld = rep(0,ncol(Z)+length(which.mu) + f)
	sigma2.mu = 42
	if(!is.null(which.mu)) 
		if(!penalize.mu)
			fixedEff <- "sumZ" 
		else
			fixedEff <- "ridge(sumZ, theta=1/sigma2.mu, scale=FALSE)"
	else fixedEff <- character(0)
	if(!is.null(X)){
		fixedEff <- paste(fixedEff, "+ X")
	}
	while((max(abs(beta-beta0ld)) > tol | max(abs(mu - mu0ld)) > tol | max(abs(sigma2 - sigma0ld)) > tol) & iter < max.iter){
		beta0ld = beta
		sigma0ld <- sigma2
		mu0ld <- mu
		formula <- formula(paste("surv ~", paste(c(sapply(1:nGroups, function(i) paste("ridge(ZZ[[",i,"]], theta=1/sigma2[",i,"], scale=FALSE)", sep="")), 
										#ifelse(!is.null(which.mu),"ridge(sumZ, theta=1/sigma.mu, scale=FALSE)","")), 
										fixedEff), 
								collapse=" + ")))
		fit <- coxph(formula)
		if(any(is.na(coef(fit)))){
			warning(paste("NA during estimation (iter: ", iter, ", coef: ", paste(which(is.na(coef(fit)[order(o)])), sep=","), ")", sep=""))
			break
		}
		if(!is.null(which.mu))
			mu[which.mu] <- coef(fit)[ncol(Z) + seq_along(which.mu)]
		if(verbose) cat("mu", mu, "\n", sep="\t")
		names(fit$df) <- c(uniqueGroups, rep("Offset", length(which.mu)>0))
		if(verbose) cat("df", fit$df,"\n", sep="\t")
		sigma2 = sapply(uniqueGroups, function(i){
					index <- which(groups==i) #& fit$coefficients > beta.thresh
					if(sigma.hat=="BLUP")
						(nu * sigma0 + sum((fit$coefficients[index])^2 ))/(nu + length(index))
					else if(sigma.hat=="df")
						(nu * sigma0 + sum((fit$coefficients[index])^2 ))/(nu + fit$df[i])
					else if(sigma.hat == "MLE")
						(nu * sigma0 + sum((fit$coefficients[index])^2 ) + sum(diag(solve(solve(fit$var)[index,index]))))/(nu + length(index))
					else if(sigma.hat == "REML")
						(nu * sigma0 + sum((fit$coefficients[index])^2 ) + sum(diag(fit$var)[index]))/(nu + length(index))
				})
		if(verbose) {
			cat("sigma2", sigma2, "\n", sep="\t")
			cat("loglik:", fit$loglik - c(0,fit$penalty[2] + 1/2 * sum(log(sigma2[groups]))),"\n", sep="\t")
		}
		if(penalize.mu){
			if(sigma.hat=="BLUP")
				sigma2.mu = (sigma0 * nu + sum((mu-0)^2)) / (nu + length(mu))
			else if(sigma.hat=="df")
				sigma2.mu = (sigma0 * nu + sum((mu-0)^2)) / (nu + fit$df["Offset"])
			else if(sigma.hat == "MLE")
				sigma2.mu = (nu * sigma0 + sum((mu-0)^2 ) + sum(diag(solve(solve(fit$var)[-(1:ncol(Z)),-(1:ncol(Z))]))))/(nu + length(mu))
			else if(sigma.hat == "REML")
				sigma2.mu = (nu * sigma0 + sum((mu-0)^2 ) + sum(diag(fit$var)[-(1:ncol(Z))]))/(nu + length(mu))
		}
		
		beta = fit$coefficients
		iter = iter+1
	}
	if(iter == max.iter)
		warning("Did not converge after ", max.iter, " iterations.")
	fit$iter[1] <- iter
	fit$sigma2 = sigma0ld
	names(fit$sigma2) <- uniqueGroups
	fit$sigma2.mu = sigma2.mu
	fit$mu = mu
	fit$Z = Z[,order(o)]
	fit$surv = surv
	C <- rbind(diag(1, ncol(Z)),t(as.matrix(MakeInteger(groups)[which.mu]))) ## map from centred to uncentred coefficients 
	if(!is.null(X))
		C <- cbind( rbind(C, matrix(0, ncol=ncol(C), nrow=f) ), rbind( matrix(0, ncol=f, nrow=nrow(C)), diag(1, f)))
	fit$groups = groups[order(o)]
	var = fit$var
	var2 = fit$var2
	colnames(var) <- rownames(var) <- colnames(var2) <- rownames(var2) <- rownames(C) <- c(colnames(Z), which.mu, colnames(X))
	colnames(C) <- c(colnames(Z), colnames(X))
	p <- ncol(Z)
	i <- 1:nrow(C)
	i[1:p] <- order(o)
	j <- 1:ncol(C)
	j[1:p] <- order(o)
	fit$C <- C[i,j]
	fit$Hinv <- var[i,i] ## Hinv 
	fit$V <- var2[i,i] ## Hinv I Hinv
	fit$z <- (fit$coefficients / sqrt(diag(var)))[i] ## z-scores of centred coefficients
	fit$z2 <- (fit$coefficients / sqrt(diag(var2)))[i] ## z-scores of centred coefficients (var2)
	fit$var = (t(C) %*% var %*% C)[j,j] ## covariance of uncentred coef
	fit$var2 = (t(C) %*% var2 %*% C)[j,j] ## covariance of uncentred coef (var2)
	w <- ncol(Z) + seq_along(which.mu)
	fit$mu.var = var[w,w] ## covariance of mean
	fit$mu.var2 = var2[w,w] ## covariance of mean (var2)
	fit$means = fit$means[1:p][j]
	fit$coefficients <- (fit$coefficients %*% C)[j]
	names(fit$means) <- names(fit$coefficients) <-  colnames(Z)[j]
	fit$terms <- fit$terms[1:length(uniqueGroups)]
	fit$penalized.loglik <- fit$loglik[2] - fit$penalty[2] - 1/2 * sum(log(fit$sigma2[groups]))
	## Fake call for predict.coxph and survfit.coxph
	call <- match.call()
	if(Z.df){
		call["data"] <- call["Z"]
		formula <- as.formula(paste(as.character(call["surv"]),"~",paste(colnames(Z)[j], collapse="+")))
	}else{
		formula <- as.formula(paste(as.character(call["surv"]),"~",as.character(call["Z"])))
	}
	attr(formula,".Environment") <- parent.frame()
	fit$formula <- formula
	call["formula"] <- call("foo",formula=formula)["formula"]
	fit$terms <- terms(formula)
	fit$call <- call
	class(fit) <- c("CoxMFX", class(fit))
	return(fit)
}

SimSurvTDNonp <- function(dataFrame, coef, time0, time1, event, timeTpl, coefTpl) {
	w <- which(timeTpl < time1 & timeTpl > time0)
	index <- c(1:nrow(dataFrame), w) 
	d <- dataFrame[index,]
	d$index <- index
	d$time0 <- c(time0, timeTpl[w])
	d$time1 <- c(pmin(time1, timeTpl, na.rm=TRUE), time1[w])
	d$transplant <- c(rep(0,nrow(dataFrame)), rep(1, length(w)))
	e <- c(event, event[w])
	e[w] <- 0
	d$event <- e
	
	i <- 1
	
	simTime1 <- simTplTimes <- simDeathTimes <- simCensTimes <- rep(NA,nrow(dataFrame))
	#H0 <- basehaz(coxph(Surv(time0, time1, event ) ~ as.matrix(dataFrame) %*% coef), centered = FALSE)
	
	for(idx in list(is.na(timeTpl), !is.na(timeTpl))){
		surv <- Surv(time0[idx], time1[idx], event[idx] )
		
		risk0 <- as.matrix(dataFrame[idx, ]) %*% coef + ifelse(i>1, coefTpl,0) 
		H0 = basehaz(coxph(surv ~ risk0), centered = FALSE)
		
		## Simulate deaths times
		FHazInv <- splinefun(c(0,H0$hazard), c(0,H0$time), method="monoH.FC")
		FHaz <- splinefun(c(0,H0$time), c(0,H0$hazard),  method="monoH.FC")
		n = length(risk0)
		
		## Simulate censoring times
		nCens <- sum(event==0, na.rm=TRUE)
		FSurv <- splinefun(exp(-H0$hazard) ~ H0$time, method="monoH.FC") ## Cum survival dist
		x <- sort(na.omit(time1[event==0])) ## observed (conditioned) censored times
		t <- table(x) / nCens ## empirical Dist of cond. censoring times
		y <- t / FSurv(unique(x)) ## Need to get unconditional distribution
		FCens <- cumsum(y)
		FCens <- FCens/max(FCens) ## Unconditioned cens.
		FCensInvDist <- splinefun(unique(x) ~ FCens, method="monoH")
		simCensTimes[idx] <- FCensInvDist(runif(n,0,1)) ## Simulate censoring times
		
		simDeathTimes[idx] = FHazInv(rexp(n, exp(risk0)))
		#w <- simDeathTimes > simTplTimes & simDeathTimes < simCensTimes
		#w <- sample(seq_along(simDeathTimes), length(tplTimes))
		#simDeathTimes[w] <- FHazInv(FHaz(simDeathTimes[w]) + rexp(length(w), exp(coefTpl + risk0[w])))
		simTime1[idx] <- pmin(simDeathTimes[idx], simCensTimes[idx])

		## TD observations
		if(i > 1){
			w <- which(!is.na(timeTpl))
			tplTimes <- na.omit(timeTpl)
			x <- sort(tplTimes) ## observed (conditioned) censored times
			t <- table(x) / length(tplTimes) ## empirical Dist of cond. censoring times
			#y <- t / FSurv(unique(x)) ## Need to get unconditional distribution
			FTpl <- cumsum(t)
			FTpl <- FTpl/max(FTpl) ## Unconditioned cens.
			FTplDist <- splinefun(FTpl ~ unique(x), method="monoH")
			FTplInvDist <- splinefun(unique(x) ~ FTpl, method="monoH")
			#simTplTimes <- FTplInvDist(runif(n,0,1)) ## Simulate censoring times
			
			#simTplTimes[w] <- pmin(FTplInvDist(runif(length(w),0,pmax(0,pmin(1,FTplDist(simTime1[w]))))), simTime1[w]-0.5)
			simTplTimes[w] <- runif(length(w), min(simTime1[w]), simTime1[w])
		}
		i <- i+1
	}
	## Put together
	index <- c(1:nrow(dataFrame), w) 
	d <- dataFrame[index,, drop=FALSE]
	d$index <- index
	d$time0 <- c(time0, simTplTimes[w])
	d$time1 <- c(pmin(simTime1, simTplTimes, na.rm=TRUE), simTime1[w])
	d$transplant <- c(rep(0,nrow(dataFrame)), rep(1, length(w)))
	evt <- simDeathTimes < simCensTimes
	e <- c(evt, evt[w])
	e[w] <- 0
	d$event <- e
	return(d)
}

DensityEstimates <- function(coxRFX, newx = range(coef(coxRFX)), n = 100){
	x <- seq(newx[1], newx[2], length.out=n)
	c <- coef(coxRFX)
	v <- diag(coxRFX$var)
	z <- mapply(function(i,j) dnorm(x, i, j), c, sqrt(v))
	sapply(levels(coxRFX$groups), function(g) rowMeans(z[,coxRFX$groups==g, drop=FALSE]))
}


#### Paralllelized stepwise forward selection
parStep = function(X, surv, criterion = "AIC", max.iter = ncol(X), mc.cores=1, verbose=FALSE){
	select = numeric()
	loglik = numeric()
	penalty = ifelse(criterion == "AIC",1,log(nrow(X))/2)
	while(length(select) < max.iter){
		k = length(select) + 1
		scope = setdiff(1:ncol(X), select)
		if(is.null(scope))
			break
		l = mclapply(scope, function(i) {
					x = as.data.frame(X[,c(select,i)])
					t = try(coxph(surv~., data=x)$loglik[2])
					ifelse(class(t)!="try-error",t,-Inf)
				}, mc.cores=mc.cores)
		logliks = unlist(l)		
		add = scope[which.max(logliks)]
		if(k>1)
			if(all(logliks <= loglik[k-1] + penalty))
				break
		loglik = c(loglik,max(logliks))
		select = c(select,add)
		if(verbose) cat(".")
	}
	if(verbose) cat("\n")
	return(data.frame(select=select,loglik=loglik,AIC = - 2*loglik + 2 * 1:length(loglik), BIC = - 2 * loglik +  1:length(loglik) * log(nrow(X)) ))
}

GetPairs <- function(names, scope){
	pairs <- strsplit(scope, "\\.")
	lapply(pairs, function(s){
				a <- grep(paste(s[1],"$", sep=""), names, value=TRUE)
				b <- grep(paste(s[2],"$", sep=""), names, value=TRUE)
				c(a,b)
			})
}


show.CoxRFX <- function(x){
	summary(x)
	cat("\nCoefficients:\n")
	WaldTest(x)
}

ConcordanceFromVariance <- function(x) {
	.cfv <- function(x){
		stopifnot(x >= 0)
		if(x == 0)
			return(.5)			
		f <- function(x, sigma2) dnorm(x, sd=sqrt(2*sigma2)) / (1 +  exp(- x))
		2*integrate(f, 0, 10*sqrt(x), sigma2=x)$value
	}
	sapply(x, .cfv)
}

SimCoef <- function(coxRFX=NULL, groups = coxRFX$groups, mu=coxRFX$mu, sigma2=coxRFX$sigma2){
	c <- rnorm(length(groups), mu[groups], sqrt(sigma2[groups]))
	names(c) <- names(groups)
	c
}