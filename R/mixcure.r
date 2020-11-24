mixcure = function(lformula, iformula, data, lmodel, imodel, postuncure = NULL,
				   emmax = 100, eps = 1e-4, savedata = FALSE, debug = FALSE)
{
	call <- match.call()

	if (missing(lformula)) stop("Latency formula not provided")
	if (missing(lmodel)) lmodel = list(fun = "coxph")
	lmodel$formula = lformula
	lfun = lmodel$fun
	lmodel$fun = NULL
	if (missing(iformula)) stop("Incidence formula not provided")
	if (missing(imodel)) imodel = list(fun = "glm", family = binomial())
	imodel$formula = iformula
	ifun = imodel$fun
	imodel$fun = NULL


	survtime = model.response(model.frame(formula = lformula,
									data = data, na.action = na.pass))
	ttime = survtime[, 1]
	cens = survtime[, 2]

	if (is.null(postuncure)) {
		postuncure = rep(1, length(ttime))
		postuncure[cens==0] = seq(1.0, 0.0, along = postuncure[cens==0])
	}

	if(debug) fcat("EM:")
	emdetail = NULL
	for (em in 1:emmax) {
		if(debug) fcat(".")

		imodel$data = cbind(data, postuncure = postuncure)
		ifit = do.call(paste("incidence.", ifun, sep = ""), imodel)

		lmodel$data = cbind(data, postuncure = postuncure)
		lfit = do.call(paste("latency.", lfun, sep = ""), lmodel)

		survprob = lfit$survprob
        postuncure = cens + (1-cens)*ifit$uncureprob *
			survprob/(1-ifit$uncureprob + ifit$uncureprob * survprob)


        newres = sum(abs(c(loglik(lfit), loglik(ifit))))
		if (em==1) oldres = 0.0
		prec = abs(1.0 - oldres/newres)

		emdetail = rbind(emdetail, c(em, prec))


		if (is.na(prec) || is.nan(prec) || prec <= eps) {
			oldres = newres
			break
		}
		oldres = newres
	}
	if(debug) fcat("\n")
	lmodel$data = NULL
	lmodel$fun = lfun
	imodel$data = NULL
	imodel$fun = ifun
	val = list(ifit = ifit, lfit = lfit, survprob = survprob,
		postuncure = postuncure, em = emdetail, imodel = imodel, lmodel = lmodel,
		emmax = emmax, eps = eps)
	val$call = call
	if (savedata) val$data = data
	class(val) = "mixcure"
	val
}


coef.mixcure = function(object, ...)
{
	list(latency = coef(object$lfit), incidence = coef(object$ifit))
}

print.mixcure = function(x, ...)
{
	if (!is.null(cl <- x$call)) {
		cat("Call:\n")
		dput(cl)
	}
	cat("\nLatency model:\n")
	print(x$lfit$fit, ...)
	cat("\nIncidence model:\n")
	print(x$ifit$fit, ...)
	cat("\nEM algorithm finished ", x$em[nrow(x$em), 1],
		"iterations with final error", x$em[nrow(x$em), 2], "\n")
	cat("\nNote: the standard errors and p-values reported may not be correct.\n")
	invisible(x)
}

summary.mixcure = function(object, R = 100, index = 1:2, ...)
{
	if (class(object)!="mixcure") stop("mixcure object is missing")

	if (is.null(object$stderr)) {
		if (is.null(object$data))
			stop("mixcure object should be fit with savedata = TRUE")

		bootfun = function(y, mixcurefit)
			unlist(coef(mixcure(
				lformula = mixcurefit$lmodel$formula,
				lmodel = mixcurefit$lmodel,
				iformula = mixcurefit$imodel$formula,
				imodel = mixcurefit$imodel,
				data = y, emmax = mixcurefit$emmax, eps = mixcurefit$eps)))

		object$varboot = censboot(object$data, bootfun, R = R, index = index,
			mixcurefit = object)
		object$stderr = sqrt(apply(object$varboot$t, 2, var))
		object$R = R
	}

	if (!is.null(cl <- object$call)) {
		cat("Call:\n")
		dput(cl)
	}

	coefficient = c(coef(object$lfit), coef(object$ifit))
	info = cbind(coefficient = coefficient,
		stderr = object$stderr,
		zscore = coefficient/object$stderr,
		pvalue = 2*(1-pnorm(abs(coefficient/object$stderr))))
	dimnames(info)[[1]] = names(coefficient)

	cat("\nLatency model for uncured:\n")
	if (length(coef(object$lfit)>0))
		print(info[1:length(coef(object$lfit)), , drop = FALSE], ...)

	cat("\nIncidence model:\n")
	print(info[(1+length(coef(object$lfit))):nrow(info), , drop = FALSE], ...)

	invisible(object)
}


predict.mixcure = function(object, newdata, times, ...)
{
	call <- match.call()
	if (is.null(newdata) || is.null(times))
		stop("newdata or times is not provided for prediction.")

	prediction = list()
	uncureprob = curepred(object$ifit, newdata = newdata)
	prediction$cure = cbind(uncureprob = uncureprob, cureprob = 1- uncureprob)

	prediction$uncuresurv = survpred(object$lfit, newdata = newdata,
					times = times)

	prediction$surv = lapply(mapply("c", uncureprob, prediction$uncuresurv,
							SIMPLIFY = FALSE), function(x)x[-1]*x[1] + 1 - x[1])

	prediction$times = times
	prediction$call = call
	class(prediction) = "predict.mixcure"
	prediction
}

print.predict.mixcure = function(x, digits = 4,...)
{
	if (!is.null(cl <- x$call)) {
		cat("Call:\n")
		dput(cl)
	}
	cat("\nEstimated cure and uncured probability:\n")
	print(x$cure, digits = digits, ...)
	cat("\nEstimated survival probability for uncured:\n")

	print(lapply(x$uncuresurv, function(y)cbind(time = x$time, prob = y)),
														digits = digits, ...)
	cat("\nEstimated unconditional survival probability:\n")
	print(lapply(x$surv, function(y)cbind(time = x$time, prob = y)),
														digits = digits, ...)
	invisible(x)
}





plot.predict.mixcure = function(x, type = "l", add = FALSE,
	which = 1:nrow(x$cure), curemark = FALSE, conditional = FALSE, xlab, ylab,
	ylim = c(0, 1), lty = seq(along = which), ...)
{
	if (missing(xlab)) xlab <- "Time"
	if (missing(ylab)) ylab <- "Survival Probability"

	for (i in 1:length(which)) {
		surv = if (conditional) x$uncuresurv[[which[i]]]
							else x$surv[[which[i]]]

		if (i==1 && !add)
			plot(x$times, surv, type = type, xlab = xlab, ylab = ylab,
						ylim = ylim, lty = lty[i], ...)
		else points(x$times, surv, type = type, xlab = xlab, ylab = ylab,
						ylim = ylim, lty = lty[i], ...)

		if (curemark) abline(h = x$cure[which[i], 2])
	}
	invisible()
}

survpred = function(lfit, newdata, times, ...) UseMethod("survpred")

curepred = function(ifit, newdata, times, ...) UseMethod("curepred")





residuals.mixcure = function(object, data, type = c("WLCH", "Cox-Snell",
								"M-Cox-Snell", "Martingale", "M-Martingale",
								"M2", "M3", "D2"),
							 	type2 = c("residuals", "partial"),
							 	model = c("latency", "incidence"), ...)
{
	type = match.arg(type)
	type2 = match.arg(type2)
	model = match.arg(model)

	if (class(object)!="mixcure") stop("mixcure object is missing")
	if (is.null(data)) stop("data is missing")

	mf = model.frame(formula = object$lmodel$formula, data = data,
												na.action = na.pass)
	survtime = model.response(mf)
	ttime = survtime[, 1]
	cens = survtime[, 2]
	uncureprob = object$ifit$uncureprob
	bt = coef(object$lfit$fit)
	X = model.matrix(terms(mf), data = data)
	if (length(bt) == ncol(X)) unculp = as.vector(X %*% bt)
	else unculp = as.vector(X[, -1, drop = FALSE] %*% bt)
	mf = model.frame(formula = object$imodel$formula, data = data,
												na.action = na.pass)
	Z = model.matrix(terms(mf), data = data)[, -1, drop = FALSE]
	gm = coef(object$ifit$fit)[-1]

	if (model=="latency") {
		if (type=="WLCH") {
			if (class(object$lfit)!="latency.survreg")
				stop("WLCH residuals is only defined for parametric models")

			fit = object$lfit$fit
			n = nrow(data)

			hazard = outer(ttime, 1:n, FUN = function(x, y) {
				dval = dsurvreg(x, mean = unculp[y], scale = fit$scale,
								distribution = fit$dist, parms = fit$parms)
				sval = 1-psurvreg(x, mean = unculp[y], scale = fit$scale,
								distribution = fit$dist, parms = fit$parms)
				uncureprob[y]*dval/(uncureprob[y]*sval+1-uncureprob[y])
			})

			etime = NULL
			val = NULL
			XZ = cbind(X[, -1, drop = FALSE], Z)
			for (i in 1:n) {
				if (cens[i]==1) {
					tt = ttime[i]
					etime = c(etime, tt)
					rset = ttime >= tt
					hrset = hazard[i, rset]
					val = rbind(val, XZ[i, , drop = FALSE] -
						colSums(sweep(XZ[rset, , drop = FALSE], 1, hrset,
									  FUN = "*"))/sum(hrset))
				}
			}
			list(eventtime = etime, residuals = val)
		}
		else if (type == "Cox-Snell") {
			residuals = -log(uncureprob*object$survprob + 1-uncureprob)
			resid.uncure = -log(object$survprob)
			list(type = type, cens = cens, time = ttime, weight = object$postuncure,
				residuals = residuals, resid.uncure = resid.uncure)
		}
		else if (type == "M-Cox-Snell") {
			residuals = -log(object$survprob)
			weight = object$postuncure
			resid.dist = (function(resids, cens, weight){
				design = svydesign(ids = ~1, weights = weight[weight>0])
				svykm(Surv(resids[weight>0], cens[weight>0]) ~ 1, design = design)
			})(residuals, cens, weight)
			list(type = type, cens = cens, time = ttime, weight = weight,
				residuals = residuals, resid.dist = resid.dist)
		}
		else if (type == "M2") {
			val = cens + log(object$survprob)
			list(type = type, data = data, residuals = cens + log(object$survprob))
		}
		else if (type == "M3") {
			val = cens + object$postuncure*log(object$survprob)
			list(type = type, data = data,
				unculp = unculp,
				residuals = val,
				deviance = sign(val)*sqrt(-2*(val + cens*log(cens-val))))
		}
		else if (type == "M-Martingale") {
			keep = object$postuncure!=0
			val = cens + object$postuncure*log(object$survprob)

			if (object$lmodel$fun == "coxph")
				val2 = residuals(object$lfit$fit, weighted = TRUE)
			else val2 = NULL

			deviance = sign(val)*sqrt(-2*(val + cens*log(cens-val)))
			list(type = type, xm = X[keep, ], keep = keep,
				unculp = unculp[keep], residuals = val[keep],
				deviance = deviance[keep], residuals2 = val2)
		}
		else if (type == "Martingale") {
			val = cens + log(uncureprob*(object$survprob - 1) + 1)
			list(type = type, xm = X,
				unculp = unculp,
				residuals = val,
				deviance = sign(val)*sqrt(-2*(val + cens*log(cens-val))))
		}
		else if (type == "D2") {
			val = cens + log(object$survprob)
			val = sign(val)*sqrt(-2*(val + cens*log(cens-val)))
			list(type = type, cens = cens, time = ttime, residuals = val)
		}
		else stop("The type of residuals are not defined yet")
	}
	else {
		residuals = object$postuncure - uncureprob
		s.resid = residuals/(uncureprob*(1 - uncureprob))
		if (type2 == "residuals")
			list(residuals = residuals, s.resid = s.resid,
				resid.true = 1-data$cure - uncureprob, data = data)
		else list(Z = Z, partial = sweep(sweep(Z, 2, gm, FUN = `*`), 1, s.resid,
																FUN = `*`))
	}
}



incidence.glm = function(formula, data, ...)
{
	ff = paste("glm(formula = postuncure ", deparse(formula),
					", data = data, ...)", collapse = "")
	fit = eval(parse(text = ff))
	
	fit$call = parse(text = ff)

	val = list(uncureprob = fitted(fit), fit = fit)
    class(val) = "incidence.glm"
    val
}

coef.incidence.glm = function(x)
{
    coef(x$fit)
}

loglik = function(fit) UseMethod("loglik")

loglik.incidence.glm = function(x)
{
    x$fit$deviance
}

curepred.incidence.glm = function(ifit, newdata)
{
	predict(ifit$fit, newdata = newdata, type = "response")
}

incidence.gam = function(formula, data, ...)
{
	ff = paste("postuncure ", paste(deparse(formula), collapse = ""))
	fit = gam(formula = as.formula(ff), data = data, ...)
	val = list(uncureprob = fitted(fit), fit = fit)
    class(val) = "incidence.gam"
    val
}

coef.incidence.gam = function(x)
{
	coef(x$fit)
}

loglik.incidence.gam = function(x)
{
	x$fit$deviance
}

curepred.incidence.gam = function(ifit, newdata)
{
	predict(ifit$fit, newdata = newdata, type = "response")
}

latency.coxph = function(formula, data, ...)
{
	ff = paste(deparse(formula), collapse = "")

	if (is.null(list(...)$init))
		fit = coxph(formula = as.formula(ff),
				weights = data$postuncure[data$postuncure > 0],
				data = data[data$postuncure > 0, ], ties = 'breslow', x = TRUE, ...)
	else {
		fit = coxph(formula = as.formula(ff),
			weights = data$postuncure[data$postuncure > 0],
			init = list(...)$init, data = data[data$postuncure > 0, ],
			ties = 'breslow', x = TRUE, ...)
		if (any(is.na(coef(fit)))) 
			fit = coxph(formula = as.formula(ff),
				weights = data$postuncure[data$postuncure > 0],
				data = data[data$postuncure > 0, ], ties = 'breslow', x = TRUE, ...)
	}

	stime = model.response(model.frame(formula = as.formula(ff),
									data = data, na.action = na.pass))
	basemaxtime = max(stime[stime[, 2]==1, 1])
	stime = stime[, 1]

	if (ncol(fit$x)>0) {
		surv = survfit(fit, newdata = data, type = "kaplan-meier")

		survprob = sapply(1:nrow(data), function(x) {
			if (stime[x] > basemaxtime) 0
			else summary(surv[x], times = stime[x])$surv
		})
	}
	else {
		surv = survfit(fit, type = "kaplan-meier")$surv
		survprob = ifelse(stime > basemaxtime, 0, surv)
	}

	val = list(fit = fit, survprob = survprob, maxuncutime = basemaxtime)
	class(val) = "latency.coxph"
	val
}

coef.latency.coxph = function(x)
{
    coef(x$fit)
}

loglik.latency.coxph = function(x)
{
    x$loglik
}

survpred.latency.coxph = function(lfit, newdata, times)
{
	surv = survfit(lfit$fit, newdata = newdata, type = "kaplan-meier")

	lapply(1:nrow(newdata), function(x) {
		sapply(times, function(y) {
			if (y > lfit$maxuncutime) 0
			else summary(surv[x], times = y)$surv
		})
	})
}

latency.wcox = function(formula, data, ...)
{
	ff = paste(deparse(formula), collapse = "")

	if (is.null(list(...)$init)) fit = coxph(formula = as.formula(ff),
				weights = data$postuncure[data$postuncure > 0],
				data = data[data$postuncure > 0, ],
				ties = 'breslow', ...)
	else {
		fit = coxph(formula = as.formula(ff),
				weights = data$postuncure[data$postuncure > 0],
				init = list(...)$init, data = data[data$postuncure > 0, ],
				ties = 'breslow', ...)
		if (any(is.na(coef(fit)))) 
			fit = coxph(formula = as.formula(ff),
				weights = data$postuncure[data$postuncure > 0],
				data = data[data$postuncure > 0, ], ties = 'breslow', ...)
	}
	if (is.null(fit$weights)) fit$weights = data$postuncure[data$postuncure > 0]

	m = model.frame(formula = as.formula(ff), data = data, na.action = na.pass)
	stime = model.response(m)
	X = model.matrix(terms(m), data = data)
	x = as.matrix(X[, -1])
	BX = drop(x %*% fit$coef)




	
	od = order(fit$y[, 1], 1-fit$y[, 2])
	time = fit$y[od, 1]
	cens = fit$y[od, 2]
	bx = (fit$linear.predictors + sum(fit$coef*fit$means))[od]
	w = fit$weights[od]

	webx = w*exp(bx)
	Yeb = apply(cbind(time, cens), 1, function(x) {
		if (x[2]==1) sum((time >= x[1]) * webx)
		else NA
	})
	tmp = ifelse(cens==1, w/Yeb, 0)

	tt = unique(time[cens==1])
	maxtt = max(tt)
	Haz = sapply(tt, function(x)sum((time <= x)*tmp))

	survprob = sapply(1:length(BX), function(y) {
		if (stime[y, 1] > maxtt) 0
		else {
			survfun = stepfun(tt, c(1, exp(-exp(BX[y])*Haz)))
			survfun(stime[y, 1])
		}
	})

	val = list(fit = fit, survprob = survprob, maxuncutime = maxtt,
			   baseHaz = cbind(time = tt, Haz = Haz))
	class(val) = "latency.wcox"
	val
}

coef.latency.wcox = function(x)
{
    coef(x$fit)
}

loglik.latency.wcox = function(x)
{
    loglik(x$fit)
}

survpred.latency.wcox = function(lfit, newdata, times)
{
	fit = lfit$fit
	lp = predict(fit, newdata = newdata, type = "lp")
	lp = lp + sum(lfit$fit$coef*lfit$fit$means)
	
	maxtt = lfit$maxuncutime

	lapply(lp, function(x) {
		sapply(times, function(y) {
			if (y > maxtt) 0
			else {
				survfun = stepfun(lfit$baseHaz[, "time"],
								  c(1, exp(-exp(x)*lfit$baseHaz[, "Haz"])))
				survfun(y)
			}
		})
	})
}

latency.cox.aalen = function(formula, data, ...)
{
	ff = paste(deparse(formula), collapse = "")

	fit = cox.aalen(formula = as.formula(ff), beta = list(...)$init,
					weights = data$postuncure[data$postuncure > 0],
					data = data[data$postuncure > 0, ], ...)

	stime = model.response(model.frame(formula = as.formula(ff),
									data = data, na.action = na.pass))
	basemaxtime = max(stime[stime[, 2]==1, 1])

	surv = diag(predict(fit, newdata = data, times = stime[, 1])$S0)
	survprob = ifelse(stime[, 1] > basemaxtime, 0, surv)

	val = list(fit = fit, survprob = survprob, maxuncutime = basemaxtime)
	class(val) = "latency.cox.aalen"
	val
}

coef.latency.cox.aalen = function(x)
{
    coef(x$fit)[, 1]
}

loglik.latency.cox.aalen = function(x)
{
    x$fit$loglike[2]
}

survpred.latency.cox.aalen = function(lfit, newdata, times)
{
	surv = predict(lfit$fit, newdata = newdata, times = times)$S0

	lapply(1:nrow(newdata), function(x) {
		ifelse(times > lfit$maxuncutime, 0, surv[x, ])
	})
}

latency.prop.odds = function(formula, data, ...)
{
	ff = paste(deparse(formula), collapse = "")
	fit = prop.odds(formula = as.formula(ff), beta = list(...)$init,
					weights = data$postuncure[data$postuncure > 0],
					data = data[data$postuncure > 0, ], ...)

	mf = model.frame(formula = as.formula(ff), data = data, na.action = na.pass)
	stime = model.response(mf)
	basemaxtime = max(stime[stime[, 2]==1, 1])

	Z = model.matrix(delete.response(as.formula(ff)), mf)[, -1, drop = FALSE]

	surv = sapply(1:nrow(Z), function(x)predict(fit, Z = Z[x, ],
												times = stime[x, 1])$S0)
	survprob = ifelse(stime[, 1] > basemaxtime, 0, surv)

	val = list(fit = fit, survprob = survprob, maxuncutime = basemaxtime,
			   call = as.formula(ff))
	class(val) = "latency.prop.odds"
	val
}

coef.latency.prop.odds = function(x)
{
    coef(x$fit)[, 1]
}

loglik.latency.prop.odds = function(x)
{
    x$fit$loglike
}

survpred.latency.prop.odds = function(lfit, newdata, times)
{
	mf = model.frame(formula = lfit$call, data = newdata, na.action = na.pass)
	Z = model.matrix(delete.response(lfit$call), mf)[, -1, drop = FALSE]
	surv = lapply(1:nrow(newdata), function(x)predict(lfit$fit, Z = Z[x, ],
												times = times)$S0)

	lapply(1:nrow(newdata), function(x) {
		ifelse(times > lfit$maxuncutime, 0, surv[x, ])
	})
}

latency.survfit = function(formula, data, ...)
{
	ff = paste(deparse(formula), collapse = "")
	fit = eval(parse(text = paste("survfit(formula =", deparse(formula),
		", weights = data$postuncure[data$postuncure > 0]",
		", data = data[data$postuncure > 0, ], ...)")))

	mf = model.frame(formula = as.formula(ff), data = data, na.action = na.pass)
	respdat = model.response(mf)
	stime = respdat[, 1]

	if (is.null(fit$strata)) {
		maxtime = max(stime[respdat[, 2]==1])
		survprob = getsurv(fit, times = stime, maxtime = maxtime)
	}
	else {
		stratdat = as.data.frame(model.matrix(delete.response(as.formula(ff)),
											  mf)[, -1, drop = FALSE])

		NewStratVerb <- survival::strata(stratdat)
		NewStrat <- interaction(stratdat, sep = " ")
		levs <- levels(summary(fit)$strata)

		maxtime = sapply(1:length(fit$strata),
						 function(x)max(summary(fit[x])$time))

		if (!all(gp <- match(NewStratVerb, levs, nomatch = FALSE)) && 
			!all(gp <- match(NewStrat, levs, nomatch = FALSE))) 
			stop("Not all strata levels occur in fit.")

		survprob = sapply(1:length(gp), function(x)getsurv(fit, gp[x],
												stime[x], maxtime[gp[x]]))
	}

	val = list(fit = fit, survprob = survprob, maxtime = maxtime)
	class(val) = "latency.survfit"
	val
}

coef.latency.survfit = function(x)NULL

loglik.latency.survfit = function(x)0

getsurv = function(fit, group = NULL, times, maxtime = NULL)
{
	sapply(times, function(x) {
		if (is.null(group)) val = summary(fit, times = x)$surv
		else val = summary(fit[group], times = x)$surv
		if (length(val)==0)
			val = ifelse(is.null(group), min(fit$surv), min(fit[group]$surv))
		if (!is.null(maxtime) && x > maxtime) val = 0
		val
	})
}

survpred.latency.survfit <- function (lfit, newdata, times) 
{
	object = lfit$fit

	if (class(object) != "survfit" || object$type != "right") 
		stop("Predictions only available for class 'survfit'")

	if (is.null(object$strata)) {
		warning("There are no strata in model fitting. So newdata is ignored")
		list(getsurv(object, times = times, maxtime = lfit$maxtime))
	}
	else {
		covars <- attr(terms(eval.parent(object$call$formula)), "term.labels")

		if (!all(match(covars, names(newdata), nomatch = FALSE))) 
			stop("Not all strata defining variables occur in newdata.")

		stratdat <- newdata[, covars, drop = FALSE]

		NewStratVerb <- survival::strata(stratdat)
		NewStrat <- interaction(stratdat, sep = " ")
		levs <- levels(summary(object)$strata)

		if (!all(gp <- match(NewStratVerb, levs, nomatch = FALSE)) && 
			!all(gp <- match(NewStrat, levs, nomatch = FALSE))) 
			stop("Not all strata levels in newdata occur in fit.")

		lapply(gp, function(x)getsurv(object, x, times, lfit$maxtime[x]))
	}
}

latency.survreg = function(formula, data, ...)
{
	ff = paste("survreg(formula =", deparse(formula),
			", weights = data$postuncure[data$postuncure > 0]",
			", data = data[data$postuncure > 0, ], ...)",
			collapse = "")
	fit = eval(parse(text = ff))

	ff = paste(deparse(formula), collapse = "")
	stime = model.response(model.frame(formula = as.formula(ff),
									data = data, na.action = na.pass))[, 1]
	lp = cbind(predict(fit, newdata = data, type = "lp"), stime)

	survprob = apply(lp, 1, function(x)
		1-psurvreg(x[2], mean = x[1], scale = fit$scale, distribution = fit$dist,
															parms = fit$parms))

	val = list(fit = fit, survprob = survprob)
	class(val) = "latency.survreg"
	val
}

coef.latency.survreg = function(x)
{
    coef(x$fit)
}

loglik.latency.survreg = function(x)
{
    x$loglik[2]
}

survpred.latency.survreg = function(lfit, newdata, times)
{
	fit = lfit$fit
	lp = as.list(predict(fit, newdata = newdata, type = "lp"))

	lp = lapply(lp, function(x)c(x, times))

	lapply(lp, function(x)
		1-psurvreg(x[-1], mean = x[1], scale = fit$scale, distribution = fit$dist,
															parms = fit$parms))
}

latency.flexsurvreg = function(formula, data, ...)
{
	ff = paste(deparse(formula), collapse = "")

	fit = flexsurvreg(formula = as.formula(ff),
					  weights = data$postuncure[data$postuncure > 0],
					  data = data[data$postuncure > 0, ], ...)

	stime = model.response(model.frame(formula = as.formula(ff),
									data = data, na.action = na.pass))[, 1]
	survprob = sapply(1:nrow(data), function(x)
		unlist(summary(fit, newdata = data[x, ], type = "survival",
					   t = stime[x])[[1]][2]))

	val = list(fit = fit, survprob = survprob)
	class(val) = "latency.flexsurvreg"
	val
}

latency.flexsurvspline = function(formula, data, ...)
{
	ff = paste(deparse(formula), collapse = "")

	fit = flexsurvspline(formula = as.formula(ff),
				weights = data$postuncure[data$postuncure > 0],
				data = data[data$postuncure > 0, ], ...)

	stime = model.response(model.frame(formula = as.formula(ff),
									data = data, na.action = na.pass))[, 1]
	survprob = sapply(1:nrow(data), function(x)
		unlist(summary(fit, newdata = data[x, ], type = "survival",
					   t = stime[x])[[1]][2]))

	val = list(fit = fit, survprob = survprob)
	class(val) = "latency.flexsurvreg"
	val
}

coef.latency.flexsurvreg = function(x)
{
    x$fit$res[, 1]
}

loglik.latency.flexsurvreg = function(x)
{
    x$fit$loglik
}

survpred.latency.flexsurvreg = function(lfit, newdata, times)
{
	fit = lfit$fit
	lapply(summary(fit, newdata = newdata, t = times, type = "survival"),
		   function(x)x[, 2])
}

fcat = function(...)
{
	cat(...)
	flush.console()
}
