chip.deconv <-
function (data, where = NA, center = NA, window = 30000, fit.res = 10, 
    max.steps = 200, post.proc.factor = 2, min.npeaks = 0, max.npeaks = 99999, 
    selection.method = "bic", quant.cutoff = "q0.85", n.boot = 1, 
    boot.sample.opt = c("residual", "resample", "case", "wild", 
        "position", "replicate")[1], max.peak = NA, boot.vary.res = F, 
    kernel = NA, tile.distance = NA, verbose = T, trace = F, 
    ...) 
{
    if (!exists("type.lars")) {
        type.lars = "lasso"
        boot.max.steps.factor <- 1.2
        lars.obj.return <- F
        matrix.return <- F
        shrink <- T
        fit.bg <- NA
        plot.status <- F
        plot.boot <- F
        interp <- T
        ls.final.do <- T
        cluster.nodes <- NA
    }
    in.args <- c(mget(names(formals()), envir = as.environment(-1)), 
        sapply(as.list(substitute({
            ...
        })[-1]), deparse))
    shrink.output <- function(y) {
        y[[1]]$args$matrix.return <- y[[1]]$matrix <- y[[1]]$args$kernel <- NULL
        if (length(y) == 1) 
            return(y)
        for (j in 2:length(y)) {
            y[[j]]$args$matrix.return <- y[[j]]$args$data <- y[[j]]$fit <- y[[j]]$matrix <- NULL
            y[[j]]$args$kernel <- y[[j]]$data <- y[[j]]$kernel <- y[[j]]$all.coeffs <- y[[j]]$out.info <- NULL
            y[[j]]$args <- list(fit.res = y[[j]]$args$fit.res)
        }
        y
    }
    data <- load.chip.data(data, verbose = verbose)
    if (is.na(center) || is.na(window)) {
        center <- round(mean(data[, 1]))
        window <- diff(round(range(data[, 1])))
    }
    wind <- round(c(center - window/2, center + window/2))
    w.expand <- round(c(-0.1, 0.1) * window) + wind
    data <- data[data[, 1] >= w.expand[1], , drop = F]
    if (length(data) <= 0 || nrow(data) <= 2) 
        return(NULL)
    data <- data[data[, 1] <= w.expand[2], , drop = F]
    if (length(data) <= 0 || nrow(data) <= 2) 
        return(NULL)
    if (!is.na(where) && where %in% rownames(data)) 
        data <- data[rownames(data) == where, , drop = F]
    if (length(data) <= 0 || nrow(data) <= 2) 
        return(NULL)
    if (any(is.na(data[, 1]))) 
        data <- data[-which(is.na(data[, 1])), ]
    if (length(data) <= 0 || nrow(data) <= 2) 
        return(NULL)
    if (is.na(max.peak)) 
        max.peak <- max(data[, 2], na.rm = T) * 2
    x <- data[, 1]
    y <- data[, 2]
    if (is.na(max.steps)) {
        max.steps <- round(diff(range(x))/200)
        if (verbose) 
            cat("Using max.steps =", max.steps, "\n")
    }
    if (is.na(quant.cutoff)) {
        quant.cutoff <- 0
    }
    else if (is.character(quant.cutoff)) {
        quant.cutoff <- as.numeric(gsub("q", "", quant.cutoff))
        quant.cutoff <- quantile(data[, 2], probs = quant.cutoff, 
            na.rm = T)
    }
    else {
        quant.cutoff <- quant.cutoff
    }
    if (verbose) 
        cat("Using", quant.cutoff, "as data cutoff!\n")
    if (!any(y >= quant.cutoff)) 
        return(NULL)
    if (is.na(tile.distance)) {
        if (is.null(rownames(data))) 
            rownames(data) <- rep("X1", nrow(data))
        tmp.posns <- sort(data[, 1])
        tile.distances <- numeric()
        for (i in unique(names(tmp.posns))) tile.distances <- c(tile.distances, 
            diff(tmp.posns[names(tmp.posns) == i]))
        tile.distance <- median(tile.distances[tile.distances > 
            1])
        if (verbose) 
            cat("MEAN PROBE SPACING =", tile.distance, "\n")
        rm(tmp.posns, tile.distances)
    }
    if (n.boot < 1) 
        n.boot <- 1
    ord <- order(x)
    x <- x[ord]
    y <- y[ord]
    orig.x <- orig.orig.x <- x
    orig.y <- orig.orig.y <- y
    posns <- round(orig.x - min(orig.x) + 1)
    cnts <- as.numeric(orig.y)
    orig.fit.res <- fit.res
    best.step.1 <- NA
    all.out <- list()
    iter.boot <- 1
    while (iter.boot <= n.boot) {
        is.boot <- FALSE
        if (iter.boot > 1) {
            is.boot <- TRUE
            if (verbose) 
                cat("*** BOOTSTRAP ITER:", iter.boot, "***\n")
            orig.x <- orig.orig.x
            orig.y <- orig.orig.y
            if (!is.na(pmatch("wild", boot.sample.opt))) {
                if (trace) 
                  cat("Resampling: wild!\n")
                if (length(all.out) > 0) 
                  fit <- all.out[[1]]$fit
                else fit <- matrix(c(0, 0), ncol = 2)
                dens <- density(orig.y, na.rm = T)
                level.subt <- dens$x[which.max(dens$y)]
                resids <- orig.y - level.subt - fit[, 2]
                if ("wild.1" %in% boot.sample.opt) 
                  resids <- resids * sample(c(-(sqrt(5) - 1)/2, 
                    (sqrt(5) + 1)/2), length(resids), prob = c((sqrt(5) - 
                    1)/(2 * sqrt(5)), (sqrt(5) + 1)/(2 * sqrt(5))), 
                    replace = T)
                else resids <- resids * sample(c(1, -1), length(resids), 
                  replace = T)
                orig.x <- fit[, 1]
                orig.y <- fit[, 2] + resids
            }
            if (!is.na(pmatch("residual", boot.sample.opt))) {
                if (trace) 
                  cat("Resampling: residual!\n")
                if (length(all.out) > 0) 
                  fit <- all.out[[1]]$fit
                else fit <- matrix(c(0, 0), ncol = 2)
                dens <- density(orig.y, na.rm = T)
                level.subt <- dens$x[which.max(dens$y)]
                resids <- orig.y - level.subt - fit[, 2]
                if ("residual.1" %in% boot.sample.opt) 
                  resids <- resids * sample(c(1, -1), length(resids), 
                    replace = T)
                else resids <- resids * sample(c(-(sqrt(5) - 
                  1)/2, (sqrt(5) + 1)/2), length(resids), prob = c((sqrt(5) - 
                  1)/(2 * sqrt(5)), (sqrt(5) + 1)/(2 * sqrt(5))), 
                  replace = T)
                orig.y <- resids
                quant.cutoff <- median(orig.y)
            }
            if (!is.na(pmatch("replicate", boot.sample.opt))) {
                if (trace) 
                  cat("Resampling: replicate!\n")
                sds <- tapply(orig.y, orig.x, sd, na.rm = T)
                mns <- tapply(orig.y, orig.x, mean, na.rm = T)
                unq <- unique(orig.x)
                new.data <- t(apply(cbind(orig.x, orig.y), 1, 
                  function(i) {
                    ind <- which(unq == i[1])
                    c(i[1], rnorm(1, mean = mns[ind], sd = sds[ind]))
                  }))
                orig.x <- new.data[, 1]
                orig.y <- new.data[, 2]
            }
            if (!is.na(pmatch("resample", boot.sample.opt))) {
                if (trace) 
                  cat("Resampling: resample!\n")
                orig.y <- sample(orig.y, replace = T)
            }
            if (!is.na(pmatch("case", boot.sample.opt))) {
                if (trace) 
                  cat("Resampling: case!\n")
                samp <- sort(sample(1:length(orig.orig.y), replace = T))
                if (!1 %in% samp) 
                  samp <- c(1, samp)
                orig.x <- orig.x[samp]
                orig.y <- orig.y[samp]
            }
            if (!is.na(pmatch("position", boot.sample.opt))) {
                if (trace) 
                  cat("Resampling: position!\n")
                offsets <- round(rnorm(length(orig.x), mean = 0, 
                  sd = tile.distance/4))
                sgn <- sign(offsets)
                offsets[abs(offsets) > tile.distance - 5] <- (tile.distance - 
                  5) * sgn[abs(offsets) > tile.distance - 5]
                orig.x <- orig.x + offsets
                ord <- order(orig.x)
                orig.y <- orig.y[ord]
                orig.x <- orig.x[ord]
            }
            if (!is.na(boot.vary.res)) {
                if (is.numeric(boot.vary.res)) {
                  fit.res <- sample(round(orig.fit.res - boot.vary.res):round(orig.fit.res + 
                    boot.vary.res), 1)
                }
                else if (is.logical(boot.vary.res) && boot.vary.res == 
                  TRUE) {
                  fit.res <- sample(round(orig.fit.res * 4/5):round(orig.fit.res * 
                    6/5), 1)
                }
                if (trace && fit.res != orig.fit.res) 
                  cat("Varying the boot resolution... current fit.res =", 
                    fit.res, "\n")
            }
        }
        posns <- round(orig.x - min(orig.x) + 1)
        cnts <- as.numeric(orig.y)
        dens <- density(cnts, na.rm = T)
        level.subtract <- dens$x[which.max(dens$y)]
        cnts <- cnts - level.subtract
        good.locs <- which(cnts > quant.cutoff - level.subtract)
        if (!exists("xxx")) 
            xxx <- NULL
        if (length(good.locs) <= 0 || length(unique(posns)) < 
            3) {
            if (n.boot <= 1) 
                return(NULL)
            out.info <- c(best.step = 0, n.coeffs = 0, n.coeffs.nr = 0, 
                dof = 0, sigma.e.sq = NA, bic = NA, aic = NA, 
                tile.distance = tile.distance, rss = NA)
            out <- list(data = cbind(orig.x, orig.y - level.subtract), 
                args = in.args, fit = cbind(orig.x, rep(0, length(orig.x))), 
                window = wind, kernel = kernel, coeffs = matrix(nrow = 0, 
                  ncol = 2, dimnames = list(c(), c("position", 
                    "intensity"))), out.info = out.info)
            attr(out, "class") <- "chip.deconv"
            if (!(length(matrix.return) == 1 && matrix.return == 
                FALSE) && exists("xxx")) 
                out$matrix <- xxx
            if (verbose) 
                cat("Number of coeffs:", 0, "...\n")
            if (n.boot <= 1 || iter.boot == 1) {
                if (n.boot <= 1) 
                  all.out <- out
                else {
                  for (i in 1:n.boot) all.out[[i]] <- out
                  if (shrink) 
                    all.out <- shrink.output(all.out)
                }
                attr(all.out, "class") <- "chip.deconv"
                return(invisible(all.out))
            }
            else {
                all.out[[length(all.out) + 1]] <- out
                iter.boot <- iter.boot + 1
                next
            }
        }
        try(if (interp && tile.distance > 300 && (!is.boot || 
            is.na(pmatch("position", boot.sample.opt)))) {
            require(zoo, quietly = T, warn.conflicts = F)
            orig.posns <- posns
            orig.cnts <- cnts
            tile.distances <- diff(posns)
            needs.fill <- tile.distances > 1.8 * tile.distance
            if (any(needs.fill)) {
                cnts.filled <- cnts[1]
                posns.filled <- posns[1]
                for (i in 2:length(posns)) {
                  if (!needs.fill[i - 1]) {
                    posns.filled <- c(posns.filled, posns[i])
                    cnts.filled <- c(cnts.filled, cnts[i])
                  }
                  else {
                    posns.filled <- c(posns.filled, rep(NA, round((posns[i] - 
                      posns[i - 1])/tile.distance) - 1), posns[i])
                    cnts.filled <- c(cnts.filled, rep(NA, round((posns[i] - 
                      posns[i - 1])/tile.distance) - 1), cnts[i])
                  }
                }
                posns.filled <- na.approx(posns.filled)
                if (tile.distance > 250) 
                  cnts.filled[is.na(cnts.filled)] <- 0
                else cnts.filled <- na.approx(cnts.filled, posns.filled)
                posns <- posns.filled
                cnts <- cnts.filled
                orig.x <- posns + min(orig.x) - 1
                orig.y <- cnts + level.subtract
                rm(cnts.filled, posns.filled)
            }
        }, silent = T)
        good.posns <- unique(round(sort(posns[good.locs])))
        if (!is.null(matrix.return) && class(matrix.return) %in% 
            c("matrix", "dgCMatrix")) 
            xxx <- matrix.return
        tmp <- round((-max(tile.distance, fit.res) - 1):(max(tile.distance, 
            fit.res) + 1))
        good.posns.hires <- unique(as.vector(sapply(good.posns, 
            function(i) i + tmp)))
        good.posns.hires <- good.posns.hires[good.posns.hires %in% 
            seq(1, diff(range(posns)) + fit.res * 10, by = fit.res)]
        if (plot.status) {
            if (iter.boot <= 1) 
                par(mfrow = c(2, 1))
            plot(orig.x, cnts, pch = 19, cex = 0.2)
            points(orig.x[good.locs], cnts[good.locs], pch = 19, 
                cex = 0.5, col = "red")
            points(good.posns.hires + min(orig.x), rep(min(cnts[good.locs], 
                na.rm = T), length(good.posns.hires)), pch = 19, 
                cex = 0.3, col = "blue")
        }
        if (!exists("xxx") || is.null(xxx) || any(!as.character(good.posns.hires) %in% 
            colnames(xxx)) || any(!as.character(round(posns)) %in% 
            rownames(xxx))) {
            tmp.posns <- round(posns)
            tmp.posns <- unique(tmp.posns)
            tmp.good.posns <- good.posns.hires
            if (iter.boot > 1 && !is.na(pmatch("position", boot.sample.opt))) 
                xxx <- NULL
            if (exists("xxx") && !is.null(xxx)) {
                tmp.good.posns <- tmp.good.posns[!as.character(tmp.good.posns) %in% 
                  colnames(xxx)]
                tmp.posns <- as.numeric(rownames(xxx))
            }
            if (length(tmp.good.posns) > 0) {
                if (is.null(kernel) || is.na(kernel[1])) 
                  stop("No kernel provided!\n")
                kernel[, 2] <- kernel[, 2]/max(kernel[, 2])
                xxx.tmp <- try(make.predictor.matrix(tmp.posns, 
                  kernel, fit.res = fit.res, good.posns.hires = tmp.good.posns, 
                  sparse = T, verbose = trace), silent = !trace)
                if (class(xxx.tmp) == "try-error") 
                  stop("Could not allocate predictor matrix!!!")
                colnames(xxx.tmp)[1:length(tmp.good.posns)] <- tmp.good.posns
                rownames(xxx.tmp) <- tmp.posns
                if (exists("xxx") && !is.null(xxx) && nrow(xxx) == 
                  nrow(xxx.tmp)) {
                  if (exists("cBind")) 
                    xxx <- cBind(xxx, xxx.tmp)
                  else xxx <- cbind(xxx, xxx.tmp)
                }
                else {
                  xxx <- xxx.tmp
                }
            }
        }
        xx <- xxx
        if (ncol(xx) != length(good.posns.hires)) {
            if (!all(as.character(round(good.posns.hires)) %in% 
                colnames(xx))) {
                if (trace) 
                  cat("ERROR1\n")
                next
            }
            xx <- xx[, as.character(round(good.posns.hires)), 
                drop = F]
        }
        if (nrow(xx) != length(posns)) {
            if (!all(as.character(round(posns)) %in% rownames(xx))) {
                if (trace) 
                  cat("ERROR2\n")
                next
            }
            xx <- xx[as.character(round(posns)), , drop = F]
        }
        colnames(xx) <- as.character(1:ncol(xx))
        if (!is.na(best.step.1)) 
            max.steps <- max(20, ceiling(best.step.1 * boot.max.steps.factor))
        lrs <- try(lars.pos(xx, cnts, type = type.lars, max.steps = max.steps, 
            use.Gram = F, positive = T, trace = trace, ...), 
            silent = !trace)
        if (class(lrs)[1] == "try-error") {
            if (trace) 
                cat("ERROR3\n")
            iter.boot <- iter.boot + 1
            next
        }
        max.steps <- min(max.steps, length(lrs$actions))
        lrs.coeff <- predict(lrs, type = "coeff")
        coeff.cutoff <- 0
        n.coeffs <- apply(lrs.coeff$coefficients, 1, function(i) sum(i > 
            coeff.cutoff)) + 1
        dof <- sapply(1:length(n.coeffs), function(i) max(which(n.coeffs == 
            n.coeffs[i])))
        if (is.character(selection.method)) {
            N.x <- length(cnts)
            if (!is.null(lrs$RSS)) {
                rss <- lrs$RSS
            }
            else {
                rss <- apply(lrs.coeff$coefficients, 1, function(i) sum((xx %*% 
                  i - cnts)^2, na.rm = T))
            }
            max.row <- which.min(rss)
            hb.fit <- (xx %*% lrs.coeff$coefficients[max.row, 
                ])[, 1]
            sigma.e.sq <- mean((hb.fit - cnts)^2, na.rm = T)
            bic <- rss/sigma.e.sq + log(N.x) * dof
            aic <- rss/sigma.e.sq + 2 * dof
            if (verbose) 
                cat("Step for min AIC:", which.min(aic), n.coeffs[which.min(aic)], 
                  "; BIC:", which.min(bic), n.coeffs[which.min(bic)], 
                  "; using:", selection.method, "\n")
            if (plot.status) {
                par(mfrow = c(2, 1))
                plot(rss, typ = "l")
                plot(log(aic), typ = "l")
                plot(log(bic), typ = "l")
            }
            aic[is.na(aic)] <- bic[is.na(bic)] <- Inf
            if (selection.method == "aic.and.bic") 
                best.step <- ceiling(mean(c(which.min(aic), which.min(bic)), 
                  na.rm = T))
            else if (selection.method == "bic") 
                best.step <- which.min(bic)
            else if (selection.method == "aic") 
                best.step <- which.min(aic)
            if (iter.boot > 1 && best.step >= max.steps) 
                warning(paste("max.steps is probably too low for model selection by", 
                  selection.method))
        }
        if (is.na(best.step) || best.step <= 1 && is.numeric(selection.method)) {
            if (selection.method %in% n.coeffs) 
                best.step <- min(which(n.coeffs == selection.method))
            else best.step <- max(which(n.coeffs <= selection.method)) - 
                1
        }
        if (best.step < 1) 
            best.step <- 1
        if (n.coeffs[best.step] < min.npeaks + 1 && any(n.coeffs >= 
            min.npeaks + 1)) 
            best.step <- min(which(n.coeffs >= min.npeaks + 1))
        if (n.coeffs[best.step] > max.npeaks + 1 && any(n.coeffs <= 
            max.npeaks + 1)) 
            best.step <- max(which(n.coeffs <= max.npeaks + 1))
        if (is.na(best.step) || best.step < 1) 
            best.step <- 1
        if (is.na(best.step.1)) 
            best.step.1 <- best.step
        if (!is.na(best.step.1) && best.step >= max.steps) 
            best.step.1 <- max(20, ceiling(best.step * boot.max.steps.factor))
        coeffs <- rep(0, ncol(xx))
        names(coeffs) <- colnames(lrs.coeff$coefficients)
        coeffs[colnames(lrs.coeff$coefficients)] <- lrs.coeff$coefficients[best.step, 
            , drop = F]
        if (verbose) 
            cat(sum(coeffs > coeff.cutoff), "coeffs at >", coeff.cutoff, 
                "for LARS step", best.step, "\n")
        good.coeffs <- coeffs[coeffs > coeff.cutoff]
        good.xx <- good.posns.hires[coeffs > coeff.cutoff] + 
            min(orig.x)
        tmp.fit <- (xx %*% coeffs)[, 1]
        if (!is.na(post.proc.factor)) {
            if (verbose) 
                cat("Number of coeffs:", sum(coeffs > 0), "... ")
            tmp.coeffs <- post.proc.coeffs(cbind(seq(along = coeffs), 
                coeffs), fit.res = 1, factor = post.proc.factor, 
                max.coef = max.peak, mean.do = F)
            tmp.coeffs2 <- coeffs * 0
            tmp.coeffs2[round(tmp.coeffs[, 1])] <- tmp.coeffs[, 
                2]
            coeffs <- tmp.coeffs2
            good.coeffs <- coeffs[coeffs > coeff.cutoff]
            good.xx <- good.posns.hires[coeffs > coeff.cutoff] + 
                min(orig.x)
            if (verbose) 
                cat("Reduced to", sum(coeffs > coeff.cutoff), 
                  "non-redundant coeffs.\n")
        }
        out.info <- c(best.step = best.step, n.coeffs = n.coeffs[best.step], 
            n.coeffs.nr = sum(coeffs > coeff.cutoff), dof = dof[best.step], 
            sigma.e.sq = sigma.e.sq, bic = bic[best.step], aic = aic[best.step], 
            tile.distance = tile.distance)
        final.rss <- sum(cnts^2, na.rm = T)
        if (ls.final.do && require(quadprog, quietly = T, warn.conflicts = F) && 
            best.step > 1 && n.coeffs[best.step] > 0 && sum(coeffs > 
            coeff.cutoff) > 0) {
            tmp.xx <- as.matrix(xx[, coeffs > coeff.cutoff, drop = F])
            Dmat <- t(tmp.xx) %*% tmp.xx
            if (require(corpcor, quietly = T, warn.conflicts = F) && 
                !is.positive.definite(Dmat)) 
                Dmat <- make.positive.definite(Dmat)
            coeffs.ls <- try(solve.QP(Dmat, t(t(cnts) %*% tmp.xx), 
                Amat = t(diag(ncol(tmp.xx))), bvec = rep(0, ncol(tmp.xx)), 
                meq = 0), silent = !trace)
            if (class(coeffs.ls) == "try-error") {
                if (trace) 
                  cat("ERROR4\n")
                iter.boot <- iter.boot + 1
                next
            }
            good.coeffs <- coeffs.ls$solution
            tmp.fit <- tmp.xx %*% good.coeffs
            coeffs[coeffs > coeff.cutoff] <- good.coeffs
            names(good.coeffs) <- colnames(tmp.xx)
            final.rss <- sum((tmp.fit - cnts)^2, na.rm = T)
        }
        out.info <- c(out.info, rss = final.rss)
        names(out.info) <- c("best.step", "n.coeffs", "n.coeffs.nr", 
            "dof", "sigma.e.sq", "bic", "aic", "tile.distance", 
            "rss")
        if (plot.status) {
            par(mfrow = c(2, 1))
            plot(orig.x, cnts, pch = 19, cex = 0.5, xlim = wind, 
                ylim = range(c(cnts, tmp.fit, coeffs)), xlab = "Genome coord.", 
                ylab = "Chip intensity")
            lines(orig.x, tmp.fit, col = "red")
            apply(cbind(good.xx, good.coeffs), 1, function(i) lines(rep(i[1], 
                2), c(0, i[2]), col = "darkgreen"))
            points(good.xx, good.coeffs[1:length(good.xx)], col = "darkgreen", 
                pch = 19, cex = 0.5)
        }
        out <- list(data = cbind(orig.x, orig.y - level.subtract), 
            fit = cbind(orig.x, tmp.fit), window = wind, kernel = kernel, 
            coeffs = cbind(position = good.xx, intensity = good.coeffs), 
            out.info = out.info)
        attr(out, "class") <- "chip.deconv"
        if (lars.obj.return) 
            out$lars.out <- list(lars.obj = lrs, lars.coeff = lrs.coeff)
        out$args <- in.args
        if (!(length(matrix.return) == 1 && matrix.return == 
            FALSE) && exists("xxx")) 
            out$matrix <- xxx
        try(rm(xx, lrs, lrs.coeff), silent = T)
        all.out[[length(all.out) + 1]] <- out
        if (plot.boot && n.boot > 1) 
            xqz <- try(plot(out, hi.res = NA, main = paste("BOOT ITER:", 
                iter.boot), ...))
        if (n.boot <= 1) 
            all.out <- out
        iter.boot <- iter.boot + 1
    }
    if (n.boot > 1) {
        if ((!is.na(pmatch("resample", boot.sample.opt)) || !is.na(pmatch("residual", 
            boot.sample.opt))) && length(all.out) > 0) {
            co <- all.out[[1]]$coeffs
            co.p <- apply(co, 1, function(i) sum(sapply(all.out, 
                function(j) any(j$coeffs[, 2] >= i[2])))/length(all.out))
            co <- cbind(co, p.value = co.p)
            colnames(co)[3] <- "p.value <"
            all.out[[1]]$coeffs.w.p.values <- co
        }
        if (shrink) 
            all.out <- shrink.output(all.out)
    }
    attr(all.out, "class") <- "chip.deconv"
    if (verbose && trace) 
        print(all.out)
    invisible(all.out)
}
convolve.func <-
function (loc, kernel, positions) 
{
    out <- rep(0, diff(range(positions)) + 1)
    k.length <- nrow(kernel)
    inds <- (loc - round(k.length/2)):(loc + round(k.length/2))
    tmp <- inds > 0 & inds <= length(out)
    out[inds[tmp]] <- kernel[(1:k.length)[tmp], 2]
    out
}
deconv.entire.genome <-
function (data, chroms = NA, window = 6000, step.by = 5500, centers = NULL, 
    quiet = F, plot = F, n.boot = 1, fit.res = 10, quant.cutoff = "q0.85", 
    max.peak = NA, verbose = F, kernel = NA, no.multicore = T, 
    ...) 
{
    if (!exists("fits.fin.only")) {
        fits.fin.only <- T
        in.fits <- NULL
        remove.near.edges <- T
        cluster.nodes <- NA
        save.progress <- NULL
    }
    start.time <- Sys.time()
    created.here <- FALSE
    data <- load.chip.data(data, verbose = verbose)
    fname <- attr(data, "filename")
    posns <- sort(data[, 1])
    if (any(is.na(data[, 1]))) 
        data <- data[-which(is.na(data[, 1])), ]
    if (length(data) <= 0 || nrow(data) <= 2) 
        return(NULL)
    orig.data <- data
    if (is.na(chroms) && !is.null(rownames(orig.data))) 
        chroms <- unique(rownames(orig.data))
    if (is.na(chroms)) 
        chroms <- 1
    chroms <- as.character(chroms)
    cat("Running on chromosomes:", chroms, "\n")
    if (is.na(quant.cutoff)) {
        quant.cutoff <- 0
    }
    else if (is.character(quant.cutoff)) {
        quant.cutoff <- as.numeric(gsub("q", "", quant.cutoff))
        quant.cutoff <- quantile(orig.data[, 2], probs = quant.cutoff, 
            na.rm = T)
    }
    else {
        quant.cutoff <- quant.cutoff
    }
    cat("Using", quant.cutoff, "as data cutoff!\n")
    if (is.null(kernel) || is.na(kernel[1])) 
        stop("No kernel provided!\n")
    kernel[, 2] <- kernel[, 2]/max(kernel[, 2])
    if (is.na(max.peak)) 
        max.peak <- max(orig.data[, 2]) * 2
    if (is.null(names(posns))) 
        names(posns) <- rep("XXX", length(posns))
    tmp.posns <- sort(posns)
    tile.distances <- numeric()
    for (i in unique(names(tmp.posns))) tile.distances <- c(tile.distances, 
        diff(tmp.posns[names(tmp.posns) == i]))
    tile.distance <- median(tile.distances[tile.distances > 1])
    cat("MEAN PROBE SPACING =", tile.distance, "\n")
    rm(tmp.posns, tile.distances)
    fits <- list()
    if (!is.null(in.fits)) 
        fits <- in.fits
    for (where in chroms) {
        if (!is.null(rownames(orig.data))) 
            data <- orig.data[rownames(orig.data) == where, , 
                drop = F]
        else data <- orig.data
        rng <- range(data[, 1])
        posns <- round(seq(rng[1] + round(window/2), rng[2], 
            by = step.by))
        if (!is.null(centers)) {
            if (!is.null(names(centers))) 
                posns <- centers[which(names(centers) == where)]
            else posns <- centers
            posns <- sort(posns)
        }
        if (plot) 
            par(mfrow = c(2, 1))
        apply.func <- lapply
        is.parallel <- !no.multicore
        if (is.parallel) 
            is.parallel <- require(multicore, quietly = T, warn.conflicts = F)
        if (is.parallel) 
            is.parallel <- !multicore:::isChild()
        if (is.parallel) 
            is.parallel <- multicore:::detectCores(all.tests = TRUE) > 
                1
        if (is.parallel) {
            require(doMC)
            registerDoMC()
            apply.func <- function(list, FUN, ...) foreach(l = list) %dopar% 
                {
                  FUN(l, ...)
                }
        }
        if (is.parallel && !quiet) 
            cat("Parallelizing deconvoluton of", where, "over", 
                multicore:::detectCores(all.tests = TRUE), "processor cores.\n")
        if (is.parallel) 
            warning("WARNING: If you are running on a Windows system, the 'multicore' option will not work.\nPlease re-start with parameter 'no.multicore=TRUE'.\n")
        out.fits <- apply.func(posns, function(pos) {
            tmp.fit <- chip.deconv(data = data, window = window, 
                center = pos, where = NA, tile.distance = tile.distance, 
                kernel = kernel, max.peak = max.peak, quant.cutoff = quant.cutoff, 
                verbose = verbose, fit.res = fit.res, n.boot = n.boot, 
                ...)
            if (is.null(tmp.fit) || length(tmp.fit) == 0 || class(tmp.fit) == 
                "try-error" || (is.list(tmp.fit) && is.character(tmp.fit[[1]]) && 
                substr(tmp.fit[[1]], 1, 5) == "Error")) {
                if (verbose && !quiet) 
                  cat("ERROR running chip.deconv() on", fname, 
                    where, which(posns == pos), length(posns), 
                    pos, posns[length(posns)], "; DON'T PANIC - probably no data in the window...", 
                    "skipping to next.\n")
                return(NULL)
            }
            if (!is.null(centers)) {
                ttmp <- which(posns >= pos - window & posns <= 
                  pos + window)
                if (length(ttmp) > 0) 
                  posns <- posns[!1:length(posns) %in% ttmp]
            }
            if (plot) 
                try(plot(tmp.fit, main = paste(fname, "chip", 
                  where, ind, pos)))
            if (!is.null(tmp.fit$coeffs)) 
                tmp.fit <- list(`1` = tmp.fit)
            ddata <- tmp.fit[[1]]$data
            for (i in 1:length(tmp.fit)) {
                q.fit <- tmp.fit[[i]]
                q.fit$args$matrix.return <- q.fit$args$data <- q.fit$fit <- q.fit$matrix <- q.fit$args$kernel <- NULL
                q.fit$args$where <- where
                if (i > 1) {
                  q.fit$data <- q.fit$kernel <- q.fit$all.coeffs <- q.fit$out.info <- NULL
                  q.fit$args <- list(fit.res = q.fit$args$fit.res)
                }
                if (nrow(q.fit$coeffs) > 0 && remove.near.edges) {
                  coe <- q.fit$coeffs
                  tmp <- coe[, 1] <= min(ddata[, 1] + window/20) | 
                    coe[, 1] >= max(ddata[, 1] - window/20)
                  if (any(tmp)) 
                    coe[tmp, 2] <- 0
                  coe <- coe[coe[, 2] > 0, , drop = F]
                  q.fit$coeffs <- coe
                }
                tmp.fit[[i]] <- q.fit
            }
            if (!quiet) 
                cat(fname, where, which(posns == pos), length(posns), 
                  pos, posns[length(posns)], "\t", "BEST.STEP =", 
                  tmp.fit[[1]]$out.info["best.step"], "COEFFS =", 
                  sum(tmp.fit[[1]]$coeffs[, 2] > 0), "\n")
            return(tmp.fit)
        })
        gc()
        fits[[where]] <- list()
        for (i in 1:length(out.fits)) if (!is.null(out.fits[[i]])) 
            fits[[where]][[length(fits[[where]]) + 1]] <- out.fits[[i]]
    }
    fits.fin <- list()
    for (where in names(fits)) {
        fits.fin[[where]] <- list()
        fts <- fits[[where]]
        if (length(fts) <= 0) 
            next
        fit <- fts[[1]]
        if (!is.null(fit$coeffs)) 
            fit = list(`1` = fit)
        for (j in 1:length(fit)) {
            if (length(fts) > 1) {
                for (i in 2:length(fts)) {
                  if (j > length(fts[[i]])) 
                    next
                  for (n in c("data", "fit", "coeffs", "all.coeffs")) if (!is.null(fts[[i]][[j]][[n]])) 
                    fit[[j]][[n]] <- rbind(fit[[j]][[n]], fts[[i]][[j]][[n]])
                }
            }
            for (n in c("data", "fit", "coeffs", "all.coeffs")) {
                if (is.null(fit[[j]][[n]])) 
                  next
                if (is.vector(fit[[j]][[n]])) 
                  fit[[j]][[n]] <- matrix(fit[[j]][[n]], ncol = 2)
                ord <- order(fit[[j]][[n]][, 1])
                fit[[j]][[n]] <- fit[[j]][[n]][ord, , drop = F]
                fit[[j]][[n]] <- fit[[j]][[n]][!is.na(fit[[j]][[n]][, 
                  1]), , drop = F]
            }
            fit[[j]]$window <- range(fit[[j]]$data[, 1])
            if (!quiet) 
                cat("\t", where, j, nrow(fit[[j]]$coeffs), "REDUNDANT COEFFS\n")
            attr(fit, "class") <- "chip.deconv"
            fits.fin[[where]][[j]] <- post.proc.deconv(fit[[j]], 
                factor = 6, fit.res = fit.res, max.coef = max.peak, 
                mean.do = T)
            colnames(fits.fin[[where]][[j]]$coeffs) <- c("position", 
                "intensity")
            attr(fits.fin[[where]][[j]], "class") <- "chip.deconv"
            cat("\t\t", where, j, nrow(fits.fin[[where]][[j]]$coeffs), 
                "NON-REDUNDANT COEFFS\n")
        }
        attr(fits.fin[[where]], "class") <- "chip.deconv"
    }
    if (fits.fin[[1]][[1]]$args$boot.sample.opt %in% c("resample", 
        "residual")) {
        all.coeffs.boot <- do.call(rbind, lapply(fits.fin, function(i) do.call(rbind, 
            lapply(i[2:length(i)], function(j) rbind(j$coeffs, 
                c(0, 0))))))
        for (i in 1:length(fits.fin)) {
            if (is.null(fits.fin[[i]]) || length(fits.fin[[i]]) <= 
                0) 
                next
            coeffs <- fits.fin[[i]][[1]]$coeffs
            pvs <- apply(coeffs, 1, function(j) (sum(all.coeffs.boot[, 
                2] >= j[2]) + 1)/(nrow(all.coeffs.boot) + 1))
            tmp <- cbind(coeffs, p.value = pvs)
            colnames(tmp)[3] <- "p.value <"
            fits.fin[[i]][[1]]$coeffs.w.p.values <- tmp
        }
    }
    end.time <- Sys.time()
    cat("Completed in", difftime(end.time, start.time, units = "mins"), 
        "minutes.\n")
    out <- list(fits.fin = fits.fin)
    if (!fits.fin.only) 
        out$fits <- fits
    attr(out, "class") <- "chip.deconv.entire.genome"
    invisible(out)
}
fit.peak.profile <-
function (data, tile.size, n.peaks = 30, n.skip = 10, in.kernel = NA, 
    fits = NULL, method = "Nelder-Mead", positions = c(0, 25, 
        50, 100, 150, seq(200, (mini.window * 1.5) + 1, by = 100)), 
    re.fit = 25, start.pars = c(shape = 7, scale = 50, bs.size = 20, 
        h.cutoff = 15), to.be.fit = c("shape", "scale", "bs.size", 
        "h.cutoff"), rnd = F, mini.window = max(5 * tile.size, 
        1300), plot = T, name = "", no.multicore = T, ...) 
{
    default.start.pars <- c(shape = 7, scale = 50, bs.size = 20, 
        h.cutoff = 15, offset = 0)
    default.start.pars[names(start.pars)] <- start.pars
    in.args <- c(mget(names(formals()), envir = as.environment(-1)), 
        sapply(as.list(substitute({
            ...
        })[-1]), deparse))
    data <- load.chip.data(data, verbose = T)
    orig.data <- data
    best.score.so.far <- 9e+09
    iter <- 1
    best.kernel <- best.params <- NULL
    apply.func <- lapply
    is.parallel <- !no.multicore
    if (is.parallel) 
        is.parallel <- require(multicore, quietly = T, warn.conflicts = F)
    if (is.parallel) 
        is.parallel <- !multicore:::isChild()
    if (is.parallel) 
        is.parallel <- multicore:::detectCores(all.tests = TRUE) > 
            1
    if (is.parallel) 
        apply.func <- mclapply
    if (is.parallel) 
        cat("Parallelizing profile fitting over", multicore:::detectCores(all.tests = TRUE), 
            "processor cores.\n")
    if (is.parallel) 
        warning("WARNING: If you are running on a Windows system, the 'multicore' option will not work.\nPlease re-start with parameter 'no.multicore=TRUE'.\n")
    get.profile <- function(par, ...) {
        par[names(default.start.pars)[!names(default.start.pars) %in% 
            names(par)]] <- default.start.pars[!names(default.start.pars) %in% 
            names(par)]
        par[c("shape", "scale")] <- abs(par[c("shape", "scale")]) + 
            0.01
        hc <- round(abs(par["h.cutoff"]))
        hc <- min(tile.size - 1, hc)
        bs.size <- round(abs(par["bs.size"]))
        offset <- par["offset"]
        kernel <- generate.binding.profile(tile.size = tile.size, 
            interp = T, bs.size = bs.size, plot = F, hybridization.prob = function(x, 
                ...) as.integer(x > hc), fragment.distrib = function(x, 
                ...) dgamma(x - offset, shape = par["shape"], 
                scale = par["scale"]), positions = positions, 
            verbose = F, no.multicore = no.multicore, ...)
        if (all(is.na(kernel[, 2]) | is.infinite(kernel[, 2]) | 
            kernel[, 2] == 0 | kernel[, 2] == 1)) 
            return(9e+09)
        kernel
    }
    pks <- datas <- NULL
    changed <- FALSE
    get.fit.score <- function(par, pks, kernel = NULL, plot = F, 
        return.all = F) {
        if (changed && !is.na(re.fit) && iter%%re.fit == 0) {
            cat(iter, "--> Re-fitting current best profile to all data to find biggest peaks...\n")
            tmp.fits <- deconv.entire.genome(data = orig.data, 
                kernel = best.kernel, plot = F, verbose = F, no.multicore = no.multicore, 
                ...)
            ppks <- get.biggest.peaks(tmp.fits$fits.fin, n.peaks = n.peaks)
            fits <<- fits
            pks <<- ppks
            datas <<- list()
            changed <<- FALSE
            for (i in 1:nrow(pks)) datas[[i]] <<- data[rownames(data) == 
                rownames(pks)[i] & data[, 1] >= pks[i, 1] - mini.window * 
                1.1 & data[, 1] <= pks[i, 1] + mini.window * 
                1.1, , drop = F]
        }
        par[names(default.start.pars)[!names(default.start.pars) %in% 
            names(par)]] <- default.start.pars[!names(default.start.pars) %in% 
            names(par)]
        cat(iter, par, " ")
        if (is.null(kernel) || class(kernel) != "matrix") 
            kernel <- get.profile(par, ...)
        new.fits <- apply.func(1:nrow(pks), function(i) {
            if (i > length(datas)) 
                return(NULL)
            if (is.null(datas[[i]])) 
                return(NULL)
            dat <- datas[[i]]
            chip.deconv(data = dat, fit.res = 30, window = NA, 
                max.steps = max.steps.2, where = NA, center = NA, 
                kernel = kernel, verbose = F, trace = F, n.boot = 1, 
                tile.distance = tile.size, max.npeaks = 1, min.npeaks = 1, 
                quant.cutoff = q.cutoff)
        })
        pks <- t(sapply(new.fits, function(i) {
            tmp <- i$coeffs[which.max(i$coeffs[, 2]), ]
            return(if (length(tmp) > 0) tmp else rep(NA, 2))
        }))
        rownames(pks) <- names(new.fits) <- sapply(new.fits, 
            function(i) rownames(i$data)[1])
        rss <- n.data <- rep(NA, length(new.fits))
        for (f in 1:length(new.fits)) if (!is.null(new.fits[[f]])) {
            rss[f] <- new.fits[[f]]$out.info["rss"]
            n.data[f] <- nrow(new.fits[[f]]$data)
        }
        if (is.na(n.skip)) 
            n.skip <- ceiling(length(rss)/6)
        bad.pks <- rep(FALSE, length(rss))
        if (n.skip > 0) 
            bad.pks <- is.na(rss) | is.na(n.data) | n.data < 
                median(n.data, na.rm = T)/2 | rank(rss/n.data/pks[, 
                2]^2, na.last = T) > length(rss) - n.skip
        out.score <- log(sum(rss[!bad.pks]/n.data[!bad.pks], 
            na.rm = T)) - log(sum(!bad.pks))
        ding <- ""
        if (out.score <= best.score.so.far) {
            best.score.so.far <<- out.score
            best.kernel <<- kernel
            best.params <<- par
            changed <<- TRUE
            ding <- "*"
        }
        cat("| N=", sum(!bad.pks), "WORST=", max(rss[!bad.pks]), 
            "|", out.score, ding, "\n")
        if (plot) {
            obj <- list(par = par, peaks = pks, new.fits = new.fits, 
                is.bad = bad.pks, kernel = kernel, start = list(kernel = in.kernel, 
                  par = start.pars), score = out.score, args = in.args)
            plot.fit.peak.profile(obj, n.peak.plot = 7)
        }
        iter <<- iter + 1
        if (!return.all) 
            return(out.score)
        else return(invisible(list(score = out.score, kernel = kernel, 
            new.fits = new.fits, is.bad = bad.pks)))
    }
    input.start.pars <- start.pars
    if (!"shape" %in% names(start.pars)) 
        start.pars["shape"] <- 7
    if (!"scale" %in% names(start.pars)) 
        start.pars["scale"] <- 50
    if (!"bs.size" %in% names(start.pars)) 
        start.pars["bs.size"] <- 20
    if (!"h.cutoff" %in% names(start.pars)) 
        start.pars["h.cutoff"] <- 15
    if (!"offset" %in% names(start.pars)) 
        start.pars["offset"] <- 0
    parscale <- start.pars * 0 + 1
    if ("shape" %in% names(parscale)) 
        parscale["shape"] <- 1
    if ("scale" %in% names(parscale)) 
        parscale["scale"] <- 5
    if ("bs.size" %in% names(parscale)) 
        parscale["bs.size"] <- 15
    if ("offset" %in% names(parscale)) 
        parscale["offset"] <- min(10, start.pars["offset"]/2)
    if ("h.cutoff" %in% names(parscale)) 
        parscale["h.cutoff"] <- min(50, start.pars["h.cutoff"]/2)
    if (rnd) 
        start.pars <- start.pars + (runif(length(start.pars)) - 
            0.5) * 2 * parscale
    cat(name, "Starting parameters:\n")
    print(start.pars)
    cat(name, "Allowing only these parameters to be optimized:", 
        to.be.fit, "\n")
    if (is.null(in.kernel) || is.na(in.kernel)) 
        in.kernel <- get.profile(start.pars[to.be.fit])
    if (is.null(fits)) 
        fits <- deconv.entire.genome(data = orig.data, kernel = in.kernel, 
            plot = F, verbose = F, no.multicore = no.multicore, 
            ...)
    pks <- get.biggest.peaks(fits$fits.fin, n.peaks = n.peaks)
    print(pks)
    max.steps.2 <- 50
    if (!is.null(list(...)$max.steps)) 
        max.steps.2 <- round(list(...)$max.steps/2)
    quant.cutoff <- "q0.98"
    q.cutoff <- 0
    if (!is.null(list(...)$quant.cutoff)) {
        quant.cutoff <- list(...)$quant.cutoff
        q.cutoff <- as.numeric(gsub("q", "", quant.cutoff))
        q.cutoff <- quantile(data[, 2], probs = q.cutoff, na.rm = T)
    }
    cat("Given cutoff =", q.cutoff, ", we have", length(unique(data[data[, 
        2] >= q.cutoff, 1])), "possible peaks. (N.peaks =", n.peaks, 
        ")\n")
    datas <- list()
    for (i in 1:nrow(pks)) datas[[i]] <- data[rownames(data) == 
        rownames(pks)[i] & data[, 1] >= pks[i, 1] - mini.window * 
        1.1 & data[, 1] <= pks[i, 1] + mini.window * 1.1, , drop = F]
    cat("Using optimization method:", method, "\n")
    if (method != "GA" && method != "SANN") {
        out.params <- optim(start.pars[to.be.fit], get.fit.score, 
            pks = pks, plot = plot, method = method, control = list(reltol = 1e-05, 
                maxit = 200, parscale = parscale[to.be.fit]))
    }
    else if (method == "SANN") {
        out.params <- optim(start.pars[to.be.fit], get.fit.score, 
            pks = pks, plot = plot, method = "SANN", control = list(reltol = 1e-05, 
                temp = 0.01, maxit = 20000, parscale = parscale[to.be.fit]))
    }
    final.pars <- out.params$par
    iter <- 1001
    if (!is.na(in.kernel)) {
        cat(name, "INPUT KERNEL: ")
        get.fit.score(start.pars, pks, kernel = in.kernel, plot = plot)
    }
    cat(name, "START PARAMS: ")
    start.output <- get.fit.score(start.pars, pks, plot = plot, 
        return.all = T)
    cat(name, "FINAL PARAMS: ")
    final.output <- get.fit.score(final.pars, pks, plot = plot, 
        return.all = T)
    final.output$par <- final.pars
    final.output$fits <- fits
    final.output$peaks <- pks
    start.output$par <- start.pars
    final.output$start <- start.output
    final.output$args <- in.args
    if (plot) {
        try(plot.fit.peak.profile(final.output, n.peak.plot = length(final.output$new.fits)))
        plot.fit.peak.profile(final.output, n.peak.plot = 0)
    }
    attr(final.output, "class") <- "fit.peak.profile"
    invisible(final.output)
}
generate.fake.data <-
function (posns = seq(1, 6001, by = 20), n.pts = 5, noise = 0.1, 
    min.pk = 0.1, plot = F, verbose = F, in.pts = NULL, kernel = NULL, 
    tile.size = 100, reps = 1, noise.func = function(data) noise * 
        (1 + sqrt(data)), ...) 
{
    if (is.null(kernel) || is.na(kernel[1])) 
        kernel <- generate.binding.profile(tile.size = 100, ...)
    x <- y <- numeric()
    if (is.null(in.pts) && !is.na(n.pts) && n.pts > 0) {
        x <- sample(min(posns):max(posns), n.pts)
        y <- runif(n.pts, min = min.pk)
        y <- y/max(y)
        if (verbose) 
            print(rbind(x, y))
    }
    else if (!is.null(in.pts)) {
        x <- in.pts[, 1]
        y <- in.pts[, 2]
    }
    data <- rep(0, length(posns))
    k.range <- min(kernel[, 1]):max(kernel[, 1])
    for (i in 1:length(x)) {
        k.inds <- which((x[i] + kernel[, 1]) %in% posns)
        p.inds <- which(posns %in% (x[i] + kernel[k.inds, 1]))
        data[p.inds] <- data[p.inds] + kernel[k.inds, 2] * y[i]
    }
    out.data <- NULL
    for (i in 1:reps) {
        tmp <- data + rnorm(length(data), mean = 0, sd = noise.func(data))
        out.data <- rbind(out.data, cbind(posns, tmp))
    }
    return(list(input = cbind(x, y), data = out.data, kernel = kernel))
}
generate.binding.profile <-
function (fragment.distrib = function(x, ...) dgamma(x, shape = 6, 
    scale = 50), bs.size = 1, tile.size = 50, min.frag.size = 0, 
    positions = seq(0, 1001, by = 50), intensity.scaling = function(x, 
        ...) x, hybridization.prob = function(x, ...) as.integer(x > 
        10), interp = T, plot = F, verbose = F, no.multicore = T, 
    ...) 
{
    if (!exists("kernel.method")) 
        kernel.method <- "mine"
    if (!0 %in% positions) 
        positions <- c(0, positions)
    positions <- sort(positions)
    max.dist <- max(positions)
    fragment.sizes <- fragment.distrib(positions, ...)
    if (min.frag.size > 1) 
        fragment.sizes[positions < min.frag.size] <- 0
    fragment.sizes[positions < bs.size] <- 0
    fragment.sizes[fragment.sizes < 1e-05] <- 0
    tile.start.end <- c(-ceiling(tile.size/2) + 1, ceiling(tile.size/2))
    bs.half.size <- ceiling(bs.size/2)
    tmp.pos <- 1:length(positions)
    out.distrib <- positions * 0
    apply.func <- lapply
    is.parallel <- !no.multicore
    if (is.parallel) 
        is.parallel <- require(multicore, quietly = T, warn.conflicts = F)
    if (is.parallel) 
        is.parallel <- !multicore:::isChild()
    if (is.parallel) 
        is.parallel <- multicore:::detectCores(all.tests = TRUE) > 
            1
    if (is.parallel) {
        require(doMC)
        registerDoMC()
        apply.func <- function(list, FUN, ...) foreach(l = list) %dopar% 
            {
                FUN(l, ...)
            }
    }
    if (is.parallel && verbose) 
        cat("Parallelizing generation of binding profile over", 
            multicore:::detectCores(all.tests = TRUE), "processor cores.\n")
    if (is.parallel) 
        warning("WARNING: If you are running on a Windows system, the 'multicore' option will not work.\nPlease re-start with parameter 'no.multicore=TRUE'.\n")
    tmp <- apply.func(which(fragment.sizes > 0), function(ind) {
        l <- positions[ind]
        if (verbose) 
            cat(l, "of", max(positions[which(fragment.sizes > 
                0)]), "\n")
        possible.fragment.locs <- seq.int(-l + bs.half.size + 
            1, -bs.half.size - 1)
        fragment.covers <- lapply(possible.fragment.locs, function(i) seq.int(i, 
            i + l))
        f.d <- fragment.sizes[ind]
        factor <- max(1, l - 2 * bs.size)
        intens <- intensity.scaling(l, ...)
        for (tile.ind in tmp.pos) {
            tile.pos <- positions[tile.ind]
            tile.s.e <- tile.start.end + tile.pos
            tmp1 <- tile.s.e[1]
            if (tmp1 > l) 
                next
            tmp2 <- tile.s.e[2]
            tile.overlaps <- sapply(fragment.covers, function(i) sum(i >= 
                tmp1 & i <= tmp2))
            tile.overlaps <- tile.overlaps[tile.overlaps > 0]
            sum.overlaps <- sum(hybridization.prob(tile.overlaps, 
                ...))
            out.distrib[tile.ind] <- out.distrib[tile.ind] + 
                sum.overlaps * intens * f.d/factor
        }
        out.distrib
    })
    for (i in tmp) out.distrib <- out.distrib + i
    out.distrib <- out.distrib/max(out.distrib, na.rm = T)
    out.distrib[is.na(out.distrib)] <- 0
    if (all(is.na(out.distrib) | is.infinite(out.distrib) | all(out.distrib == 
        0))) 
        out.distrib[1:length(out.distrib)] <- 1
    if (interp && !all(min(positions):max(positions) %in% seq(positions))) {
        ss <- predict(smooth.spline(positions, out.distrib, keep.data = F), 
            min(positions):max(positions))
        out.distrib <- ss$y
        positions <- ss$x
    }
    out.distrib[out.distrib > 1] <- 1
    out.distrib[out.distrib < 0] <- 0
    out.distrib <- c(rev(out.distrib), out.distrib[-1])
    posns <- unique(c(-rev(positions), positions))
    if (plot) 
        plot(posns, out.distrib, typ = "l")
    invisible(cbind(posns, out.distrib))
}
get.chip.boot.probs <-
function (obj, boot.results = "prob", window = NULL, hi.res = NA, 
    smooth = T) 
{
    is.boot <- is.null(obj$args) && !is.null(obj[[1]]$args) && 
        obj[[1]]$args$n.boot > 1
    if (!is.boot && !is.null(obj$args)) {
        tmp <- list(obj = obj)
        obj = tmp
    }
    coeffs <- NULL
    for (ii in 1:length(obj)) coeffs <- rbind(coeffs, obj[[ii]]$coeffs)
    if (is.null(window)) 
        window <- obj[[1]]$window
    w.expand <- window + c(-1000, 1000)
    coeffs <- coeffs[!is.na(coeffs[, 1]) & coeffs[, 1] >= w.expand[1] & 
        coeffs[, 1] <= w.expand[2], , drop = F]
    if (nrow(coeffs) <= 0) {
        coeffs <- cbind(w.expand[1] - 100, 0)
    }
    rng <- range(coeffs[, 1])
    tmp.hi.res <- hi.res
    if (is.null(tmp.hi.res) || is.na(tmp.hi.res)) 
        tmp.hi.res <- obj[[1]]$args$fit.res
    for (boot.res in boot.results) {
        if (smooth) {
            scal <- 1
            if (boot.res == "prob") {
                dens <- density(coeffs[, 1], bw = tmp.hi.res, 
                  n = diff(range(coeffs[, 1])))
                scal <- sum(abs(coeffs[, 1] - dens$x[which.max(dens$y)]) <= 
                  tmp.hi.res/2)/length(obj)
            }
            else if (boot.res == "scaled.prob" || boot.res == 
                "prob.scaled") {
                wts <- coeffs[, 2]
                wts[wts < 0] <- 0
                if (sum(wts) == 0) 
                  wts[wts == 0] <- 1e-05
                dens <- density(coeffs[, 1], bw = tmp.hi.res, 
                  weights = wts/sum(wts), n = diff(range(coeffs[, 
                    1])))
                scal <- median(coeffs[which(abs(coeffs[, 1] - 
                  dens$x[which.max(dens$y)]) <= tmp.hi.res/2), 
                  2])
            }
            else if (boot.res == "scale") {
                mns <- tapply(coeffs[, 2], round(coeffs[, 1]/tmp.hi.res) * 
                  tmp.hi.res, mean, na.rm = T)
                mns[mns < 0] <- 0
                dens <- density(as.numeric(names(mns)), weights = mns/sum(mns), 
                  bw = tmp.hi.res, n = diff(range(coeffs[, 1])))
                scal <- max(coeffs[which(abs(coeffs[, 1] - dens$x[which.max(dens$y)]) <= 
                  tmp.hi.res/2), 2])
            }
            hst.mid <- dens$x
            hst.cnt <- dens$y/max(dens$y) * scal
            ccoeffs <- cbind(mids = hst.mid, counts = hst.cnt, 
                hst.cnt)
            colnames(ccoeffs)[3] <- boot.res
            return(ccoeffs)
        }
        hst <- hist(coeffs[, 1], breaks = seq(floor(rng[1]/tmp.hi.res) * 
            tmp.hi.res, ceiling(rng[2]/tmp.hi.res) * tmp.hi.res, 
            by = tmp.hi.res), plot = F)
        hst.mid <- hst$mids
        hst.cnt <- hst$counts
        hst.prob <- hst.cnt/length(obj)
        ccoeffs <- cbind(mids = hst.mid, counts = hst.cnt)
        tmp.hi.res.over.two <- tmp.hi.res/2
        if (boot.res == "prob") {
            hst.prob[hst.prob > 1] <- 1
            ccoeffs <- cbind(ccoeffs, prob = hst.prob)
        }
        else if (boot.res == "scaled.prob" || boot.res == "prob.scaled" || 
            boot.res == "scale") {
            ttmp.ccoeffs <- rbind(c(0, 0), ccoeffs, c(9e+09, 
                0))
            scal <- rep(0, nrow(ttmp.ccoeffs) - 2)
            if (boot.res == "scaled.prob" || boot.res == "prob.scaled") {
                for (i in 2:(nrow(ttmp.ccoeffs) - 1)) {
                  tt <- ttmp.ccoeffs[i, 1]
                  scal[i - 1] <- sum(coeffs[abs(coeffs[, 1] - 
                    tt) <= tmp.hi.res.over.two, 2], na.rm = T)
                }
                ccoeffs <- cbind(ccoeffs, scaled.prob = scal/length(obj))
            }
            else if (boot.res == "scale") {
                for (i in 2:(nrow(ttmp.ccoeffs) - 1)) {
                  tt <- ttmp.ccoeffs[i, 1]
                  scal[i - 1] <- mean(coeffs[abs(coeffs[, 1] - 
                    tt) <= tmp.hi.res.over.two, 2], na.rm = T)
                }
                scal[is.na(scal)] <- 0
                ccoeffs <- cbind(ccoeffs, scale = scal)
            }
        }
    }
    ccoeffs
}
plot.fit.peak.profile <-
function (x, n.peak.plot = 7, plot.spline = F, ...) 
{
    obj <- x
    pks <- obj$peaks
    new.fits <- obj$new.fits
    if (n.peak.plot < -1) 
        cat("HERE\n")
    else if (n.peak.plot <= 1) 
        par(mfrow = c(1, 1))
    else if (n.peak.plot < 3) 
        par(mfrow = c(2, 2))
    else par(mfrow = c(3, 3))
    rang <- 1:max(x$args$positions)
    if (n.peak.plot > 0) {
        n.plotted <- 0
        for (i in 1:length(new.fits)) {
            if (is.null(new.fits[[i]])) 
                next
            try(plot(new.fits[[i]], main = sprintf("%d: %s RSS = %.2f %s", 
                i, rownames(pks)[i], new.fits[[i]]$out.info["rss"], 
                ifelse(obj$is.bad[i], "*", "")), pch = 19, cex = 0.5))
            n.plotted <- n.plotted + 1
            if (n.plotted >= n.peak.plot) 
                break
        }
        par <- abs(obj$par)
        par2 <- abs(obj$start$par)
        plot(rang, dgamma(rang, shape = par["shape"], scale = par["scale"]), 
            typ = "l", xlab = "Fragment length", ylab = "Fragment length distrib.", 
            main = sprintf("FINAL: mean = %.3f\nSTART: mean = %.3f", 
                par["shape"] * par["scale"], par2["shape"] * 
                  par2["scale"]))
        lines(rang, dgamma(rang, shape = par2["shape"], scale = par2["scale"]), 
            col = "green", lty = 3)
    }
    wheres <- unlist(sapply(new.fits, function(i) rownames(i$data)[1]))
    tmp <- list()
    for (w in unique(wheres)) {
        tmp[[w]] <- list(list(data = NULL, coeffs = NULL))
        ind <- 1
        for (ww in which(wheres == w)) {
            if (obj$is.bad[ww]) 
                next
            tmp[[w]][[1]]$data <- rbind(tmp[[w]][[1]]$data, new.fits[[ww]]$data)
            tmp[[w]][[1]]$coeffs <- rbind(tmp[[w]][[1]]$coeffs, 
                new.fits[[ww]]$coeffs)
            ind <- ind + 1
        }
    }
    tmp <- get.peak.profile(tmp, n.peaks = nrow(pks), out.to = x$args$mini.window, 
        filter = F)
    plot(tmp, pch = 20, cex = 0.5, main = c(paste(sprintf("%.3f", 
        abs(obj$par)), collapse = " "), sprintf("%.4f", obj$score)), 
        xlab = "Offset (bp)", ylab = "Scaled intensity")
    lines(obj$kernel, col = "red", lwd = 3)
    if (plot.spline) {
        lines(obj$start$kernel, col = "green", lty = 3)
        ss <- predict(smooth.spline(tmp, keep.data = F), min(tmp[, 
            1]):max(tmp[, 1]))
        lines(ss$x, ss$y, col = "blue", lty = 3)
    }
}
get.peak.profile <-
function (obj, n.peaks = 10, out.to = 1000, filter = T, plot.do = F) 
{
    if (!is.null(obj$coeffs)) 
        obj <- list(obj = list(obj))
    best.coeffs <- get.biggest.peaks(obj, n.peaks)
    out.data <- NULL
    done.coords <- numeric()
    done.where <- character()
    for (i in 1:nrow(best.coeffs)) {
        where <- rownames(best.coeffs)[i]
        coord <- best.coeffs[i, 1]
        height <- best.coeffs[i, 2]
        if (length(done.coords) > 0 && sum(done.where == where) > 
            0 && any(abs(done.coords[done.where == where] - coord) < 
            10)) 
            next
        done.coords <- c(done.coords, coord)
        done.where <- c(done.where, where)
        data <- obj[[where]]$data
        if (is.null(data)) 
            data <- obj[[where]][[1]]$data
        tmp <- data[abs(data[, 1] - coord) <= out.to, ]
        tmp[, 1] <- tmp[, 1] - coord
        tmp[, 2] <- tmp[, 2]/height
        out.data <- rbind(out.data, tmp)
    }
    if (filter) 
        out.data <- out.data[abs(out.data[, 1]) < out.to/2 | 
            out.data[, 2] < 0.4, ]
    x <- out.data[, 1]
    y <- out.data[, 2]
    if (filter) {
        bad.inds <- which(abs(x) > 600 & y > median(y[abs(x) > 
            600], na.rm = T) + 3 * mad(y[abs(x) > 600], na.rm = T))
        while (length(bad.inds) > 0) {
            cat("Removing", length(bad.inds), "points.\n")
            x[bad.inds] <- NA
            y[bad.inds] <- NA
            bad.inds <- which(abs(x) > 600 & y > median(y[abs(x) > 
                600], na.rm = T) + 3 * mad(y[abs(x) > 600], na.rm = T))
        }
        x <- x[!is.na(x)]
        y <- y[!is.na(y)]
        y <- y - median(y[x > 650], na.rm = T)
        if (plot.do) 
            plot(x, y, pch = 20, cex = 0.7)
    }
    cbind(x, y)
}
get.biggest.peaks <-
function (obj, n.peaks = 20) 
{
    if (!is.null(obj$coeffs)) 
        obj <- list(obj = list(obj))
    tmp <- NULL
    wheres <- character()
    for (where in names(obj)) {
        coeffs <- obj[[where]][[1]]$coeffs
        if (is.null(coeffs)) 
            next
        tmp <- rbind(tmp, coeffs)
        wheres <- c(wheres, rep(where, nrow(coeffs)))
    }
    rownames(tmp) <- wheres
    tmp <- tmp[order(tmp[, 2], decreasing = T), , drop = F]
    if (nrow(tmp) <= n.peaks) 
        return(tmp)
    tmp <- tmp[tmp[, 2] > 0.1, , drop = F]
    if (nrow(tmp) <= n.peaks) 
        return(tmp)
    best.coeffs <- tmp[1, , drop = F]
    np <- 2
    while (sum(!duplicated(round(best.coeffs[, 1]/1000) * 1000)) < 
        n.peaks && nrow(best.coeffs) < n.peaks * 100) {
        np <- np + 1
        best.coeffs <- tmp[1:np, , drop = F]
    }
    best.coeffs <- best.coeffs[!duplicated(as.integer(best.coeffs[, 
        1]/10) * 10), , drop = F]
    best.coeffs
}
make.predictor.matrix <-
function (posns, kernel, fit.res = 1, verbose = F, sparse = T, 
    good.posns.hires) 
{
    posn.range <- diff(range(posns))
    sp.string <- if (sparse) 
        "SPARSE"
    else ""
    if (verbose) 
        cat("TRYING TO ALLOCATE A", sp.string, "PREDICTOR MATRIX OF SIZE", 
            length(posns), "x", length(good.posns.hires), "\n")
    if (length(posns) < 10000) {
        xx <- matrix(0, nrow = length(posns), ncol = length(good.posns.hires))
        for (i in 1:ncol(xx)) xx[, i] <- convolve.func(good.posns.hires[i], 
            kernel, posns)[posns]
        xx[is.na(xx)] <- 0
    }
    else {
        step.by <- 250
        xxx <- matrix(0, nrow = length(posns), ncol = step.by)
        if (require(Matrix, quietly = T, warn.conflicts = F)) 
            xx <- Matrix(0, nrow = length(posns), ncol = length(good.posns.hires), 
                sparse = T)
        else xx <- matrix(0, nrow = length(posns), ncol = length(good.posns.hires))
        tmp.rows <- 1:length(posns)
        ii <- 1
        for (i in 1:ncol(xx)) {
            xxx[, (i%%step.by + 1)] <- convolve.func(good.posns.hires[i], 
                kernel, posns)[posns]
            if (i%%step.by == 0 || i == ncol(xx)) {
                xxx[is.na(xxx)] <- 0
                if (require(Matrix, quietly = T, warn.conflicts = F)) 
                  xx[tmp.rows, ii:(ii + (i - 1)%%step.by)] <- Matrix(xxx[tmp.rows, 
                    1:(1 + (i - 1)%%step.by)], sparse = T)
                else xx[tmp.rows, ii:(ii + (i - 1)%%step.by)] <- xxx[tmp.rows, 
                  1:(1 + (i - 1)%%step.by)]
                ii <- ii + step.by
                if (i >= ncol(xx)) 
                  break
            }
        }
        rm(xxx)
        gc()
    }
    colnames(xx) <- as.character(1:ncol(xx))
    if (sparse && class(xx) == "matrix" && require(Matrix, quietly = T, 
        warn.conflicts = F)) 
        xx <- Matrix(xx, sparse = T)
    if (verbose) 
        cat("DONE MAKING", sp.string, "PREDICTOR MATRIX:", dim(xx), 
            object.size(xx), "\n")
    xx
}
plot.chip.deconv <-
function (x, boot.results = c("scaled.prob", "prob", "scale", 
    "conf=95", "NONE")[1], where = NA, center = NA, window = NULL, 
    verbose = F, plot.genes = F, org = NA, hi.res = NA, quants = c(0.95, 
        0.5, 0.05), smooth = T, ylim = NULL, ...) 
{
    if (!exists("scale.boot.results")) 
        scale.boot.results <- NA
    in.args <- c(mget(names(formals()), envir = as.environment(-1)), 
        sapply(as.list(substitute({
            ...
        })[-1]), deparse))
    obj <- x
    is.boot <- is.null(obj$args) && !is.null(obj[[1]]$args) && 
        obj[[1]]$args$n.boot > 1
    if (!is.boot && !is.null(obj$args)) {
        tmp <- list(obj = obj)
        obj = tmp
    }
    if (is.null(window)) 
        window <- obj[[1]]$window
    if (length(window) == 1) {
        if (is.na(center)) 
            center <- mean(obj[[1]]$window)
        window <- round(c(-window/2, window/2)) + center
    }
    if (is.null(hi.res) || is.na(hi.res)) 
        hi.res <- obj[[1]]$args$fit.res
    fit.bg <- NA
    w.expand <- window + c(-1000, 1000)
    dat <- obj[[1]]$data
    dat <- dat[dat[, 1] >= w.expand[1] & dat[, 1] <= w.expand[2], 
        , drop = F]
    y.range <- range(c(0, dat[, 2]))
    yr <- y.range
    if (plot.genes) 
        y.range[1] <- y.range[1] - diff(y.range)/10
    x.range <- round(range(dat[, 1]))
    coeffs <- NULL
    for (ii in 1:length(obj)) coeffs <- rbind(coeffs, obj[[ii]]$coeffs)
    spline.coeffs <- coeffs[grep("spline", rownames(coeffs), 
        fixed = T), , drop = F]
    if (nrow(coeffs) > 1) 
        coeffs <- coeffs[!is.na(coeffs[, 1]) & coeffs[, 1] >= 
            w.expand[1] & coeffs[, 1] <= w.expand[2], , drop = F]
    if (nrow(coeffs) > 1) 
        coeffs <- coeffs[order(coeffs[, 1]), , drop = F]
    coeffs[, 1] <- coeffs[, 1] - x.range[1]
    coe <- obj[[1]]$coeffs
    coe <- coe[coe[, 1] >= w.expand[1] & coe[, 1] <= w.expand[2], 
        , drop = F]
    out.scale <- NULL
    posns <- dat[, 1]
    preds <- t(rep(0, length(posns)))
    if (!is.na(hi.res) && hi.res > 0 && !is.na(quants) && quants != 
        "NONE" && nrow(coeffs) > 0) {
        posns <- seq(1, diff(range(dat[, 1])), by = hi.res)
        xx <- make.predictor.matrix(posns, obj[[1]]$kernel, fit.res = hi.res, 
            good.posns.hires = unique(round(coeffs[, 1])), sparse = F, 
            verbose = verbose)
        colnames(xx) <- as.character(unique(round(coeffs[, 1])) + 
            x.range[1])
        preds <- matrix(0, nrow = length(obj), ncol = nrow(xx))
        for (i in 1:length(obj)) {
            ccoeffs <- obj[[i]]$coeffs
            if (nrow(ccoeffs) <= 0) 
                next
            ccoeffs <- ccoeffs[ccoeffs[, 1] >= w.expand[1] & 
                ccoeffs[, 1] <= w.expand[2], , drop = F]
            if (nrow(ccoeffs) <= 0) 
                next
            ccoeffs <- ccoeffs[order(ccoeffs[, 1]), , drop = F]
            tmp.pos <- unique(round(ccoeffs[, 1]))
            whch <- !as.character(tmp.pos) %in% colnames(xx)
            if (any(whch)) {
                tmp.pos[whch] <- tmp.pos[whch] + 1
                whch <- !as.character(tmp.pos) %in% colnames(xx)
                if (any(whch)) {
                  tmp.pos <- unique(round(ccoeffs[, 1]))
                  tmp.pos[whch] <- tmp.pos[whch] - 2
                  whch <- !as.character(tmp.pos) %in% colnames(xx)
                }
                if (any(whch)) {
                  tmp.pos <- unique(round(ccoeffs[, 1]))
                  tmp.pos <- tmp.pos[!whch]
                }
            }
            whch <- !as.character(tmp.pos) %in% colnames(xx)
            if (any(whch) || !any(!whch)) 
                next
            tmp.coe <- cbind(tmp.pos, ccoeffs[!whch, 2])
            xxx <- xx[, as.character(round(tmp.coe[, 1])), drop = F]
            preds[i, ] <- (xxx %*% tmp.coe[, 2])[, 1]
            if (i > 1 && (obj[[1]]$args$boot.sample.opt %in% 
                c("resample", "residual"))) {
                preds <- preds[1, , drop = F]
                break
            }
        }
        if (nrow(preds) > 1) {
            preds <- t(sapply(quants, function(q) apply(preds, 
                2, quantile, probs = q, na.rm = T)))
        }
    }
    else {
        if (!is.null(obj[[1]]$fit)) 
            preds <- t(obj[[1]]$fit[, 2])
    }
    boot.res.names <- c(prob = "Probability", scaled.prob = "Scaled prob.", 
        scale = "Scale", `conf=95` = "95% Confidence")
    y.range <- range(y.range, preds, coe[, 2], 0)
    pch <- in.args$pch
    if (is.null(pch)) 
        pch <- 20
    cex <- in.args$cex
    if (is.null(cex)) 
        cex <- 0.5
    if (is.null(ylim)) 
        ylim <- y.range
    old.par <- par(pch = pch, cex = cex, ...)
    on.exit(par(old.par))
    plot(dat, xlim = window, ylim = ylim, xlab = "Genome coord.", 
        ylab = "Chip intensity", ...)
    apply(coe, 1, function(i) lines(rep(i[1], 2), c(0, i[2]), 
        col = "darkgreen"))
    points(coe, col = "darkgreen", pch = 20)
    if (!is.null(preds)) 
        for (i in 1:nrow(preds)) lines(posns + x.range[1], preds[i, 
            ], col = "red", lty = i + 1)
    if (is.boot && obj[[1]]$args$boot.sample.opt %in% c("case", 
        "wild", "position", "replicate")) {
        tmp.hi.res <- hi.res
        if (is.null(tmp.hi.res) || is.na(tmp.hi.res)) 
            tmp.hi.res <- obj[[1]]$args$fit.res
        if (is.null(tmp.hi.res) || is.na(tmp.hi.res)) 
            tmp.hi.res <- 10
        for (boot.res in boot.results) {
            if (substr(boot.res, 1, 4) == "conf") {
                coeffs <- NULL
                for (ii in 1:length(obj)) coeffs <- rbind(coeffs, 
                  obj[[ii]]$coeffs)
                coeffs <- coeffs[!is.na(coeffs[, 1]) & coeffs[, 
                  1] >= w.expand[1] & coeffs[, 1] <= w.expand[2], 
                  , drop = F]
                rng <- range(coeffs[, 1])
                tmp1 <- seq(floor(rng[1]), ceiling(rng[2]), by = 1)
                tmp2 <- rep(0, length(tmp1))
                tmp.lst <- list()
                for (ii in 1:length(obj)) {
                  c <- obj[[ii]]$coeffs
                  ttmp <- tmp2
                  ttmp[tmp1 %in% round(c[, 1])] <- c[, 2]
                  tmp.lst[[ii]] <- cbind(tmp1, ttmp)
                }
                coeffs <- do.call(rbind, tmp.lst)
                coeffs <- coeffs[coeffs[, 1] >= w.expand[1] & 
                  coeffs[, 1] <= w.expand[2], , drop = F]
                conf.int <- 0.95 * length(obj)
                if (substr(boot.res, 5, 5) == "=") 
                  conf.int <- as.integer(strsplit(boot.res, "=")[[1]][2])/100 * 
                    length(obj)
                ccoeffs <- coeffs[coeffs[, 2] > 0, , drop = F]
                tile.dist <- max(obj[[1]]$out.info["tile.distance"], 
                  100)
                for (i in 1:nrow(coe)) {
                  coord <- coe[i, 1]
                  dists <- abs(ccoeffs[, 1] - coord)
                  for (dist in seq(1, tile.dist, by = tmp.hi.res)) {
                    frac <- sum(dists <= dist)
                    if (frac >= conf.int) 
                      break
                  }
                  if (frac >= conf.int) {
                    hts <- ccoeffs[dists < dist, 2]
                    hts.low.hi <- quantile(hts, c(1 - conf.int/length(obj), 
                      conf.int/length(obj)), na.rm = T)
                    lines(rep(coord - dist, 2), c(0, hts.low.hi[2]), 
                      col = "green")
                    lines(rep(coord + dist, 2), c(0, hts.low.hi[2]), 
                      col = "green")
                    if (dist < tile.dist) 
                      dist <- tile.dist
                    lines(rep(coord, 2), rep(hts.low.hi[2], 2), 
                      col = "green")
                    lines(c(coord - dist, coord + dist), rep(hts.low.hi[1], 
                      2), col = "green")
                    lines(c(coord - dist, coord + dist), rep(hts.low.hi[2], 
                      2), col = "green")
                  }
                }
            }
            else {
                tmp <- get.chip.boot.probs(obj, boot.results = boot.res, 
                  window = w.expand, hi.res = tmp.hi.res, smooth = smooth)
                ccoeffs <- tmp[, c(1, 3), drop = F]
                coef.rng <- range(ccoeffs[, 2])
                if (is.na(scale.boot.results)) 
                  scale.boot.results <- y.range[2]
                out.scale <- scale.boot.results/max(ccoeffs[, 
                  2])
                ccoeffs[, 2] <- ccoeffs[, 2] * out.scale
                lines(rbind(cbind(c(-999, min(ccoeffs[, 1]) - 
                  tmp.hi.res), c(0, 0)), ccoeffs, cbind(c(max(ccoeffs[, 
                  1]) + tmp.hi.res, 9e+09), c(0, 0))), col = "green")
                lines(ccoeffs, col = "green")
                tmp1 <- seq(0, max(coef.rng), length = 4)
                tmp2 <- seq(0, scale.boot.results, length = length(tmp1))
                axis(4, tmp2, labels = sprintf("%.1f", tmp1), 
                  ...)
                m.cex <- 1
                mtext(boot.res.names[boot.res], 4, line = 0, 
                  padj = 1.2, outer = F, cex = m.cex)
            }
        }
    }
    if (is.na(where)) 
        where <- obj[[1]]$args$where
    if (plot.genes) 
        plot.genes.in.window(mean(w.expand), where, diff(range(w.expand)), 
            new.plot = F, yoff = y.range[1] + 0.4, org = org, 
            yscale = abs(yr[1] - y.range[1])/2, ...)
    if (!is.null(out.scale)) 
        return(invisible(out.scale))
}
plot.chip.deconv.entire.genome <-
function (x, where = NA, center = NA, window = NULL, ...) 
{
    fit <- x
    wheres <- where
    if (is.na(wheres)) 
        wheres <- names(fit$fits.fin)
    for (where in wheres) {
        par(mfrow = c(2, 1), mar = c(2, 2, 2, 2), mgp = c(3, 
            1, 0) * 0.3, lwd = 2)
        f <- fit$fits.fin[[where]]
        wind.size <- window
        if (is.null(window)) 
            wind.size <- f[[1]]$window
        if (is.na(center)) {
            if (!is.null(f$data)) 
                rng <- range(f$data[, 1])
            else rng <- range(f[[1]]$data[, 1])
            for (i in seq(rng[1] + wind.size/2, rng[2], by = wind.size * 
                2/3)) {
                wind <- round(i + c(-wind.size, wind.size)/2)
                cat(where, i, wind, "\n")
                plot(f, where = where, window = wind, main = paste(where, 
                  ": ", wind[1], "-", wind[2], sep = ""), ...)
            }
        }
        else {
            for (cen in center) plot(f, where = where, window = wind.size, 
                center = cen, ...)
        }
    }
}
plot.genes.in.window <-
function (center, where, window, new.plot = T, yoff = 0, yscale = 1, 
    org = NA, gene.names = T, main = "", gene.coords = NULL, 
    ...) 
{
    if (is.null(gene.coords)) {
        if (exists("gene.coords", envir = .GlobalEnv)) 
            coords <- get("gene.coords", envir = .GlobalEnv)
    }
    window <- round(c(center - window/2, center + window/2))
    w.expand <- round(c(window[1] - diff(window)/5, window[2] + 
        diff(window)/5))
    in.window <- function(x, w) {
        x >= w[1] & x <= w[2]
    }
    genes.in <- coords[tolower(as.character(coords[, "where"])) == 
        tolower(where) & (in.window(as.integer(as.vector(coords[, 
        "Start"])), w.expand) | in.window(as.integer(as.vector(coords[, 
        "Stop"])), w.expand)), , drop = F]
    if (nrow(genes.in) <= 0) 
        return()
    genes.x <- apply(cbind(as.integer(as.vector(genes.in[, "Start"])), 
        as.integer(as.vector(genes.in[, "Stop"]))), 1, mean)
    genes.is.fwd <- as.character(genes.in[, "Orientation"]) == 
        "For"
    if (new.plot) {
        plot(0, 0, xlim = window, ylim = c(-1, 1), ann = F, xaxt = "n", 
            yaxt = "n", bty = "n")
        axis(side = 1, pos = 0, mgp = c(0, 0.5, 0))
    }
    for (i in 1:nrow(genes.in)) {
        start <- as.integer(as.vector(genes.in[i, "Start"]))
        end <- as.integer(as.vector(genes.in[i, "Stop"]))
        name <- genes.in[i, "Gene_Name"]
        if (genes.in[i, "Orientation"] == "For") {
            rect(start, yoff + yscale * 0.2, end, yoff + yscale * 
                1, col = "yellow", border = "black", lwd = 1)
            if (gene.names) 
                text(genes.x[i], yoff + yscale * 0.6, labels = name, 
                  adj = c(0.5, 0.5), col = "black", ...)
        }
        else if (genes.in[i, "Orientation"] == "Rev") {
            rect(start, yoff + yscale * -1, end, yoff + yscale * 
                -0.2, col = "orange", border = "black", lwd = 1)
            if (gene.names) 
                text(genes.x[i], yoff + yscale * -0.6, labels = name, 
                  adj = c(0.5, 0.5), col = "black", ...)
        }
    }
}
post.proc.coeffs <-
function (coeffs, fit.res, max.coef = NA, factor = 6, mean.do = F) 
{
    coeffs <- coeffs[coeffs[, 2] > 0, , drop = F]
    if (nrow(coeffs) <= 1) 
        return(coeffs)
    if (is.na(max.coef)) 
        max.coef <- 100
    coeffs[coeffs[, 2] > max.coef, 2] <- max.coef
    ind <- 1
    old.row <- coeffs[1, ]
    tmp <- c(old.row, ind)
    for (i in 2:nrow(coeffs)) {
        row <- coeffs[i, ]
        if (row[1] - old.row[1] > fit.res * factor || (row[1] - 
            old.row[1] > fit.res * factor * 2 && (row[2] > old.row[2] * 
            5 || old.row[2] >= row[2] * 5))) 
            ind <- ind + 1
        tmp <- rbind(tmp, c(row, ind))
        old.row <- row
    }
    coeffs.out <- NULL
    for (i in unique(tmp[, 3])) {
        rows <- tmp[tmp[, 3] == i, , drop = F]
        if (!mean.do) 
            out.row <- c(weighted.mean(rows[, 1], rows[, 2]), 
                sum(rows[, 2]))
        else out.row <- c(weighted.mean(rows[, 1], rows[, 2]), 
            mean(rows[, 2]))
        coeffs.out <- rbind(coeffs.out, out.row)
    }
    rownames(coeffs.out) <- NULL
    coeffs.out
}
post.proc.deconv <-
function (obj, ...) 
{
    coeffs <- obj$coeffs
    if (nrow(coeffs) <= 1) {
        warning("Only zero or one coefficient -- no post.proc necessary!")
        return(obj)
    }
    coeffs <- coeffs[coeffs[, 2] > 0, , drop = F]
    if (nrow(coeffs) <= 1) {
        warning("Only zero or one coefficient -- no post.proc necessary!")
        return(obj)
    }
    obj$coeffs <- post.proc.coeffs(coeffs, ...)
    obj
}
print.chip.deconv <-
function (x, ...) 
{
    obj <- x
    is.boot <- is.null(obj$args) && !is.null(obj[[1]]$args) && 
        obj[[1]]$args$n.boot > 1
    if (!is.boot && !is.null(obj$args)) {
        tmp <- list(obj = obj)
        obj = tmp
    }
    cat("MeDiChI deconvolution over window:", obj[[1]]$window, 
        "\n")
    cat("Input arguments:\n")
    args <- c("center", "where", "fit.res", "max.steps", "n.boot", 
        "boot.sample.opt", "quant.cutoff")
    for (i in args) {
        if (!is.null(obj[[1]]$args[[i]])) 
            cat("\t", i, "=", obj[[1]]$args[[i]])
        if (which(args == i)%%2 == 0 || which(args == i) == length(args)) 
            cat("\n")
        else cat("\t")
    }
    if (!is.null(obj[[1]]$out.info)) {
        cat("Info:\n")
        print(obj[[1]]$out.info, ...)
    }
    cat("Coeffs:\n")
    print(coef(obj[[1]], ...), ...)
}
print.chip.deconv.entire.genome <-
function (x, ...) 
{
    fit <- x
    cat("MeDiChI entire genome deconvolution:\n")
    for (chr in names(fit$fits.fin)) {
        cat("Chromosome", chr, ":\n")
        print(fit$fits.fin[[chr]])
    }
}
coef.chip.deconv <-
function (object, ...) 
{
    obj <- object
    is.boot <- is.null(obj$args) && !is.null(obj[[1]]$args) && 
        obj[[1]]$args$n.boot > 1
    if (!is.boot && !is.null(obj$args)) {
        tmp <- list(obj = obj)
        obj = tmp
    }
    if (!is.null(obj[[1]]$coeffs.w.p.values)) 
        return(obj[[1]]$coeffs.w.p.values)
    return(obj[[1]]$coeffs)
}
coef.chip.deconv.entire.genome <-
function (object, ...) 
{
    obj <- object
    out <- NULL
    for (chr in names(obj$fits.fin)) {
        tmp <- coef(obj$fits.fin[[chr]], ...)
        rownames(tmp) <- rep(chr, nrow(tmp))
        out <- rbind(out, tmp)
    }
    out
}
load.chip.data <-
function (data, verbose = T, filter.spurious.probes = F) 
{
    if (is.matrix(data)) {
        fname <- ""
    }
    else {
        fname <- data
        if (is.character(data) && file.exists(data)) {
            cat("Loading data file:", data, "\n")
            data <- read.table(data)
        }
        else if (any(class(data) == "connection") && summary(data)[["can read"]] == 
            "yes") {
            cat("Loading data from connection... file:", summary(data)$description, 
                "\n")
            fname <- summary(data)$description
            data <- read.table(data)
        }
        if (is.character(fname) && length(grep(".gff", fname, 
            fixed = T)) > 0) {
            cat("Data is a GFF file - extracting relevant information.\n")
            gff <- data
            colnames(gff) <- c("seqname", "source", "feature", 
                "start", "end", "score", "strand", "frame")
            data <- as.matrix(cbind(apply(gff[, c("start", "end")], 
                1, mean), gff$score))
            rownames(data) <- as.character(gff$seqname)
            if (min(data[, 2]) < 0) {
                cat("Found negative scores - assuming log_2 intensities; unlogging.\n")
                data[, 2] <- 2^data[, 2]
            }
        }
        if (is.data.frame(data)) {
            if (ncol(data) == 2) 
                data <- as.matrix(data)
            else if (ncol(data == 3)) {
                tmp <- as.matrix(data[, 2:3])
                rownames(tmp) <- as.character(data[, 1])
                data <- tmp
                rm(tmp)
                gc()
            }
            fname <- ""
        }
        if (is.null(rownames(data))) 
            rownames(data) <- rep("XXX", nrow(data))
        if (verbose) {
            cat("Got", nrow(data), "measurements; ranges:", range(data[, 
                1]), "; ", range(data[, 2]), "; chromosomes:\n")
            print(unique(rownames(data)))
        }
        attr(data, "filename") <- fname
    }
    data
}
my.drop <-
function (x) 
{
    q = dim(x)
    if (length(q) == 0) 
        x
    else if (q[1] == 1) 
        x[1, ]
    else if (q[2] == 1) 
        x[, 1]
    else x
}
lars.pos <-
function (x, y, type = c("lasso", "lar", "forward.stagewise"), 
    trace = FALSE, Gram, eps = .Machine$double.eps, max.steps, 
    use.Gram = TRUE, positive = FALSE, ...) 
{
    call <- match.call()
    type <- match.arg(type)
    TYPE <- switch(type, lasso = "LASSO", lar = "LAR", forward.stagewise = "Forward Stagewise")
    if (trace) 
        cat(paste(TYPE, "sequence\n"))
    nm <- dim(x)
    n <- nm[1]
    m <- nm[2]
    im <- inactive <- seq(m)
    one <- rep(1, n)
    vn <- dimnames(x)[[2]]
    meanx <- (one %*% x)[1, ]/n
    x <- scale(x, meanx, FALSE)
    normx <- sqrt((one %*% (x^2))[1, ])
    nosignal <- normx/sqrt(n) < eps
    if (any(nosignal)) {
        ignores <- im[nosignal]
        inactive <- im[-ignores]
        normx[nosignal] <- eps * sqrt(n)
        if (trace) 
            cat("LARS Step 0 :\t", sum(nosignal), "Variables with Variance < eps; dropped for good\n")
    }
    else ignores <- NULL
    names(normx) <- NULL
    x <- scale(x, FALSE, normx)
    if (use.Gram & missing(Gram)) {
        if (m > 500 && n < m) 
            cat("There are more than 500 variables and n<m;\nYou may wish to restart and set use.Gram=FALSE\n")
        if (trace) 
            cat("Computing X'X .....\n")
        Gram <- t(x) %*% x
    }
    mu <- mean(y)
    y <- my.drop(y - mu)
    Cvec <- my.drop(t(y) %*% x)
    ssy <- sum(y^2)
    residuals <- y
    if (missing(max.steps)) 
        max.steps <- 8 * min(m, n - 1)
    beta <- matrix(0, max.steps + 1, m)
    Gamrat <- NULL
    arc.length <- NULL
    R2 <- 1
    RSS <- ssy
    first.in <- integer(m)
    active <- NULL
    actions <- as.list(seq(max.steps))
    drops <- FALSE
    Sign <- NULL
    R <- NULL
    k <- 0
    while ((k < max.steps) & (length(active) < min(m - length(ignores), 
        n - 1))) {
        action <- NULL
        k <- k + 1
        C <- Cvec[inactive]
        if (!positive) 
            Cmax <- max(abs(C))
        else Cmax <- max(C)
        if (!any(drops)) {
            if (!positive) 
                new <- abs(C) >= Cmax - eps
            else new <- C >= Cmax - eps
            C <- C[!new]
            new <- inactive[new]
            for (inew in new) {
                if (use.Gram) {
                  R <- updateR.faster(Gram[inew, inew], R, my.drop(Gram[inew, 
                    active]), Gram = TRUE, eps = eps)
                }
                else {
                  R <- updateR.faster(x[, inew], R, x[, active], 
                    Gram = FALSE, eps = eps)
                }
                if (attr(R, "rank") == length(active)) {
                  nR <- seq(length(active))
                  R <- R[nR, nR, drop = FALSE]
                  attr(R, "rank") <- length(active)
                  ignores <- c(ignores, inew)
                  action <- c(action, -inew)
                  if (trace) 
                    cat("LARS Step", k, ":\t Variable", inew, 
                      "\tcollinear; dropped for good\n")
                }
                else {
                  if (first.in[inew] == 0) 
                    first.in[inew] <- k
                  active <- c(active, inew)
                  if (!positive) 
                    Sign <- c(Sign, sign(Cvec[inew]))
                  else Sign <- c(Sign, abs(sign(Cvec[inew])))
                  action <- c(action, inew)
                  if (trace) 
                    cat("LARS Step", k, ":\t Variable", inew, 
                      "\tadded\n")
                }
            }
        }
        else action <- -dropid
        Gi1 <- backsolve(R, backsolvet(R, Sign))
        dropouts <- NULL
        if (type == "forward.stagewise") {
            directions <- Gi1 * Sign
            if (!all(directions > 0)) {
                if (use.Gram) {
                  nnls.object <- nnls.lars.faster(active, Sign, 
                    R, directions, Gram[active, active], trace = trace, 
                    use.Gram = TRUE, eps = eps)
                }
                else {
                  nnls.object <- nnls.lars.faster(active, Sign, 
                    R, directions, x[, active], trace = trace, 
                    use.Gram = FALSE, eps = eps)
                }
                positive <- nnls.object$positive
                dropouts <- active[-positive]
                action <- c(action, -dropouts)
                active <- nnls.object$active
                Sign <- Sign[positive]
                Gi1 <- nnls.object$beta[positive] * Sign
                R <- nnls.object$R
                C <- Cvec[-c(active, ignores)]
            }
        }
        A <- 1/sqrt(sum(Gi1 * Sign))
        w <- A * Gi1
        if (!use.Gram) 
            u <- (x[, active, drop = FALSE] %*% w)[, 1]
        if (length(active) >= min(n - 1, m - length(ignores))) {
            gamhat <- Cmax/A
        }
        else {
            if (use.Gram) {
                a <- my.drop(w %*% Gram[active, -c(active, ignores), 
                  drop = FALSE])
            }
            else {
                a <- (u %*% x[, -c(active, ignores), drop = FALSE])[1, 
                  ]
            }
            if (!positive) 
                gam <- c((Cmax - C)/(A - a), (Cmax + C)/(A + 
                  a))
            else gam <- (Cmax - C)/(A - a)
            gamhat <- min(gam[gam > eps], Cmax/A)
        }
        if (type == "lasso") {
            dropid <- NULL
            b1 <- beta[k, active]
            z1 <- -b1/w
            zmin <- min(z1[z1 > eps], gamhat)
            if (zmin < gamhat) {
                gamhat <- zmin
                drops <- z1 == zmin
            }
            else drops <- FALSE
        }
        beta[k + 1, ] <- beta[k, ]
        beta[k + 1, active] <- beta[k + 1, active] + gamhat * 
            w
        if (use.Gram) {
            Cvec <- Cvec - gamhat * Gram[, active, drop = FALSE] %*% 
                w
        }
        else {
            residuals <- residuals - gamhat * u
            Cvec <- (t(residuals) %*% x)[1, ]
        }
        Gamrat <- c(Gamrat, gamhat/(Cmax/A))
        arc.length <- c(arc.length, gamhat)
        if (type == "lasso" && any(drops)) {
            dropid <- seq(drops)[drops]
            for (id in rev(dropid)) {
                if (trace) 
                  cat("Lasso Step", k + 1, ":\t Variable", active[id], 
                    "\tdropped\n")
                R <- downdateR(R, id)
            }
            dropid <- active[drops]
            beta[k + 1, dropid] <- 0
            active <- active[!drops]
            Sign <- Sign[!drops]
        }
        if (!is.null(vn)) 
            names(action) <- vn[abs(action)]
        actions[[k]] <- action
        inactive <- im[-c(active, ignores)]
    }
    beta <- beta[seq(k + 1), ]
    dimnames(beta) <- list(paste(0:k), vn)
    if (trace) 
        cat("Computing residuals, RSS etc .....\n")
    residuals <- y - x %*% t(beta)
    beta <- scale(beta, FALSE, normx)
    RSS <- apply(residuals^2, 2, sum)
    R2 <- 1 - RSS/RSS[1]
    Cp <- ((n - k - 1) * RSS)/rev(RSS)[1] - n + 2 * seq(k + 1)
    object <- list(call = call, type = TYPE, R2 = R2, RSS = RSS, 
        Cp = Cp, actions = actions[seq(k)], entry = first.in, 
        Gamrat = Gamrat, arc.length = arc.length, Gram = if (use.Gram) Gram else NULL, 
        beta = beta, mu = mu, normx = normx, meanx = meanx)
    class(object) <- c("lars", "lars.pos")
    object
}
nnls.lars.faster <-
function (active, Sign, R, beta, Gram, eps = 1e-10, trace = FALSE, 
    use.Gram = TRUE) 
{
    if (!use.Gram) 
        x <- Gram
    M <- m <- length(active)
    im <- seq(m)
    positive <- im
    zero <- NULL
    while (m > 1) {
        zero.old <- c(m, zero)
        R.old <- downdateR(R, m)
        beta0 <- backsolve(R.old, backsolvet(R.old, Sign[-zero.old])) * 
            Sign[-zero.old]
        beta.old <- c(beta0, rep(0, length(zero.old)))
        if (all(beta0 > 0)) 
            break
        m <- m - 1
        zero <- zero.old
        positive <- im[-zero]
        R <- R.old
        beta <- beta.old
    }
    while (TRUE) {
        while (!all(beta[positive] > 0)) {
            alpha0 <- beta.old/(beta.old - beta)
            alpha <- min(alpha0[positive][(beta <= 0)[positive]])
            beta.old <- beta.old + alpha * (beta - beta.old)
            dropouts <- match(alpha, alpha0[positive], 0)
            for (i in rev(dropouts)) R <- downdateR(R, i)
            positive <- positive[-dropouts]
            zero <- im[-positive]
            beta0 <- backsolve(R, backsolvet(R, Sign[positive])) * 
                Sign[positive]
            beta <- beta.old * 0
            beta[positive] <- beta0
        }
        if (use.Gram) {
            w <- 1 - Sign * my.drop(Gram %*% (Sign * beta))
        }
        else {
            jw <- x %*% (Sign * beta)
            w <- 1 - Sign * my.drop(t(jw) %*% x)
        }
        if ((length(zero) == 0) || all(w[zero] <= 0)) 
            break
        add <- order(w)[M]
        if (use.Gram) {
            R <- updateR.faster(Gram[add, add], R, my.drop(Gram[add, 
                positive]), Gram = TRUE, eps = eps)
        }
        else {
            R <- updateR.faster(x[, add], R, x[, positive], Gram = FALSE, 
                eps = eps)
        }
        positive <- c(positive, add)
        zero <- setdiff(zero, add)
        beta0 <- backsolve(R, backsolvet(R, Sign[positive])) * 
            Sign[positive]
        beta[positive] <- beta0
    }
    if (trace) {
        dropouts <- active[-positive]
        for (i in dropouts) {
            cat("NNLS Step:\t Variable", i, "\tdropped\n")
        }
    }
    list(active = active[positive], R = R, beta = beta, positive = positive)
}
updateR.faster <-
function (xnew, R = NULL, xold, eps = .Machine$double.eps, Gram = FALSE) 
{
    xtx <- if (Gram) 
        xnew
    else sum(xnew^2)
    norm.xnew <- sqrt(xtx)
    if (is.null(R)) {
        R <- matrix(norm.xnew, 1, 1)
        attr(R, "rank") <- 1
        return(R)
    }
    Xtx <- if (Gram) 
        xold
    else (t(xnew) %*% xold)[1, ]
    r <- backsolvet(R, Xtx)
    rpp <- norm.xnew^2 - sum(r^2)
    rank <- attr(R, "rank")
    if (rpp <= eps) 
        rpp <- eps
    else {
        rpp <- sqrt(rpp)
        rank <- rank + 1
    }
    R <- cbind(rbind(R, 0), c(r, rpp))
    attr(R, "rank") <- rank
    R
}
.onLoad <-
function( libname, pkgname ) { ##.onAttach
    packageStartupMessage( "Loading MeDiChI, version ", VERSION, " (", DATE, ")\n", sep="" )
  }

VERSION <-
"0.4.1"
DATE <-
"Tue Jul 26 08:29:25 2011"
