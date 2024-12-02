options(width=160)
library("glmnet")

score_file = "degron_train.csv.gz"
print(sprintf("Loading degron scores from: %s",score_file))
degron = read.csv(gzfile(score_file))


###
### Amino acid composition
###
print(sprintf("Get amino acid composition of %d training tiles", nrow(degron)))

# the 20 natural amino acids in my favorite order
aa_one = strsplit("CDEKRHNQAGSTVMLIFYWP","")[[1]]

# get amino acid composition as a matrix
get_composition = function(s, aa_set=aa_one, aa_prefix=NA) {
    sl = strsplit(s, "")
    tl = lapply(sl, function(v){ table(v)[aa_set]  })
    if (is.na(aa_prefix)) {
        aa_colnames = aa_set
    } else {
        stopifnot( length(aa_prefix) == 1 )
        aa_colnames = sprintf("%s_%s",aa_prefix,aa_set)
    }
    cm = matrix(unlist(tl), ncol=20, byrow=T, dimnames=list(NULL,aa_colnames))
    cm[is.na(cm)] = 0
    return(cm)
}

# composition of full tiles
aa_comp = get_composition(degron$aa)

# check that all tiles have 30 natural amino acids
stopifnot(all( apply(aa_comp, MARGIN=1, sum) == 30 ))


###
### Logistic Regression
###
# composition of positions 1-25 in tiles
aa_comp25 = get_composition(substr(degron$aa,1,25))
stopifnot(all( apply(aa_comp25, MARGIN=1, sum) == 25 ))

score_col = "abundance_score"
degron_binary = degron[,score_col] < median(degron[,score_col])

lambdas = 10^seq(2, -5, -0.1)
# Alpha=0 is ridge/L2 and aplha=1 is lasso regularization
# For 'binomial', results are returned only for the class corresponding to the second level of the factor response
# Measure 'deviance' (actual deviance) is default, alternative 'mse' (squared loss), 'class' (miscalssification) or 'auc'
cv_lr25 = cv.glmnet(aa_comp25, degron_binary, alpha=0, family='binomial', lambda=lambdas, type.measure="auc", nfolds=10)

# select a regularization strength to use in the final model
lambda = 1e-5

plot(cv_lr25)
abline(v=log(lambda), lty=2, lwd=2)
plot(cv_lr25$glmnet.fit, xvar="lambda")
abline(v=log(lambda), lty=2, lwd=2)

lr25 = glmnet(aa_comp25, degron_binary, alpha=0, family='binomial', lambda=lambda, type.measure="auc")
# m1 = matrix(c(coef(cv_lr25,s=lambda,exact=T)[aa_one,],coef(lr25)[aa_one,]), ncol=length(aa_one), byrow=T, dimnames=list(c("cv_lr25","lr25"),aa_one))
# barplot(m1, beside=T, legend=T, las=2, args.legend=list(x="top"))
# quartz.save("human25_remake.png", type="png")

cv25_fit_resp = predict(lr25, newx=aa_comp25, s=lambda, exact=T, type="response")
lr25_fit_resp = predict(lr25, newx=aa_comp25, type="response")
rp_deg = cor(lr25_fit_resp, degron[,"degron_score"], method="pearson")
rp_facs = cor(lr25_fit_resp, degron[,"abundance_score"], method="pearson")
print(sprintf("Logistic regression pearson to degron score %.4f and facs score %.4f", rp_deg, rp_facs))


###
### Regression
###

plot_residuals = function(residuals, categories) {
    stopifnot(length(residuals) == length(categories))
    cat_level = unique(categories); cat_level = cat_level[order(cat_level, decreasing=F)]
    ncat = length(cat_level)
    ymax = 12/(max(residuals) - min(residuals))
    m_res = mean(residuals); sd_res = sd(residuals)
    hr = hist(residuals, breaks=100, freq=F, ylim=c(0,ymax), xlim=c(-6*sd_res,6*sd_res))
    xx = seq(min(residuals),max(residuals),0.01)
    lines(xx, dnorm(xx, m_res, sd_res), lwd=2)
    for (ci in seq(ncat)) {
        h = hist(residuals[which(categories==cat_level[ci])], breaks=100, plot=F)
        lines(h$mids, h$density, col=ci+1, lwd=2)
    }
    pop = table(categories)/length(categories)
    legend("topright", c(sprintf("Avg %.2f sd %.2f",m_res,sd_res),sprintf("Bin %d %.1f%%",seq(ncat),pop*100)), col=seq(ncat+1), lty=1, lwd=2)
    # legend("topright", sprintf("Q%d: %.3g",seq(5),quantile(residuals)))
}

plot_glm = function(y, fit, categories=NA) {
    # quartz(height=8, width=8)
    # par(mfrow=c(2,2))
    quartz(height=4, width=16)
    par(mfrow=c(1,4))
    rmse = sqrt(sum((y-fit$fitted.values)^2))
    rp = cor(fit$fitted.values, y, method="pearson", use="complete.obs")
    rs = cor(fit$fitted.values, y, method="pearson", use="complete.obs")
    s = sprintf("RMSE=%.2f, Pearson=%.3f Spearman=%.3f", rmse, rp, rs)
    breaks = seq(10,300,10)
    
    h1 = hist2d(fit$fitted.values, y, nbins=200, col=c("white", rev(gray.colors(length(breaks)-1)),"black"), breaks=c(0,breaks,1e12),
                xlab="Predicted score", ylab="Assay score", main=s)
    abline(c(0,1), col=2, lwd=2)

    s = sprintf("GLM with %s %s", fit$family$family, fit$family$link)
    h2 = hist2d(fit$linear.predictors, y, nbins=200, col=c("white", rev(gray.colors(length(breaks)-1)),"black"), breaks=c(0,breaks,1e12),
                xlab="Linear dimension", ylab="Assay score", main=s)
    points(fit$linear.predictors, fit$fitted.values, col=2, pch=20, cex=.2)

    param = coef(fit)
    stopifnot(names(param)[1] == "(Intercept)")
    s = sprintf("Intercept %.2f", param[1])
    barplot(param[2:length(param)], las=2, main=s)

    if (any(is.na(categories))) { categories = rep(1,length(y)) }
    plot_residuals(y-fit$fitted.values, categories)
}

###
### Regression with regularization
###

plot_glmnet = function(x, y, fit, i_lambda, mark_param=NA, categories=NA) {
    stopifnot(i_lambda <= length(fit$lambda))
    lambda = fit$lambda[i_lambda]
    
    quartz(height=4, width=16)
    par(mfrow=c(1,4))
    
    s = sprintf("GLMnet alpha=%.2f, opt=%s", fit$call$alpha, fit$call$type.measure)
    # plot(fit, xvar="lambda", label=T)
    # abline(v=log(lambda), lty=2)
    plot(fit, xvar="dev", label=T)
    abline(v=fit$dev.ratio[i_lambda], lty=2)
    title(main=s, line=2.5)
    legend("topleft", sprintf("Lambda = %.2g",lambda), lty=2)

    # error measures
    fitted.values = predict(fit, newx=x, s=lambda, exact=T)
    y_mean = mean(y)
    tss = sum((y-y_mean)^2)
    rss = sum((y - fitted.values)^2)
    r2 = 1 - rss/tss
    rp = cor(fitted.values, y, method="pearson", use="complete.obs")
    rs = cor(fitted.values, y, method="spearman", use="complete.obs")
    s = sprintf("R2=%.3f, RMSE=%.2f, Pearson=%.3f, Spearman=%.3f", r2, sqrt(rss), rp, rs)
    breaks = seq(10,300,10)
    
    h1 = hist2d(fitted.values, y, nbins=200, col=c("white", rev(gray.colors(length(breaks)-1)),"black"), breaks=c(0,breaks,1e12),
                xlab="Predicted score", ylab="Assay score", main=s)
    abline(c(0,1), col=2, lwd=2)

    param = coef(fit, s=lambda, exact=T)[,1]
    stopifnot(names(param)[1] == "(Intercept)")
    s = sprintf("Intercept %.2f", param[1])
    i_act = which(abs(param[2:length(param)]) > 1e-6) +1
    i_ina = setdiff(seq_along(param),c(1,i_act))
    print(sprintf("Using %d (%.0f%%):", length(i_act), length(i_act)/length(c(i_act,i_ina))*100))
    print(names(param)[i_act])
    print(sprintf("Not using %d:", length(i_ina)))
    print(names(param)[i_ina])
    barplot(param[i_act], main=s, las=2)
    # mark the offset of the non-regualized model
    if (mark_param %in% names(param) ) { i=which(mark_param==names(param)); abline(h=param[i,1], lty=2) }

    if (any(is.na(categories))) { categories = rep(1,length(y)) }
    plot_residuals(y-fitted.values, categories)

    return( names(param)[i_act] )
}


###
### PSSM terms
###
print(sprintf("Get PSSM of %d training tiles", nrow(degron)))

# plot heatmap function
heatmap = function(score_mat, mark_mat=NA, grad_range_min=-0.25, grad_range_max=0.25, n_scale_tics=11, filename="heatmap.jpg") {
    ns = abs(grad_range_min)*1000
    nd = abs(grad_range_max)*1000
    col_destab = "red"
    col_neutral = "white"
    col_stab = "cyan"
    col_grad = c(colorRampPalette( c(col_stab, col_neutral), space="rgb")(ns), colorRampPalette( c(col_neutral,col_destab), space="rgb")(nd))
    col_breaks = c(-100, seq(grad_range_min, grad_range_max, length.out=ns+nd-1), 100)

    # quartz(width=6, height=8)
    jpeg(filename, width=9, height=12, units="cm", res=300, pointsize=8)
    layout(matrix(c(1,1,3,2),ncol=2,byrow=F), width=c(5,1), height=c(2,3))
    par(mar=c(3,5,3,1)+.1)
    image(score_mat, xaxt="n", yaxt="n", zlim=c(grad_range_min,grad_range_max), col=col_grad, breaks=col_breaks)
    axis(1, seq(0,20,length=20)/20, aa_one, las=1, gap.axis=0)
    axis(3, seq(0,20,length=20)/20, aa_one, las=1, gap.axis=0)
    axis(2, seq(0,30,length=30)/30, labels=paste("Res ",seq(-30,-1,1)), las=2)

    if (! any(is.na(mark_mat))) {
        stopifnot(all( dim(score_mat) == dim(mark_mat) ))
	coords = expand.grid(seq(0,1,length.out=nrow(mark_mat)),seq(0,1,length.out=ncol(mark_mat)))
	points(coords[which(mark_mat),1], coords[which(mark_mat),2], pch=16)
    }

    par(mar=c(3,0.5,1,4))
    image(t(col_breaks), zlim=c(grad_range_min,grad_range_max), col=col_grad, breaks=col_breaks, xaxt="n", yaxt="n")
    axis(4, seq(0,1,length=n_scale_tics), seq(grad_range_min,grad_range_max,length=n_scale_tics), las=2)
    dev.off()
}

# calculate mean score of peptides that have a given aa at given position
nres = 30
score_col = "abundance_score"
splitlist = strsplit(degron$aa, "")
mean_score_mat = matrix(NA, nrow=length(aa_one), ncol=nres, dimnames=list(aa_one,seq(nres)))
n_mat = matrix(NA, nrow=length(aa_one), ncol=nres, dimnames=list(aa_one,NULL))
index_list = list()
for (pos in seq(nres)) {
    aa_pos = sapply(splitlist, "[", pos)
    index_list[[pos]] = list()
    for (aa in aa_one) {
        di = which(aa_pos == aa)
	# di = intersect(di, ihr_i2_deg)
        index_list[[pos]][[aa]] = di
        mean_score_mat[aa,pos] = mean(degron[di,score_col])
        n_mat[aa,pos] = length(di)
    }
}
avg_score = mean(degron[,score_col])

heatmap(mean_score_mat-avg_score, grad_range_min=-.1, grad_range_max=.1)

# function to count each AA ant each position
get_pssm = function(seqs, resi, aa=aa_one) {
    cns = paste0( rep(aa, times=length(resi)), rep(resi,each=length(aa)) )
    sl = strsplit(seqs, "")
    pssm = matrix(0, nrow=length(sl), ncol=length(aa)*length(resi), dimnames=list(NULL, cns))
    for (ri in resi) {
        d = data.frame(ri=seq(length(sl)), aa=sapply(sl, "[", ri))
        d = d[which(d$aa %in% aa),]
        d$aa_pos = paste0(d$aa, ri)
        d$ci = match(d$aa_pos,cns)
        pssm[matrix(c(d$ri,d$ci),nrow=nrow(d),byrow=F)] = 1
    }
    return(pssm)
}

pssm = get_pssm(degron$aa, seq(30))


###
### Pair couplings and features
###

# function to calculate coupling between pairs of amino acids at specific positions
calc_pair_coupling = function(aaij, scores, posi, posj, n_cut) {
    aai = aaij[1]
    aaj = aaij[2]
    pepi = index_list[[posi]][[aai]]
    pepj = index_list[[posj]][[aaj]]
    pepij = intersect(pepi, pepj)
    pepi = setdiff(pepi, pepij)
    pepj = setdiff(pepj, pepij)
    avg_score = mean(scores)
    if (length(pepij) > n_cut & length(pepi) > n_cut & length(pepj) > n_cut) {
        return( mean(scores[pepij])*avg_score / (mean(scores[pepi])*mean(scores[pepj])) )
    } else {
        return( NA )
    }
}

# all position pairs 30 x 29 / 2 = 435
resi_first = 1
resi_last = 30
pos_pairs = expand.grid(seq(resi_first,resi_last), seq(resi_first,resi_last))
pos_pairs = pos_pairs[which(pos_pairs$Var1 < pos_pairs$Var2),]

# all amino acid pairs 20 x 20 = 400
pair_couplings_df = expand.grid(aa_one,aa_one)

# calculate all pair couplings, evaluate 400 x 435 = 174,000 pairs
cut = 100
for (ppi in seq(nrow(pos_pairs))) {
    pp_id = paste0(pos_pairs[ppi,],collapse=":")
    pair_couplings_df[,pp_id] = apply(pair_couplings_df[,c(1,2)], MARGIN=1, calc_pair_coupling, scores=degron$abundance_score, posi=pos_pairs[ppi,1], posj=pos_pairs[ppi,2], n_cut=cut)
}

# get all pairs with calculated couplings, i.e. occurences above threshold
coupled_pairs = data.frame(name=NULL, coupling=NULL)
for (ppi in seq(nrow(pos_pairs))) {
    pp_id = paste0(pos_pairs[ppi,],collapse=":")
    name_df = cbind(matrix(pos_pairs[ppi,c("Var1","Var2")], byrow=T, ncol=2, nrow=400), pair_couplings_df[,c("Var1","Var2")])
    pair_names = apply(name_df, MARGIN=1, function(v){ sprintf("%s%s:%s%s",v[3],v[1],v[4],v[2]) })
    df = data.frame(name=pair_names, coupling=pair_couplings_df[,pp_id])
    coupled_pairs = rbind(coupled_pairs, df[which(! is.na(df$coupling)),])
}
coupled_pairs = coupled_pairs[order(coupled_pairs$coupling),]
write.csv(coupled_pairs, file="coupled_pairs.csv", quote=F, row.names=F)

save(pair_couplings_df, coupled_pairs, file="couplings_pair.rda")
# load("cache/couplings_pair.rda")
print(sprintf("Got %d (%.1f%%) calculated pairs couplings from positions %d-%d with more than %d occurences", nrow(coupled_pairs),
              nrow(coupled_pairs)/(nrow(pos_pairs)*nrow(pair_couplings_df))*100, resi_first, resi_last, cut))

pair_coupling_threshold = 0.2
degron_pairs = coupled_pairs[which(coupled_pairs$coupling < (1.0-pair_coupling_threshold)),"name"]
stable_pairs = coupled_pairs[which(coupled_pairs$coupling > (1.0+pair_coupling_threshold)),"name"]

get_pair_data = function(seqs, pair_names) {
    stopifnot(all(grepl(":", pair_names)))
    psl = strsplit(pair_names, ":")
    aal = lapply(psl, substr, 1, 1)
    posl = lapply(psl, function(v){ as.numeric(substr(v,2,nchar(v))) })

    pairs = matrix(0, nrow=length(seqs), ncol=length(pair_names), dimnames=list(NULL,pair_names))
    sl = strsplit(seqs, "")
    for (pi in seq_along(pair_names)) {
        aa1_matching = sapply(sl,"[",posl[[pi]][1]) == aal[[pi]][1]
        aa2_matching = sapply(sl,"[",posl[[pi]][2]) == aal[[pi]][2]
        si = which(aa1_matching & aa2_matching)
	# pairs[matrix(c(si,rep(pi,times=length(si))), nrow=length(si), ncol=2, byrow=F)] = 1
	if (length(si) > 0) {
	    pairs[si,pi] = 1
	} else {
	    print(sprintf("WARNING: Pair %s not found",pair_names[pi]))
	}
    }
    return(pairs)
}

print(sprintf("Get pair features for %d degron and %d stable out of %d (%.1f%%) calculated pairs with couplings stronger than %.1f%%",
              length(degron_pairs), length(stable_pairs), nrow(coupled_pairs), length(c(degron_pairs,stable_pairs))/nrow(coupled_pairs)*100,
	      pair_coupling_threshold*100))
pairs = get_pair_data(degron$aa, c(degron_pairs,stable_pairs))


###
### Triplet couplings and features
###

# function to calculate coupling between pairs of amino acids at specific positions
calc_triplet_coupling = function(aaijk, scores, posi, posj, posk, n_cut) {
    aai = aaijk[1]
    aaj = aaijk[2]
    aak = aaijk[3]
    pepi = index_list[[posi]][[aai]]
    pepj = index_list[[posj]][[aaj]]
    pepk = index_list[[posk]][[aak]]
    pepijk = intersect(pepi, intersect(pepj, pepk))
    pepi = setdiff(pepi, pepijk)
    pepj = setdiff(pepj, pepijk)
    pepk = setdiff(pepk, pepijk)
    avg_score = mean(scores)
    if (length(pepijk) > n_cut & length(pepi) > n_cut & length(pepj) > n_cut & length(pepk) > n_cut) {
        return( mean(scores[pepijk])*avg_score*avg_score / (mean(scores[pepi])*mean(scores[pepj])*mean(scores[pepk])) )
    } else {
        return( NA )
    }
}


# all triplets in CT
resi_first = 24
resi_last = 30
pos_trip = expand.grid(seq(resi_first,resi_last), seq(resi_first,resi_last), seq(resi_first,resi_last))
pos_trip = pos_trip[which(pos_trip$Var1 < pos_trip$Var2 & pos_trip$Var2 < pos_trip$Var3),]

# all amino acid triplets 20 x 20 x 20 = 8000
triplet_couplings_df = expand.grid(aa_one,aa_one,aa_one)

# calculate couplings of all triplets
cut = 100
for (ppi in seq(nrow(pos_trip))) {
    pp_id = paste0(pos_trip[ppi,],collapse=":")
    triplet_couplings_df[,pp_id] = apply(triplet_couplings_df[,c(1,2,3)], MARGIN=1, calc_triplet_coupling, scores=degron$abundance_score, posi=pos_trip[ppi,1], posj=pos_trip[ppi,2], posk=pos_trip[ppi,3], n_cut=cut)
    print(sprintf("Done triplet %d of %d %s with %d couplings calculated", ppi, nrow(pos_trip), pp_id, sum(! is.na(triplet_couplings_df[,pp_id]))))
}

# get all triplets with calculated couplings, i.e. occurences above threshold
coupled_triplets = data.frame(name=NULL, coupling=NULL)
for (ppi in seq(nrow(pos_trip))) {
    pp_id = paste0(pos_trip[ppi,],collapse=":")
    name_df = cbind(matrix(pos_trip[ppi,c("Var1","Var2","Var3")], byrow=T, ncol=3, nrow=8000), triplet_couplings_df[,c("Var1","Var2","Var3")])
    trip_names = apply(name_df, MARGIN=1, function(v){ sprintf("%s%s:%s%s:%s%s",v[4],v[1],v[5],v[2],v[6],v[3]) })
    df = data.frame(name=trip_names, coupling=triplet_couplings_df[,pp_id])
    coupled_triplets = rbind(coupled_triplets, df[which(! is.na(df$coupling)),])
}
coupled_triplets = coupled_triplets[order(coupled_triplets$coupling),]
write.csv(coupled_triplets, file="coupled_triplets.csv", quote=F, row.names=F)

save(triplet_couplings_df, coupled_triplets, file="couplings_trip.rda")
# load("cache/couplings_trip.rda")

n_all_trip = nrow(pos_trip)*nrow(triplet_couplings_df)
print(sprintf("Got %d of %d (%.1f%%) calculated triplet couplings from positions %d-%d with more than %d occurences", nrow(coupled_triplets),
              n_all_trip, nrow(coupled_triplets)/n_all_trip*100, resi_first, resi_last, cut))

trip_coupling_threshold = 0.2
degron_triplets = coupled_triplets[which(coupled_triplets$coupling < (1.0-trip_coupling_threshold)),"name"]
stable_triplets = coupled_triplets[which(coupled_triplets$coupling > (1.0+trip_coupling_threshold)),"name"]

get_triplet_data = function(seqs, triplet_names) {
    stopifnot(all(grepl(":", triplet_names)))
    tsl = strsplit(triplet_names, ":")
    aal = lapply(tsl, substr, 1, 1)
    posl = lapply(tsl, function(v){ as.numeric(substr(v,2,nchar(v))) })

    pairs = matrix(0, nrow=length(seqs), ncol=length(triplet_names), dimnames=list(NULL,triplet_names))
    ssl = strsplit(seqs, "")
    for (ti in seq_along(triplet_names)) {
        aa1_matching = sapply(ssl,"[",posl[[ti]][1]) == aal[[ti]][1]
        aa2_matching = sapply(ssl,"[",posl[[ti]][2]) == aal[[ti]][2]
        aa3_matching = sapply(ssl,"[",posl[[ti]][3]) == aal[[ti]][3]
        si = which(aa1_matching & aa2_matching & aa3_matching)
	# pairs[matrix(c(si,rep(ti,times=length(si))), nrow=length(si), ncol=2, byrow=F)] = 1
	if (length(si) > 0) {
	    pairs[si,ti] = 1
	} else {
	    print(sprintf("WARNING: Pair %s not found", triplet_names[ti]))
	}
    }
    return(pairs)
}

print(sprintf("Get triplet features for %d degron and %d stable out of %d (%.1f%%) calculated triplets with coupling stronger than %.1f%%",
              length(degron_triplets), length(stable_triplets), nrow(coupled_triplets), length(c(degron_triplets,stable_triplets))/nrow(coupled_triplets)*100,
	      trip_coupling_threshold*100))
triplets = get_triplet_data(degron$aa, c(degron_triplets,stable_triplets))

save(aa_comp, pssm, pairs, triplets, file="features.rda")
# load("cache/features.rda")

###
### Lasso regression
###
library("glmnet")

x = cbind(aa_comp, pssm, pairs, triplets)
y = degron[,"abundance_score"]
lambdas = 10^seq(-1, -6, -0.01)

lasso = glmnet(x, y, alpha=1, family='gaussian', lambda=lambdas, type.measure="mse")
save(lasso, file="lasso_fit.rda")
# load("cache/lasso_fit.rda")

# sort param in order of apprearence in lasso regression
param_mat = lasso$beta
first_li = apply(param_mat, MARGIN=1, function(v){ v=c(v,1); min(ifelse(abs(v) > 1e-9, seq_along(v), NA), na.rm=T) })
param_mat = param_mat[order(first_li, decreasing=F),]

# select after importance until insignificant
feat_sig = rownames(param_mat)[1:179]
fit_lm = lm(y ~ x[,feat_sig])
sum_lm = summary(fit_lm)
coef_lm = coef(fit_lm)[2:length(coef(fit_lm))]
names(coef_lm) = sub("x[, feat_sig]","",names(coef_lm), fixed=T)
rp = cor(y, fit_lm$fitted.values, method="pearson", use="complete.obs")
mse = mean((y-fit_lm$fitted.values)^2)
print(sprintf("Using %d parameters (excl. intercept) R-sq %.4f, Pearson %.4f and RMSE %.4f", length(coef_lm), sum_lm$r.squared, rp, sqrt(mse)))

# select feature for a human30 linear model
feat_h30 = rownames(param_mat)[1:55]
fit_h30 = lm(y ~ x[,feat_h30])
sum_h30 = summary(fit_h30)
coef_h30 = coef(fit_h30)[2:length(coef(fit_h30))]
names(coef_h30) = sub("x[, feat_h30]","",names(coef_h30), fixed=T)
rp = cor(y, fit_h30$fitted.values, method="pearson")
mse = mean((y-fit_h30$fitted.values)^2)
print(sprintf("Using %d parameters (excl. intercept) R-sq %.4f, Pearson %.4f and RMSE %.4f", length(coef_h30), sum_h30$r.squared, rp, sqrt(mse)))

param_df = data.frame(feature=rownames(param_mat), coef_lasso=param_mat[,ncol(param_mat)])
param_df$coef_lm=NA
param_df[match(names(coef_lm),param_df$feature),"coef_lm"] = coef_lm
param_df$coef_h30=NA
param_df[match(names(coef_h30),param_df$feature),"coef_h30"] = coef_h30
write.csv(param_df, file="ordered_features.csv", quote=F, row.names=F)

