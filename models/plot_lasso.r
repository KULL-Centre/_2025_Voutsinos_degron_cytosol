options(width=160)
library("gplots")

feat = read.csv("ordered_features.csv")
degron = read.csv(gzfile("degron_score_pred.csv.gz"))

print(sprintf("PAP     model C-degron mean %.3f, sd %.3f and median %.3f", mean(degron$pap_ct,na.rm=T), sd(degron$pap_ct,na.rm=T), median(degron$pap_ct,na.rm=T)))
print(sprintf("Human30 model C-degron mean %.3f, sd %.3f and median %.3f", mean(degron$h30_ct,na.rm=T), sd(degron$h30_ct,na.rm=T), median(degron$h30_ct,na.rm=T)))

library("glmnet")
load("cache/lasso_fit.rda")
h30_rsq = 0.656

# varianbles for plotting
rownames(feat) = feat$feature
coef_mat = t(as.matrix(feat[1:179,c("coef_lm","coef_h30")]))
coef_max = max(coef_mat, na.rm=T)
coef_min = min(coef_mat, na.rm=T)

# quartz(width=15, height=10)
jpeg("lasso_analysis.jpg", width=18, height=12, units="cm", res=300, pointsize=7)
layout(matrix(c(1,2,3,3,4,4), ncol=2, byrow=T), width=c(2,2), height=c(1.2,1,1))

# consider redoing with more param ticks above and non-significant parameters in gray in the back and most important in front
par(mar=c(5,5,2,3)+.1)
plot(0,0,col=0, xlim=c(0.0,0.67), ylim=c(-0.25,0.07), xlab="R-squared", ylab="Coefficients")

for (feature in rownames(lasso$beta)) {
    lines(lasso$dev.ratio, lasso$beta[feature,], col=8, lwd=3)
}

fl     = c("E","D","L","I","V","F","Y","W","E29:E30","K28:E29:E30","P27:P29:A30")
fl_col = c( 2,  2,  4,  4,  4,  5,  5,  5,        3,            3,            6)
for (fi in seq(length(fl))) {
    lines(lasso$dev.ratio, lasso$beta[fl[fi],], col=fl_col[fi], lwd=1)
}
legend("bottomleft", c("E, D","L, I, V","F, Y, W","E29:E30, K28:E29:E30","P27:P29:A30"), col=c(2,4,5,3,6), lty=1, lwd=1)
abline(v=h30_rsq, lty=2)
text(-0.09, 0.08, "A", xpd=T, font=2, cex=2)

# this should contain human30 ct and internal terms also
par(mar=c(5,5,2,3)+.1)
breaks = seq(-1,2,0.025)
i = which(degron$set=="train")
hs = hist(degron[i,"abundance_score"], breaks=breaks, plot=F)
hp = hist(degron[i,"h30"], breaks=breaks, plot=F)
hi = hist(degron[i,"h30"]-degron[i,"h30_ct"], breaks=breaks, plot=F)
hc = hist(degron[i,"h30_ct"], breaks=breaks, plot=F)
plot(0,0,col=0, xlim=c(-.3,1), ylim=c(0,7.5), xlab="Score", ylab="Density")
lines(hs$mids, hs$density, lwd=2, col=4)
lines(hi$mids, hi$density, lwd=2, col=3)
lines(hc$mids, hc$density*0.5, lwd=2, col=2)
lines(hp$mids, hp$density, lwd=2, col="darkorange")
legend("topright", c("FACS score","Predicted score","Composition","C-degron"), lty=1, lwd=2, col=c(4,"darkorange",3,2))
text(-.45, 7.5, "B", xpd=T, font=2, cex=2)

par(mar=c(2,4,5,2)+.1)
bp = barplot(coef_mat[,1:90], names.arg=rep(NA,90), beside=T, ylim=c(coef_min,coef_max))
bp_x = apply(bp, MARGIN=2, mean)
text(bp_x, coef_max, "|", xpd=T, cex=.7)
text(bp_x, coef_max+0.011, colnames(coef_mat)[1:90], adj=0, srt=60, xpd=T, cex=.9)
abline(h=0, lwd=.5)
text(-18, 0.15, xpd=T, "C", font=2, cex=2)

bp = barplot(coef_mat[,91:ncol(coef_mat)], names.arg=rep(NA,ncol(coef_mat)-90), beside=T, ylim=c(coef_min,coef_max))
bp_x = apply(bp, MARGIN=2, mean)
text(bp_x, coef_max, "|", xpd=T, cex=.7)
text(bp_x, coef_max+0.011, colnames(coef_mat)[91:ncol(coef_mat)], adj=0, srt=60, xpd=T, cex=.9)
abline(h=0, lwd=.5)

dev.off()
# quartz.save("lasso_analysis.png", type="png")
