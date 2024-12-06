options(width=160)
degron = read.csv(gzfile("degron_score_pred.csv.gz"))

i = which(degron$set=="test")

plot_scat = function(x, y, i_mark=c(), ...) {
    rp = cor(x, y, method="pearson")
    rs = cor(x, y, method="spearman")
    print(sprintf("Pearson %.4f, Spearman %.4f",rp,rs))
    plot(x, y, pch=16, cex=.2, ...)
    if (length(i_mark) > 0) {
        points(x[i_mark], y[i_mark], pch=16, cex=1.6, col=c(2,3,4))
	text(x[i_mark], y[i_mark], seq(length(i_mark)), cex=.7)
    }
}

sel = c(1734, 553, 2973)
pap_int = degron[i[sel],"pap"]-degron[i[sel],"pap_ct"]
cbind(degron[i[sel],c("name","aa","abundance_score","abundance_std","pap_ct","pap","h30_ct","h30")], pap_int)

# quartz(width=12, height=6)
jpeg("pap.jpg", width=12, height=6, units="cm", res=300, pointsize=7)

par(mfrow=c(1,2), mar=c(5,5,1,1)+.1)
plot_scat(degron[i,"pap"], degron[i,"abundance_score"], xlab="Predicted score", ylab="Abundance score", i_mark=sel)
plot_scat(degron[i,"pap"]-degron[i,"pap_ct"], degron[i,"abundance_score"], xlab="Predicted internal score", ylab="Abundance score", i_mark=sel)

dev.off()