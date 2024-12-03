options(width=160)

# load library info
lib = read.csv(gzfile("../library/cytosol_tiles.csv.gz"))

# load degron scores
degron = read.csv("../score/degrons.csv")

# annotate scored tiles with library info
i = match(degron$name,lib$name)
degron$set = lib[i,"set"]
# degron$resi_first = lib[i,"resi_first"]
# degron$resi_last = lib[i,"resi_last"]

# Split into test and training data
cns = c("name","dna","aa","n","degron_score","degron_std","abundance_score","abundance_std")
train = degron[which(degron$set=="train"),cns]
test  = degron[which(degron$set=="test"), cns]

# Check that AA tile are unique 
stopifnot( length(train$aa) == length(unique(train$aa)) )
stopifnot( length(test$aa) == length(unique(test$aa)) )

# Check no overlap
stopifnot(! any( test$aa %in% train$aa ))

# Check that no gene id's are overlapping using tile name
splitlist = strsplit(train$name, "_")
train_prot = unique(sapply(splitlist, "[", 1))
splitlist = strsplit(test$name, "_")
test_prot  = unique(sapply(splitlist, "[", 1))
stopifnot(! any(test_prot %in% train_prot) )

print(sprintf("Loaded %d measurements. Using %d unique AA tiles for training and %d (%.2f%%) for testing/hold-out (%d or %.2f%% of proteins)",
              nrow(degron), nrow(train), nrow(test), nrow(test)/(nrow(train)+nrow(test))*100,
	      length(test_prot), length(test_prot)/(length(test_prot)+length(train_prot))*100))


# Dump fasta file of all measured tiles for prediction
fasta_file = file("all_data.fasta")
write(sprintf(">%s\n%s",degron$name,degron$aa), fasta_file)
close(fasta_file)

# # Dump test data in fasta format
# fasta_file = file("test_data.fasta")
# i = which(degron$set == "test")
# write(sprintf(">%s\n%s",degron[i,"name"],degron[i,"aa"]), fasta_file)
# close(fasta_file)

# Dump training data in csv format
i = which(degron$set == "train")
write.csv(degron[i,], row.names=F, quote=F, file="degron_train.csv")

# Load predictions of all data
cn = c("name", "seq", "score_ct", "score", "resi_first", "resi_last", "aa", "resi")
pap = read.table(gzfile("all_data_pap.txt.gz")); colnames(pap) = cn
h30 = read.table(gzfile("all_data_h30.txt.gz")); colnames(h30) = cn
h25 = read.table(gzfile("all_data_h25.txt.gz")); colnames(h25) = cn
i_pap = match(pap$name,degron$name)
degron[i_pap,"pap_ct"] = pap$score_ct
degron[i_pap,"pap"] = pap$score
i_h30 = match(h30$name,degron$name)
degron[i_h30,"h30_ct"] = h30$score_ct
degron[i_h30,"h30"] = h30$score
i_h25 = match(h25$name,degron$name)
degron[i_h25,"h25"] = h25$score

# Dump training data in csv format
write.csv(degron, row.names=F, quote=F, file="degron_score_pred.csv")


# # Plot scores
# quartz(height=5, width=15)
# par(mfrow=c(1,3))

# ha =  hist(degron$degron_score, breaks=100, plot=F)
# htr = hist(train$degron_score, breaks=100, plot=F)
# hte = hist(test$degron_score, breaks=100, plot=F)
# plot(0,0,col=0, xlim=c(-0.2,1.2), ylim=c(0,4), xlab="Degron score", ylab="Density")
# lines(ha$mids, ha$density, lwd=2, col=1)
# lines(hte$mids, hte$density, lwd=2, col=2)
# lines(htr$mids, htr$density, lwd=2, col=3)
# legend("top", c("All", "Test", "Train"), lwd=2, col=c(1,2,3))

# ha =  hist(degron$abundance_score, breaks=100, plot=F)
# htr = hist(train$abundance_score, breaks=100, plot=F)
# hte = hist(test$abundance_score, breaks=100, plot=F)
# plot(0,0,col=0, xlim=c(0,1), ylim=c(0,7), xlab="FACS score", ylab="Density")
# lines(ha$mids, ha$density, lwd=2, col=1)
# lines(hte$mids, hte$density, lwd=2, col=2)
# lines(htr$mids, htr$density, lwd=2, col=3)
# legend("top", c("All", "Test", "Train"), lwd=2, col=c(1,2,3))

# ha =  hist(log2(degron$abundance_score), breaks=100, plot=F)
# htr = hist(log2(train$abundance_score), breaks=100, plot=F)
# hte = hist(log2(test$abundance_score), breaks=100, plot=F)
# plot(0,0,col=0, xlim=c(-13,0), ylim=c(0,.5), xlab="log-FACS score", ylab="Density")
# lines(ha$mids, ha$density, lwd=2, col=1)
# lines(hte$mids, hte$density, lwd=2, col=2)
# lines(htr$mids, htr$density, lwd=2, col=3)
# legend("top", c("All", "Test", "Train"), lwd=2, col=c(1,2,3))

# quartz.save("test_train_distributions.png", type="png")

