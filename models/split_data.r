options(width=160)

# load library info
lib = read.csv(gzfile("../library/cytosol_tiles.csv.gz"))

# load degron scores
degron = read.csv("../score/degrons.csv")

# assign data set info to all measured tiles

i = match(degron$name,lib$name)
degron$set = lib[i,"set"]
degron$resi_first = lib[i,"resi_first"]
degron$resi_last = lib[i,"resi_last"]

# Split into test and training data
cns = c("name","dna","aa","n","degron_score","degron_std","abundance_score","abundance_std")
train = degron[which(degron$set=="train"),cns]
test  = degron[which(degron$set=="test"), cns]

# Check that AA tile are unique 
stopifnot( length(train$aa) == length(unique(train$aa)) )
stopifnot( length(test$aa) == length(unique(test$aa)) )

# Check no overlap
stopifnot(! any( test$aa %in% train$aa ))

# Get gene id's from tile name
splitlist = strsplit(train$name, "_")
train_prot = unique(sapply(splitlist, "[", 1))
splitlist = strsplit(test$name, "_")
test_prot  = unique(sapply(splitlist, "[", 1))
stopifnot(! any(test_prot %in% train_prot) )

print(sprintf("Loaded %d measurements and using %d unique AA tiles for training and storing %d (%.2f%%) for testing/hold-out (%d or %.2f%% of proteins)",
              nrow(degron), nrow(train), nrow(test), nrow(test)/(nrow(train)+nrow(test))*100,
	      length(test_prot), length(test_prot)/(length(test_prot)+length(train_prot))*100))

# Dump fasta file of both test and training data for feature generation
fasta_file = file("model_data.fasta")
i = which(degron$set %in% c("test","train"))
write(sprintf(">%s\n%s",degron[i,"name"],degron[i,"aa"]), fasta_file)
close(fasta_file)

# Dump test data in fasta format
fasta_file = file("test_data.fasta")
i = which(degron$set == "test")
write(sprintf(">%s\n%s",degron[i,"name"],degron[i,"aa"]), fasta_file)
close(fasta_file)

# Dump training data in csv format
write.csv(train, row.names=F, quote=F, file="degron_train.csv")
# write.csv(test,  row.names=F, quote=F, file="degron_test.csv")


# # Dump sequences
# write.table(test[,c("name","aa")], file="degron_test.seq", row.names=F, col.names=F, quote=F)
# write.table(train[,c("name","aa")], file="degron_train.seq", row.names=F, col.names=F, quote=F)

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

