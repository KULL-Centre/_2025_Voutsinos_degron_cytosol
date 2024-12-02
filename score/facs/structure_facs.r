options(width=160)

args = commandArgs(trailingOnly=TRUE)
if (interactive()) {
    infiles = c("data/export_81835_Q1 BV421-A- , PE-Texas Red-A+.csv")
} else if (length(args) < 1) {
    print("")
    print("usage: Rscript structure_facs.r  <facs1.csv>  [facs2.csv  ...]")
    quit(save="no")
} else {
    infiles = args
}
# print(sprintf("Args: %s",paste0(args, collapse=" ")))


parse_facs_csv = function(filename, n_cols=20) {
    d = read.table(filename, sep=",", row.names=NULL, fill=T, col.names=sprintf("col%.02d",seq(n_cols)))

    # keep reading until an excess of columns are read (read.table(fill=T)) only looks at first 5 lines)
    while (! all(is.na(d[,n_cols]) | n_cols > 200) ) {
        n_cols = n_cols+10
        d = read.table(filename, sep=",", row.names=NULL, fill=T, col.names=sprintf("col%.02d",seq(n_cols)))
    }

    facs = list()
    facs$sample = d[which(grepl("TUBE", d$col01)), "col02"]
    facs$date = d[which(grepl("DATE", d$col01)), "col02"]

    # find the number of rows with settings assuming a header row. Look for the row from which everything in column 1 is numeric
    # print("Search end of settings lines")
    settings_rows = 0
    while (suppressWarnings( any(is.na(as.numeric(d[(settings_rows+2):nrow(d),"col01"]))) & settings_rows < 500)) {
        settings_rows = settings_rows +1
    }
    stopifnot(settings_rows < 500)

    facs$events = read.table(filename, sep=",", skip=settings_rows, header=T)
    return(facs)
}


##
## Settings
##
settings = list()
settings$infiles = infiles

facs_events = list()
facs_table = list()
for (infile in infiles) {
    print(sprintf("Reading %s",infile))
    facsfile = parse_facs_csv(infile)
    print(sprintf("    Found sample %s from %s with %d events and %d channels",facsfile$sample,facsfile$date,nrow(facsfile$events),ncol(facsfile$events)))
    name = facsfile$sample
    if ( name %in% names(facs_events)) {
        i = 1
	name = sprintf("%s_%02d", facsfile$sample, i)
	while (name %in% names(facs_events)) {
	    i = i+1
	    name = sprintf("%s_%02d", facsfile$sample, i)
	}
	print(sprintf("    Renaming %s to %s", facsfile$sample, name))
    }
    facs_table[[name]] = c(facsfile$sample, facsfile$date, nrow(facsfile$events), ncol(facsfile$events))
    facs_events[[name]] = facsfile$events
}

# save(facs_events, facs_table, settings, file="facs_raw.rda")

# channels we are using
settings$cn_gfp = "GFP.A"
settings$cn_cherry = "PE.Texas.Red.A"
settings$cn_ratio = "GFP.mCherry....GFP.A.PE.Texas.Red.A"

# only store used channels
cns = c(settings$cn_gfp, settings$cn_cherry, settings$cn_ratio)
for (sample in names(facs_events)) {
    stopifnot(all( cns %in% colnames(facs_events[[sample]]) ))
    facs_events[[sample]] = facs_events[[sample]][,cns]
}

# reformat as table
facs_table = data.frame(name = names(facs_table),
                        sample = sapply(facs_table, "[[", 1),
                        date = sapply(facs_table, "[[", 2),
                        events = sapply(facs_table, "[[", 3))
write.csv(facs_table, quote=F, row.names=F, file="facs_table.csv")

# used profiles
settings$sel = list()
settings$sel[["even1"]] = "E1"      # File 26-10-2022_E1_007.fcs    with 4.32 E5 events - alt E1b, E1b_02
settings$sel[["even2"]] = "E2"      # File 07-12-2022_E2_007.fcs    with 4.51 E5 events - alt E2b_01, E2b_02
settings$sel[["odd3"]]  = "O3b_01"  # File 26-10-2022_O3b_006.fcs   with 395055  events - alt O3b, O3_01
settings$sel[["odd4"]]  = "O4"      # File 27-10-2022_O4_005.fcs    with 4.15 E5 events - alt O4b, O4_02, O4b_02
settings$sel[["ct1"]]   = "CT1 I"   # File 16-12-2022_CT1 I_005.fcs with 4.67 E5 events - alt CT1 IIb

facs_set = settings
save(facs_events, facs_table, facs_set, file="facs.rda")

# plot all distributions
hl = list()
x_min = 0; x_max = 0; y_max = 0
for (sample in names(facs_events)) {
    h = hist(facs_events[[sample]][,settings$cn_ratio], breaks=200, plot=F)
    x_min = min(c(x_min, h$mids))
    x_max = max(c(x_max, h$mids))
    y_max = max(c(y_max, h$density))
    hl[[sample]] = h
}
quartz(width=12, height=8)
plot(0,0,col=0, xlim=c(-0.1,1.1), ylim=c(0,10), xlab="GFP/mCherry", ylab="Density")
i = 0
for (sample in names(facs_events)) {
    i = i+1
    lines(hl[[sample]]$mids, hl[[sample]]$density, col=i, lty=i%/%8+1)
}
legend("topright", names(hl), col=seq(i), lty=seq(i)%/%8+1, ncol=4, cex=.8)
quartz.save("all_distributions.png", type="png")

