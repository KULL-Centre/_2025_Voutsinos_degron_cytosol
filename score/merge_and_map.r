options(width=160, digits=4, stringsAsFactors=F)

args = commandArgs(trailingOnly=TRUE)
if (interactive()) {
    samples_file = "samples.csv"
    files = c("O3r1_S11_L001_counts.txt","O3r1_S11_L002_counts.txt","O3r1_S11_L003_counts.txt")
    # files = c("FL-pDNA-4_S70_L001_counts.txt", "FL-pDNA-4_S70_L002_counts.txt", "FL-pDNA-4_S70_L003_counts.txt", "FL-pDNA-4_S70_L004_counts.txt")
} else if (length(args) < 2) {
    print("")
    print("usage: Rscript merge_and_map.r  <samples.csv>  <counts1.txt>  [counts2.txt  ...]")
    quit(save="no")
} else {
    samples_file = args[1]
    files = args[2:length(args)]
}
print(sprintf("Args: %s",paste0(args, collapse=" ")))


##
## Settings
##
settings = list()

# number of read counts required to consider DNA covered in summary table
settings$coverage_min_counts = 20

# store unmapped reads for later analysis
settings$store_unmapped = TRUE
settings$unmapped_min_reads = 5
settings$unmapped_good_lengths = c(90,72,66)


###
###  Libraries
###

lib = read.csv(gzfile("../library/cytosol_tiles.csv.gz"))
raw = list()
raw[["even1"]] = lib[which(lib$even1),c("name","dna_rm_restrict","aa")]
raw[["even2"]] = lib[which(lib$even2),c("name","dna_rm_restrict","aa")]
raw[["odd3"]]  = lib[which(lib$odd3),c("name","dna_rm_restrict","aa")]
raw[["odd4"]]  = lib[which(lib$odd4),c("name","dna_rm_restrict","aa")]
raw[["ct1"]]   = lib[which(lib$ct1),c("name","dna_rm_restrict","aa")]

# check that DNA library is unique and translate to amino acid sequences
for (libname in names(raw)) {
    # DNA sequences should be unique per library
    colnames(raw[[libname]]) = c("name","dna","aa")
    stopifnot(length(raw[[libname]][,"dna"]) == length(unique(raw[[libname]][,"dna"])))    
    print(sprintf("Check if library %s contains unique AA sequences: %s",libname,length(raw[[libname]][,"aa"]) == length(unique(raw[[libname]][,"aa"]))))
}


###
###  Control sequences
###
i_ctrl = which(lib$pool == "control")
ctrl = data.frame(names=lib[i_ctrl,"name"], dna=lib[i_ctrl,"dna"], aa=lib[i_ctrl,"aa"])

# check that all library files have the same number of columns
name_cols = colnames(raw[[1]])
stopifnot(all(sapply(raw,ncol) == length(name_cols)))

# Read information about samples
samples = read.csv(samples_file)

###
### Read count files
###
raw_read_counts = list()
if (settings$store_unmapped) { unmapped_raw = list() }

# for each file, map counts to library and put them in a data frame
for (file in files) {
    # check filename and extract file_id
    fl = strsplit(file, "/")[[1]]
    stopifnot(substr(fl[length(fl)], nchar(fl[length(fl)])-10, nchar(fl[length(fl)])) == "_counts.txt")
    file_id = substr(fl[length(fl)], 1, nchar(fl[length(fl)])-11)

    # determine sample index
    if ( grepl("_L00",file_id) ) {
        si = which(samples$file == substr(file_id, 1, nchar(file_id)-5))
    } else {
	si = which(samples$file == file_id)
    }
    stopifnot(length(si) == 1)

    # library of sample
    libname = samples[si,"lib"]
    stopifnot(libname %in% names(raw))
    stopifnot(! file_id %in% colnames(raw[[libname]]))

    # read counts file: col 1 should be the DNA sequence, col 1 the counts
    cf = read.table(file)
    colnames(cf) = c("dna","counts","rpm")

    # map file read counts to library
    i_mapped = match(raw[[libname]][,"dna"],cf$dna)
    raw[[libname]][,file_id] = cf[i_mapped,"counts"]
    
    # set un-observed variants to zero counts
    i_na = is.na( raw[[libname]][,file_id] )
    raw[[libname]][i_na,file_id] = 0

    # store number of unique variants and total read counts
    n_unq_dna = nrow(cf)
    n_raw_counts = sum(cf$counts)
    raw_read_counts[[file_id]] = c(n_unq_dna, n_raw_counts)

    # Report
    n_mapped_dna = sum(raw[[libname]][,file_id] > 0)
    n_mapped_counts = sum(raw[[libname]][,file_id])
    # n_filtered_counts = sum(cf$counts)
    # print(sprintf("Mapped %d dna seq covering %.2f%% of library and %d of %d length-filtered read counts (%.2f%%)",
    #               n_mapped_dna, n_mapped_dna/nrow(raw[[libname]])*100,
    # 		  n_mapped_counts, n_filtered_counts, n_mapped_counts/n_filtered_counts*100))

    print(sprintf("Mapped %d out of %d raw counts (%.2f%%) from %s",n_mapped_counts,n_raw_counts,n_mapped_counts/n_raw_counts*100,file_id))

    # store unmapped redas for later analysis
    if (settings$store_unmapped) {
        i_above_threshold = which(cf$counts >= settings$unmapped_min_reads)
        unmapped_raw[[file_id]] = cf[setdiff(i_above_threshold,i_mapped),c("dna","counts")]
    }
    
    # print(table(unlist(strsplit(cf[,"dna"],""))))
}

# merge lane counts
counts = list()
unmapped_counts = list()
lane_names = list()
samples$obs = FALSE
for (libname in names(raw)) {
    if (ncol(raw[[libname]]) <= length(name_cols)) {
        print(sprintf("No data for %s - library will not be considered further",libname))
	next
    }
    
    # init counts data frame
    counts[[libname]] = raw[[libname]][,name_cols]
    
    # column names of columns with reads
    cn = colnames(raw[[libname]])[(length(name_cols)+1):ncol(raw[[libname]])]
    
    # column names with lane info removed
    cn_nolane = unique( sapply(cn, function(s){ if (grepl("_L00",s)) {substr(s,1,nchar(s)-5)} else {s} }) )
    
    for (sample_name in cn_nolane) {
        stopifnot(sample_name %in% samples$file)
	si = which(samples$file==sample_name)
	stopifnot(! samples[si,"obs"])
	samples[si,"obs"] = TRUE
	
        sample_lane_names = cn[which(grepl(sample_name, cn))]
	counts[[libname]][,sample_name] = apply(raw[[libname]][,sample_lane_names], MARGIN=1, sum)

        # consider merging raw_read_counts[[file_id]] = c(n_unq_dna, n_raw_counts)

	# report correlations of raw[[libname]][,sample_lane_names]
	lane_names[[sample_name]] = sample_lane_names
	print(sprintf("Merging %d lanes into %s", length(sample_lane_names), sample_name))

        # merge lanes for unmapped read counts
	if (settings$store_unmapped) {
            sample_unmapped_dna = unique(unlist(sapply(unmapped_raw[sample_lane_names], "[", "dna")))

            unmapped_counts[[sample_name]] = data.frame(dna = sample_unmapped_dna)
	    for (sn in sample_lane_names) {
	        unmapped_counts[[sample_name]][,sn] = unmapped_raw[[sn]][match(sample_unmapped_dna,unmapped_raw[[sn]][,"dna"]),"counts"]
	    }
	    # print(sprintf("Merged lanes for %s resulting in %d unmapped reads from %d unique DNA sequences",
	    #               sample_name, sum(unmapped_counts[[sample_name]][,"sum"]), length(sample_unmapped_dna)))
	}
    }
}
# print(sprintf("Done lane merging, distribution of lanes per sample:"))
# print(table(unlist(n_lanes)))

# counts per sample, info on lanes and unmapped reads not present
save(counts, samples, name_cols, ctrl, settings, file="counts.rda")

# data on unmapped reads and lanes
if (settings$store_unmapped) {
    save(unmapped_counts, raw_read_counts, lane_names, file="counts_unmapped.rda")
}

# table of read counts
all_lib_dna = unique(unlist(sapply(raw, "[", "dna")))
counts_summary = data.frame(sample=NULL, lib=NULL, coverage = NULL, ctrl_cover=NULL,
                            all=NULL, mapped=NULL, controls=NULL, mapped_other_lib=NULL, unmapped=NULL, unmapped_bad_len=NULL)

for (sn in sort(names(unmapped_counts))) {
    si = which(sn == samples$file)
    libname = samples[si,"lib"]
    sample_lane_names = lane_names[[sn]]

    i_ctrl = which(counts[[libname]][,"name"] %in% ctrl$name)
    stopifnot(length(i_ctrl) == nrow(ctrl))

    mapped_other_lib = 0
    unmapped = 0
    unmapped_bad_len = 0

    if (settings$store_unmapped) {
        # unmapped DNA that mappes to other libs
        imo = which(unmapped_counts[[sn]][,"dna"] %in% all_lib_dna)
        unmapped_counts[[sn]][,"sum"] = apply(unmapped_counts[[sn]][,sample_lane_names], MARGIN=1, sum, na.rm=T)
        mapped_other_lib = sum(unmapped_counts[[sn]][imo,"sum"])

        # completely unmapped DNA
        ium = setdiff(seq(nrow(unmapped_counts[[sn]])), imo)
        completely_unmapped = sum(unmapped_counts[[sn]][ium,"sum"])
        unmapped = sum(unmapped_counts[[sn]][ium,"sum"])

        # unmapped DNA with bad length, i.e. length not in
        unmapped_counts[[sn]][,"length"] = sapply(unmapped_counts[[sn]][,"dna"], nchar)
        ibl = ium[ which(! unmapped_counts[[sn]][ium,"length"] %in% settings$unmapped_good_lengths) ]
        unmapped_bad_len = sum(unmapped_counts[[sn]][ibl,"sum"])
    }
    
    df = data.frame(sample = sn, lib = libname,
                    coverage = sum(counts[[libname]][,sn] >= settings$coverage_min_counts)/nrow(counts[[libname]])*100,
		    ctrl_cover = sum(counts[[libname]][i_ctrl,sn] >= settings$coverage_min_counts)/nrow(ctrl)*100,
                    all = sum(sapply(raw_read_counts[sample_lane_names], "[", 2)),
                    mapped = sum(counts[[libname]][,sn]),
	            controls = sum(counts[[libname]][i_ctrl,sn]),
	            mapped_other_lib = mapped_other_lib,
		    unmapped = unmapped,
		    unmapped_bad_len = unmapped_bad_len)
    counts_summary = rbind(counts_summary, df)
}
print("")
print("Read counts summary")
print(counts_summary)
write.csv(counts_summary, file="counts_summary.csv")

pct_summary = counts_summary
cns = colnames(counts_summary)[5:ncol(counts_summary)]
for (ri in seq(nrow(pct_summary))) {
    pct_summary[ri,cns] = pct_summary[ri,cns]/pct_summary[ri,"all"]*100
}
print("")
print("Read percent summary")
print(pct_summary)

