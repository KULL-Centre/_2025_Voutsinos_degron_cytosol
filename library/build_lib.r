options(width=160)

# Find a golden-ish transcript
# This is neither the golden nor canonical transcript but in most cases it should be and in all cases it should be a relevant transcript
get_canonical_transcript = function(gene) {
    i = which(d$Gene.stable.ID == gene)
    if (length(i) > 1) {
        # APPRIS can be principal1 to 5, alternative1 or 2, or nothing
    	appris = d[i,"APPRIS.annotation"]
	# the 2 first also get the principal score
	# P1 should beat everything, 
        appris_score = 40 * (appris == "principal1") + 10 * (appris == "principal2") +
	               10 * (substr(appris, 1, 9) == "principal") + 1  * (substr(appris, 1, 9) == "alternati")

        # Transcript support level 1 to 3 where tsl3 is only observed once, tsl2 more, and tsl1 has a lot of observations
	# tsl1 beats tsl2+refseq which is not always the golden transcript, 
        tsl = substr(d[i,"Transcript.support.level..TSL."], 1, 4)
        tsl_score = 25 * (tsl=="tsl1") + 10 * (tsl=="tsl2") + 1 * (tsl=="tsl3")

        # Gencode basic means that the transcript is complete
        gencode_score = 10 * (d[i,"GENCODE.basic.annotation"] == "GENCODE basic")

        # Match to a RefSeq entry is also a good sign, chance being the MANE select transcript
        refseq_score = 10 * (d[i,"RefSeq.match.transcript"] != "")
	
	# Transcript length score
	# If everything else match, choose the longest
	length_score = rank(d[i,"Transcript.length..including.UTRs.and.CDS."], ties.method="average")
	length_score = 1 * (length_score/max(length_score))

        score = appris_score + tsl_score + gencode_score + refseq_score + length_score
	isort = order(score, decreasing=TRUE)
        i = i[isort]
	score = score[isort]
	if (score[1] == score[2]) {
	    print(sprintf("WARNING: Transcripts with same score for %s",gene))
	    print(cbind(d[i,c(1,3,5,6,7,8,10)],score))
	}
	i = i[1]
    }
    # print(sprintf("Using transcript %s for gene %s",d[i,"Transcript.stable.ID"],gene))
    d[i,"Transcript.stable.ID"]
}


tile_seq = function(sequence, tile_len, tile_stride, make_ct_tile=T) {
    n = nchar(sequence)
    if (n < tile_len) { return(NULL) }
    n_tiles = (n-tile_len) %/% tile_stride +1
    first = seq(0,n_tiles-1)*tile_stride+1
    last = first+tile_len-1
    # if whole tiles cannot cover the sequence, make a terminal tile if requested
    if (make_ct_tile & (n_tiles-1)*tile_stride+tile_len < n) {
        first = c(first, n-tile_len+1)
	last = c(last, n)
    }
    tiles = substring(sequence, first, last)
    
    name = names(sequence)
    if (is.null(name)) {
        names(tiles) = paste(seq_along(first), first, last, sep=".")
    } else {
        names(tiles) = paste(name, seq_along(first), first, last, sep=".")
    }
    return(tiles)
}

codon_table = list()
# T           T                            C                           A                           G
codon_table[["TTT"]] = "F"; codon_table[["TCT"]] = "S"; codon_table[["TAT"]] = "Y"; codon_table[["TGT"]] = "C" # T
codon_table[["TTC"]] = "F"; codon_table[["TCC"]] = "S"; codon_table[["TAC"]] = "Y"; codon_table[["TGC"]] = "C" # C
codon_table[["TTA"]] = "L"; codon_table[["TCA"]] = "S"; codon_table[["TAA"]] = "*"; codon_table[["TGA"]] = "*" # A
codon_table[["TTG"]] = "L"; codon_table[["TCG"]] = "S"; codon_table[["TAG"]] = "*"; codon_table[["TGG"]] = "W" # G
# C
codon_table[["CTT"]] = "L"; codon_table[["CCT"]] = "P"; codon_table[["CAT"]] = "H"; codon_table[["CGT"]] = "R" # T
codon_table[["CTC"]] = "L"; codon_table[["CCC"]] = "P"; codon_table[["CAC"]] = "H"; codon_table[["CGC"]] = "R" # C
codon_table[["CTA"]] = "L"; codon_table[["CCA"]] = "P"; codon_table[["CAA"]] = "Q"; codon_table[["CGA"]] = "R" # A
codon_table[["CTG"]] = "L"; codon_table[["CCG"]] = "P"; codon_table[["CAG"]] = "Q"; codon_table[["CGG"]] = "R" # G
# A
codon_table[["ATT"]] = "I"; codon_table[["ACT"]] = "T"; codon_table[["AAT"]] = "N"; codon_table[["AGT"]] = "S" # A
codon_table[["ATC"]] = "I"; codon_table[["ACC"]] = "T"; codon_table[["AAC"]] = "N"; codon_table[["AGC"]] = "S" # C
codon_table[["ATA"]] = "I"; codon_table[["ACA"]] = "T"; codon_table[["AAA"]] = "K"; codon_table[["AGA"]] = "R" # A
codon_table[["ATG"]] = "M"; codon_table[["ACG"]] = "T"; codon_table[["AAG"]] = "K"; codon_table[["AGG"]] = "R" # G
# G
codon_table[["GTT"]] = "V"; codon_table[["GCT"]] = "A"; codon_table[["GAT"]] = "D"; codon_table[["GGT"]] = "G" # A
codon_table[["GTC"]] = "V"; codon_table[["GCC"]] = "A"; codon_table[["GAC"]] = "D"; codon_table[["GGC"]] = "G" # C
codon_table[["GTA"]] = "V"; codon_table[["GCA"]] = "A"; codon_table[["GAA"]] = "E"; codon_table[["GGA"]] = "G" # A
codon_table[["GTG"]] = "V"; codon_table[["GCG"]] = "A"; codon_table[["GAG"]] = "E"; codon_table[["GGG"]] = "G" # G

translate = function(dna, allow_truncation=FALSE, non_nat_char="x") {
    n = nchar(dna)
    codons = substring(dna, seq(1, n-2, by=3), seq(3, n, by=3))
    if (length(codons)*3 != n & ! allow_truncation) return(NA)
    aas = ifelse(codons %in% names(codon_table), codon_table[codons], non_nat_char)
    paste0(aas, collapse="")
}

replace_restriction_site = function(query_seq, target_seq, replace_frame0, replace_frame1, replace_frame2) {
    # target_seq            GCGGCCGC
    #                 ___   ___   __
    # replace_frame0  ATG...GCaGCCGC : Ala GCG -> GCA
    #                  ___   ___   _
    # replace_frame1   ATG..GCGcCCGC : Arg CGG -> CGC
    #                   ___   ___   
    # replace_frame2    ATG.GCGGaCGC : Gly GGC -> GGA
    query_seq = toupper(query_seq)
    target_seq = toupper(target_seq)
    # Make a replace list where the modulus of 3 may be used as index for fast replacement
    replace_list = c(toupper(replace_frame0), toupper(replace_frame2), toupper(replace_frame1))
    # use perl regular expressions with look-ahead assertions (see ?regex) to find overlapping matches
    # e.g. gregexpr("(?=CC)", "ACCCAAGACCA")
    grep_list = gregexpr(paste0("(?=",target_seq,")"), query_seq, perl=TRUE)
    new_query_seq = rep("", length(query_seq))
    ns = 0; nr = 0
    for (i in seq_along(query_seq)) {
        match_vec = grep_list[[i]]
        s = query_seq[i]
        if (match_vec[1] > 0) {
	    for (j in seq_along(match_vec)) {
                # An overlapping site may already be changed by prev. match so test if it is still there
	        if (substr(s, match_vec[j], match_vec[j]+nchar(target_seq)-1) == target_seq) {
		    # Replace
                    s = paste0(substr(s, 1, match_vec[j]-1),
                               replace_list[(match_vec[j]-1) %% 3 +1],
	                       substr(query_seq[i], match_vec[j]+nchar(target_seq), nchar(query_seq[i])) )
                    nr = nr+1
		}
            }
	    # Test if amino acid sequence is the same
	    # stopifnot(dna2aa(query_seq[i]) == dna2aa(toupper(s)))
	    stopifnot(translate(query_seq[i]) == translate(toupper(s)))
	    # print(paste(translate(query_seq[i]),query_seq[i],i,collapse=", "))
	    ns = ns+1
        }
        new_query_seq[i] = s
    }
    print(sprintf("Replaced %d restriction sites in %d of %d sequences",nr,ns,length(query_seq)))
    return(new_query_seq)
}


# load genes of cytosolic proteins and find a transcript for each
d = read.csv(gzfile("biomart_cytosol.csv.gz"))
genes = unique(d$Gene.stable.ID)
transcripts = sapply(genes, get_canonical_transcript)

# load sequences of all transcripts
s = read.table(gzfile("biomart_cytosol_seq.txt.gz"))
colnames(s) = c("name","seq")
name_split_list = strsplit(s$name, "|", fixed=TRUE)
uniprot = sapply(name_split_list, function(v){if (length(v) >4) return(v[5]) else "NoGene"})
seq_df = data.frame(i=seq_along(uniprot), gene=sapply(name_split_list, '[[', 1), transcript=sapply(name_split_list, '[[', 3), uniprot=uniprot, seq=s$seq)

# make a dataframe of cytosolic proteins using selected transcripts
i_t2c = match(transcripts,seq_df$transcript)
cytosol = data.frame(gene=genes, transcript=transcripts, uniprot=seq_df[i_t2c,"uniprot"], seq=seq_df[i_t2c,"seq"])

# truncate sequence that have stop codons
cytosol$seq = substr(cytosol$seq, 1, nchar(cytosol$seq)-3)
cytosol$n_base = nchar(cytosol$seq)

# select sequences longer than 90 nucleotides and remove stop codon
i_rm = which(cytosol$n_base < 90)
print(sprintf("Remove %d sequence(s) that are shorter than 90 bases (a tile): %s", length(i_rm), paste(cytosol[i_rm,"gene"],collapse=",")))
cytosol = cytosol[-i_rm,]

# check that genes and transcripts are unique
stopifnot(length(cytosol$transcript) == length(unique(cytosol$transcript)))
stopifnot(length(cytosol$gene) == length(unique(cytosol$gene)))
print(sprintf("Selected %d unique genes and transcripts of cytosolic proteins", nrow(cytosol)))

# make a new vector of sequences that can have names to be returned by 'lapply' and expanded by 'unlist'
cyt_seq = cytosol$seq
names(cyt_seq) = cytosol$transcript
tiles = unlist(lapply(cyt_seq, tile_seq, 90, 45))

n_s_unq = length(unique(cyt_seq))
n_t_unq = length(unique(tiles))
print(sprintf("Sliced %d sequences (%d or %.0f%% unique) into %d tiles (%d or %.0f%% unique)",
              length(cyt_seq), n_s_unq, n_s_unq/length(cyt_seq)*100,
	      length(tiles),   n_t_unq, n_t_unq/length(tiles)*100))

# check that these are identical to original set
old_tiles = read.table("../../library/biomart_cytosol/cytosol90.seq")
stopifnot(all(tiles == old_tiles$V2))

name_split_list = strsplit(names(tiles), ".", fixed=T)
tile_lib = data.frame(meta = names(tiles),
                      name =NA,
                      transcript = sapply(name_split_list, '[[', 1),
                      gene = NA,
		      uniprot=NA,
                      tile_number = as.numeric( sapply(name_split_list, '[[', 2) ),
                      first_resi_dna = as.numeric( sapply(name_split_list, '[[', 3) ), 
                      last_resi_dna = as.numeric( sapply(name_split_list, '[[', 4) ),
		      dna = tiles)
i_c2t = match(tile_lib$transcript,cytosol$transcript)
tile_lib$gene = cytosol[i_c2t,"gene"]
tile_lib$uniprot = cytosol[i_c2t,"uniprot"]
tile_lib$name = sprintf("%s_tile%04d", tile_lib$gene, tile_lib$tile_number)
rownames(tile_lib) = NULL

# Make libraries of even, odd and CT fragments
print("Split library into 3")
tile_lib$lib = NA
tile_lib[which(tile_lib$tile_number %% 2 == 0),"lib"] = 1
tile_lib[which(tile_lib$tile_number %% 2 == 1),"lib"] = 2
tile_lib[which(diff(tile_lib$tile_number) != 1),"lib"] = 3
tile_lib[length(tile_lib$tile_number),"lib"] = 3  # last fragment is also CT
print(sprintf("Library size %d split into lib1_even (%d), lib2_odd (%d), and lib3_ct (%d)",
              length(tile_lib$lib), sum(tile_lib$lib==1), sum(tile_lib$lib==2), sum(tile_lib$lib==3)))
stopifnot(length(which(is.na(tile_lib$lib))) == 0)

find_nearest_tile_n = function(i,indices,n) {
    i_prev = i
    while (tile_lib[indices[i_prev],"tile_number"] != n & i_prev > 0)  i_prev = i_prev-1
    i_post = i
    while (tile_lib[indices[i_post],"tile_number"] != n & i_post < length(indices))  i_post = i_post+1
    print(tile_lib[indices[i_prev:i_post], c("name","tile_number")])
    if (i_post-i <= i-i_prev & i_post < length(indices)) {
        print(sprintf("Change split index %d (%d) to %d (%d) with tile number %d (%d)",
                      i,indices[i],i_post,indices[i_post],n,tile_lib[indices[i_post],"tile_number"]))
        i = i_post
    } else if (i_post-i > i-i_prev & i_prev > 0) {
        print(sprintf("Change split index %d (%d) to %d (%d) with tile number %d (%d)",
                      i,indices[i],i_prev,indices[i_prev],n,tile_lib[indices[i_prev],"tile_number"]))
        i = i_prev
    } else {
        print(sprintf("Keep index %d (%d) with tile number %d (target %d)",i,indices[i],tile_lib[indices[i],"tile_number"],n))
    }
    return(i)
}

# Make random sublibraries
print("Split libraries into 2 sub-libraries")
tile_lib$sublib = NA

# Split lib 1 in 2 by random
i_sub1 = which(tile_lib$lib == 1)
i_split1 = find_nearest_tile_n(length(i_sub1) %/% 2, i_sub1, 2)
# i_sub12 = sample(i_sub1, length(i_sub1) %/% 2)
i_sub12 = i_sub1[i_split1:length(i_sub1)]
tile_lib[i_sub1,"sublib"] = 1
tile_lib[i_sub12,"sublib"] = 2
print(sprintf("Lib1_even size %d split into sublib1 (%d) and sublib2 (%d)",
              length(i_sub1),length(i_sub1)-length(i_sub12),length(i_sub12)))

# Split lib 2 in 2 by random
i_sub2 = which(tile_lib$lib == 2)
i_split2 = find_nearest_tile_n(length(i_sub2) %/% 2, i_sub2, 1)
# i_sub24 = sample(i_sub2, length(i_sub2) %/% 2)
i_sub24 = i_sub2[i_split2:length(i_sub2)]
tile_lib[i_sub2,"sublib"] = 3
tile_lib[i_sub24,"sublib"] = 4
print(sprintf("Lib2_odd size %d split into sublib3 (%d) and sublib4 (%d)",
              length(i_sub2),length(i_sub2)-length(i_sub24),length(i_sub24)))

# Assign all of lib3 to the same sublib 
i_sub3 = which(tile_lib$lib == 3)
tile_lib[i_sub3,"sublib"] = 1
print("Lib3 assigned to sublib1")

stopifnot(length(which(is.na(tile_lib$sublib))) == 0)

dna = tile_lib$dna
print("Replace NotI restriction site GCGGCCGC")
dna = replace_restriction_site(dna, "GCGGCCGC", "GCaGCCGC", "GCGcCCGC", "GCGGaCGC")

print("Replace BsiWI restriction site CGTACG")
dna_rm_restrict = replace_restriction_site(dna, "CGTACG", "CGaACG", "CGTtCG", "CGTAtG")
tile_lib$dna_rm_restrict = dna_rm_restrict

# Adapters, nt is 5 prime and ct is 3 prime
print("Append adapters and sub-library adapters on all sequences")
adapt = list()
adapt[['all']]['nt'] = "GCGGCCGCGTTCTAGAGGCAGCGGAGCCACC"
adapt[['all']]['ct'] = "TAGTAACTTAAGAATTCACCGGTCTGACCTCGTACG"
adapt[['sub1']]['nt'] = "GTTCGTATGACTACGCTC"
adapt[['sub1']]['ct'] = "GCAACGAGATGAGTGATC"
adapt[['sub2']]['nt'] = "GGTCAATTCTAGTGAACG"
adapt[['sub2']]['ct'] = "GTCGTCCTGATACTTTGG"
adapt[['sub3']]['nt'] = "ATCCGTCATGTATGAGAG"
adapt[['sub3']]['ct'] = "TCATGTTTCCGTTGTAGG"
adapt[['sub4']]['nt'] = "GGTATTTGCTCGACAATC"
adapt[['sub4']]['ct'] = "AGTGAATAATCTTGGCGC"

# firstly, add common adapters to sequences with rescriction sites removed
dna = paste0(adapt[['all']]['nt'], tile_lib$dna_rm_restrict, adapt[['all']]['ct'])

# secondly, add sub-library specific adapters
#   build vectors with N-terminal (left or 5') and C-terminal (right or 3') sub-library adapters
adapt_nt = rep("", length(dna))
adapt_nt[which(tile_lib$sublib == 1)] = adapt[['sub1']]['nt']
adapt_nt[which(tile_lib$sublib == 2)] = adapt[['sub2']]['nt']
adapt_nt[which(tile_lib$sublib == 3)] = adapt[['sub3']]['nt']
adapt_nt[which(tile_lib$sublib == 4)] = adapt[['sub4']]['nt']

adapt_ct = rep("", length(dna))
adapt_ct[which(tile_lib$sublib == 1)] = adapt[['sub1']]['ct']
adapt_ct[which(tile_lib$sublib == 2)] = adapt[['sub2']]['ct']
adapt_ct[which(tile_lib$sublib == 3)] = adapt[['sub3']]['ct']
adapt_ct[which(tile_lib$sublib == 4)] = adapt[['sub4']]['ct']

# assemble and check lengths
stopifnot(length(adapt_nt) == length(dna))
stopifnot(length(adapt_ct) == length(dna))
tile_lib$dna_adapt = paste0(adapt_nt, dna, adapt_ct)
stopifnot(all(nchar(tile_lib$dna_adapt) == (90+nchar(adapt[['all']]['nt'])+nchar(adapt[['all']]['ct'])+2*18 )))

# Historically, the following code is added much later

# control sequences
ctrl = read.table("control_sequences.csv", header=T, sep=";")
colnames(ctrl) = c("meta","dna","last_resi_dna")
ctrl$first_resi_dna = 1
ctrl$pool = "control"
ctrl$sublib = "all"
ctrl$tile_number = seq(nrow(ctrl))
# add underscore to control names to sort them last
ctrl$name = sprintf("_%s",ctrl$meta)
i_underscore = regexpr("_", ctrl$meta)
ctrl$gene = paste0("_", substr(ctrl$meta, 1, i_underscore-1))
ctrl$dna_rm_restrict = ctrl$dna

# merge tile library with control sequences
tile_lib$pool = NA
tile_lib[which(tile_lib$lib==1),"pool"] = "even"
tile_lib[which(tile_lib$lib==2),"pool"] = "odd"
tile_lib[which(tile_lib$lib==3),"pool"] = "ct"
stopifnot(all(! is.na(tile_lib$pool)))

tile_lib$sublib_int = tile_lib$sublib
tile_lib$sublib = sprintf("sublib%d",tile_lib$sublib_int)
stopifnot(all(! is.na(tile_lib$sublib)))

cns = c("name","gene","tile_number","dna","dna_rm_restrict","pool","sublib")
lib = rbind(tile_lib[,cns], ctrl[,cns])
lib$aa = sapply(lib$dna, translate)
stopifnot( sum(is.na(lib$aa)) == 0 )

# this should be the order
resi_first_dna = c(tile_lib$first_resi_dna, ctrl$first_resi_dna)
resi_last_dna = c(tile_lib$last_resi_dna, ctrl$last_resi_dna)
# all should have an integer number of codons
stopifnot(all((resi_last_dna-resi_first_dna+1) %% 3 == 0))
lib$resi_first = (resi_first_dna-1) %/% 3 +1
lib$resi_last  = (resi_last_dna-1)  %/% 3 +1

# check DNA indices
out_of_frame = unique(c(which((resi_first_dna-1) %% 3 != 0),
                        which((resi_last_dna) %% 3 != 0)))
if (length(out_of_frame) > 0) {
    print(sprintf("There are %d tile DNA indices that seems out of frame",length(out_of_frame)))
    print(cbind(lib[out_of_frame,c("name","pool","sublib","aa")],
                data.frame(resi_first_dna=resi_first_dna[out_of_frame],resi_last_dna=resi_last_dna[out_of_frame])))
}

# order library
lib = lib[order(lib$name),]
rownames(lib) = NULL

# Assign gene copy to genes with identical DNA sequences
cytosol$dna_copyof = NA
t = table(cytosol$seq)
redundant_dna = names(t[t>1])
for (dna in redundant_dna) {
    # find primary copy by alphabetic ordring
    i_dna = which(cytosol$seq == dna)
    stopifnot(length(i_dna) > 1)
    i_dna = i_dna[order(cytosol[i_dna,"gene"], cytosol[i_dna,"transcript"])]
    stopifnot( all(is.na(cytosol[i_dna,"dna_copyof"])) )
    cytosol[i_dna[2:length(i_dna)],"dna_copyof"] = i_dna[1]
}

# Assign unique uniprot to genes without uniprot map
i_unmapped = which(cytosol$uniprot == "NoGene")
ii_notcopy = which(is.na(cytosol[i_unmapped,"dna_copyof"]))
cytosol[i_unmapped[ii_notcopy],"uniprot"] = sprintf("NoUnip%02d",seq_along(ii_notcopy))

ii_copy = which(! is.na(cytosol[i_unmapped,"dna_copyof"]))
cytosol[i_unmapped[ii_copy],"uniprot"] = cytosol[ cytosol[i_unmapped[ii_copy],"dna_copyof"] , "uniprot" ]

stopifnot(all(cytosol$uniprot != "NoGene"))

# try to remove gene copies upfront
cytosol_redun = cytosol
cytosol = cytosol_redun[which(is.na(cytosol_redun$dna_copyof)),]
cytosol$dna_copyof = NULL


# amino acid sequences
cytosol$seq_aa = sapply(cytosol$seq, translate)

# DNA sequences with a non-integer number of codons
for (i in which(is.na(cytosol$seq_aa))) {
    excess_bases = cytosol[i,"n_base"] %% 3
    cytosol[i,"seq_aa"] = translate(substr(cytosol[i,"seq"], 1, cytosol[i,"n_base"] - excess_bases))
    stopifnot(! is.na(cytosol[i,"seq_aa"]))
    print(sprintf("WARNING: Removed %d excess bases from %s %s in order to translate to %d amino acids. CT tile may be out of frame",
          excess_bases, cytosol[i,"gene"], cytosol[i,"uniprot"], nchar(cytosol[i,"seq_aa"])))
    print(sprintf("CT tile in lib   %s  %s", lib[which(lib$gene==cytosol[i,"gene"] & lib$pool=="ct"),"dna"],
                  translate(lib[which(lib$gene==cytosol[i,"gene"] & lib$pool=="ct"),"dna"]) ))
    print(sprintf("CT tile inframe  %s  %s", substr(cytosol[i,"seq"], nchar(cytosol[i,"seq"])-89-excess_bases, nchar(cytosol[i,"seq"])-excess_bases),
                  substr(cytosol[i,"seq_aa"],nchar(cytosol[i,"seq_aa"])-29,nchar(cytosol[i,"seq_aa"]))))
}

# copies at aa-level - are there more copies here compared to DNA level?
t = table(cytosol$seq_aa)
print(sprintf("Of %d proteins without DNA copies, %d (%.2f%%) have copies at amino acid level",nrow(cytosol),sum(t>1),sum(t>1)/nrow(cytosol)*100))
# Example is HSP70: P0DMV8 and P0DMV9 with 8 synonymous copies

cytosol$aa_copyof = NA
t = table(cytosol$seq_aa)
redundant_aa = names(t[t>1])
for (aa in redundant_aa) {
    # find primary copy by alphabetic ordring - these can have DNA copies
    i_aa = which(cytosol$seq_aa == aa)
    # i_aa = which(cytosol$seq_aa == aa & is.na(cytosol$aa_copyof))
    # ordering ensures that a dna-main is selected as aa-main tile when selecting index 1
    i_aa = i_aa[order(cytosol[i_aa,"gene"], cytosol[i_aa,"transcript"])]
    stopifnot(length(i_aa) > 1)
    stopifnot( all(is.na(cytosol[i_aa,"aa_copyof"])) )
    cytosol[i_aa[2:length(i_aa)],"aa_copyof"] = i_aa[1]

    # # check that all but one protein is annotated as copies
    # i_aa = which(cytosol$seq_aa == aa)
    # stopifnot( sum(is.na(cytosol[i_aa,"aa_copyof"])) == 1 )
}

# does unique aa-sequences have unique uniprot entries
t = table(cytosol[which(is.na(cytosol$aa_copyof)),"uniprot"])
nonuniq_unip = names(t[t>1])
print("Possible isoforms and variants")
for (unip in nonuniq_unip) {
    i_unip = which(cytosol$uniprot==unip)
    i_uniq_aa = i_unip[which(is.na(cytosol[i_unip,"aa_copyof"]))]
    print(sprintf("    Uniprot %s occurs %d times with %d unique aa-sequences of lengths %s",
                  unip, length(i_unip), length(i_uniq_aa), paste(nchar(cytosol[i_uniq_aa,"seq_aa"]),collapse=",")))
}

# annotate each protein with a list of genes that are a copy of this gene on aa-level
i = which(! is.na(cytosol$aa_copyof))
agg = aggregate(cytosol[i,"gene"], by=list(cytosol[i,"aa_copyof"]), paste0, collapse=":")
cytosol$copies = NA
cytosol[agg[,1],"copies"] = agg[,2]

# check that all genes appears as either main (aa_copyof==NA) or in list of copies
genes = c(cytosol[which(is.na(cytosol$aa_copyof)),"gene"], na.omit(unlist(strsplit(cytosol$copies,":"))))
stopifnot( length(genes)==length(unique(genes)) )
stopifnot( setequal(cytosol$gene, genes) )


# === Clustering of full amino acid sequences
# dump sequences for MMseqs2 search for clustering (no controls in cytosol dataframe)
i_uniq = which(is.na(cytosol$aa_copyof))
# write.table(cytosol[i_uniq,c("gene","seq_aa")], file="cytosol_uniq.seq", row.names=F, col.names=F, quote=F)

# Please cite: M. Steinegger and J. Soding. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nature Biotechnology, doi:10.1038/nbt.3988 (2017).
# MMseqs2 Version: f5b588be54de35876eff90f3694b67e1d00856ef
# /projects/prism/people/skr526/mmseqs/bin/mmseqs createdb cytosol_uniq.fasta cytosol_uniq_db
# /projects/prism/people/skr526/mmseqs/bin/mmseqs cluster -s 7.5 -e 0.001 -c 0.6 --cluster-mode 1 cytosol_uniq_db cytosol_uniq_clust_c60_e-3 tmp
# /projects/prism/people/skr526/mmseqs/bin/mmseqs createtsv cytosol_uniq_db cytosol_uniq_db cytosol_uniq_clust_c60_e-3 cytosol_uniq_clust_c60_e-3.tsv

# read mmseqs clusters
cytosol_cluster = read.table("cytosol_uniq_clust_c60_e-3.tsv", sep="\t")
colnames(cytosol_cluster) = c("cluster_rep","cluster_member")
stopifnot( setequal(cytosol[i_uniq,"gene"], cytosol_cluster$cluster_member) )
print(sprintf("Read seq clustering of %d unique proteins in %d clusters",length(unique(cytosol_cluster$cluster_member)),length(unique(cytosol_cluster$cluster_rep))))

# merge into cytosol dataframe, amino acid sequences that are identical are also similar
cytosol$aa_similarto = cytosol$aa_copyof
# aditionally, assign cluster members with their respective cluster representatives
cytosol[match(cytosol_cluster$cluster_member,cytosol$gene),"aa_similarto"] = match(cytosol_cluster$cluster_rep, cytosol$gene)
cytosol[match(cytosol_cluster$cluster_member,cytosol$gene),"aa_similarto"] = cytosol_cluster$cluster_rep

# annotate clustered sequences with cluster size (incl. cluster representative)
cytosol$cluster_size = NA
clust_tab = table(cytosol_cluster$cluster_rep)
cytosol[match(names(clust_tab),cytosol$gene),"cluster_size"] = clust_tab

# cluster representatives with low similarity to everything else
nredun_genes = cytosol[which(cytosol$aa_similarto == cytosol$gene),"gene"]
solo_genes = cytosol[which(cytosol$cluster_size == 1),"gene"]
print(sprintf("Library is made from %d (%.2f%%) AA unique genes, %d (%.2f%%) non-redundant and %d (%.2f%%) with cluster size 1",
              length(i_uniq), length(i_uniq)/nrow(cytosol)*100,
	      length(nredun_genes), length(nredun_genes)/nrow(cytosol)*100,
	      length(solo_genes), length(solo_genes)/nrow(cytosol)*100))
print("Distribution of cluster sizes")
print(table(clust_tab))


# === mark redundant tiles at DNA level with a primary tile ===
# This is used to make a set of non-redundant tile DNA sequences used for illumina read alignment
# Start from scratch with tile sequences because tile DNA sequences have restriction sites removed, DNA sequences in cytosol data frame does not.
print("Identfy identical tiles on DNA level")
lib$dna_copyof = NA

t = table(lib$dna)
redun_tile_dna = names(t[t>1])
print(sprintf("Found %d (%.2f%%) DNA-redundant tiles (%d unique)", sum(t[t>1]), sum(t[t>1])/nrow(lib)*100, length(redun_tile_dna)))

primary_genes = c()
# primary_genes = cytosol[which(is.na(cytosol$dna_copyof)),"gene"]
for (dna in redun_tile_dna) {
    i_dna = which(dna == lib$dna)
    # i_dna = i_dna[order(lib[i_dna,"gene"], lib[i_dna,"tile"])]
    i_dna = i_dna[order(lib[i_dna,"name"])]
    stopifnot(is.na(lib[i_dna,"dna_copyof"]))
    i_primary = i_dna[ which(lib[i_dna,"gene"] %in% primary_genes) ]
    if (length(i_primary) < 1) {
        i_primary = i_dna[1]
    } else {
        i_primary = i_primary[1]
    }
    i_copies = setdiff(i_dna, i_primary)
    lib[i_copies,"dna_copyof"] = i_primary
    primary_genes = c(primary_genes, lib[i_primary,"gene"])
}

# number of complete genes
i_nr = which(is.na(lib$dna_copyof))
# i_nr = i_nr[order(lib[i_nr,"gene"], lib[i_nr,"tile"])]
stopifnot(length(unique(lib[i_nr,"dna"])) == length(i_nr))
# tiles_per_gene = aggregate(lib$tile, by=list(lib$gene), length)
tiles_per_gene = aggregate(lib$name, by=list(lib$gene), length)
nr_tiles_per_gene = aggregate(lib[i_nr,"name"], by=list(lib[i_nr,"gene"]), length)

df = data.frame(gene=nr_tiles_per_gene[,1], nr_tiles=nr_tiles_per_gene[,2], tot_tiles=tiles_per_gene[match(nr_tiles_per_gene[,1], tiles_per_gene[,1]),2])
df = df[order(df$tot_tiles-df$nr_tiles),]

t = table(df$tot_tiles-df$nr_tiles)
print(sprintf("Total genes in non-redundant library: %d of %d (%.2f%%). %d complete (%.2f%%)",
              nrow(df), nrow(tiles_per_gene), nrow(df)/nrow(tiles_per_gene)*100, t[1], t[1]/nrow(tiles_per_gene)*100))


# mark redundant tiles at AA-level with a primary tile
print("Identfy identical tiles on AA level")

lib$aa = sapply(lib$dna, translate)

# are there redundant aa-sequences among DNA-unique tiles
i_uniq_dna = which(is.na(lib$dna_copyof))
t = table(lib[i_uniq_dna,"aa"])
print(sprintf("Of %d unique DNA tiles, %d (%.2f%%) have redundant AA sequence (%d unique)", length(i_uniq_dna), sum(t[t>1]), sum(t[t>1])/length(i_uniq_dna)*100, sum(t>1)))


# if a protein is a copy, then tiles should also be, so start with proteins
copied_genes = cytosol[which(! is.na(cytosol$aa_copyof)),"gene"]
# which tiles belongs to these proteins
il_cpgenes = which(lib$gene %in% copied_genes)
# cytosol index of copies for all tiles
ic_cp = cytosol[match(lib[il_cpgenes,"gene"],cytosol$gene),"aa_copyof"]
# tiles should have the same primary/main gene/protein as in the cytosol data frame
lib$aa_genecopyof = NA
lib[il_cpgenes,"aa_genecopyof"] = cytosol[ic_cp, "gene"]

lib$aa_copyof = match(paste(lib$aa_genecopyof, lib$tile, sep="_"), lib$name)
lib[which(is.na(lib$aa_genecopyof)),"aa_copyof"] = NA

# then, annotate a primary tile to the rest
t = table(lib[which(is.na(lib$aa_copyof)),"aa"])
redun_tile_aa = names(t[t>1])
print(sprintf("From %d copied proteins, %d copied tiles are annotated with same primary gene and %d additional redundant tiles (%d unique) with copied protein if possible",
              length(copied_genes), length(il_cpgenes), sum(t[t>1]), sum(t>1)))
for (aa in redun_tile_aa) {
    i_aa = which(aa == lib$aa)
    # some of these may already be annotated as copies of another tile
    primary_genes = na.omit(unique(lib[i_aa,"aa_genecopyof"]))
    if (length(primary_genes) < 1) {
        # None of the tiles are annotated so order according to gene and choose the first as primary
        # i_aa = i_aa[order(lib[i_aa,"gene"], lib[i_aa,"tile"])]
        i_aa = i_aa[order(lib[i_aa,"name"])]
        i_primary = i_aa[1]
    } else if (length(is.na(lib[i_aa,"aa_copyof"])) == 1) {
        # if just one is not annotated as being a copy, use this as primary
        i_primary = i_aa[which(is.na(lib[i_aa,"aa_copyof"]))]
        print(sprintf("Tile %s with %d occur. already has primary tile %s",
	              aa, length(i_aa), lib[i_primary,"name"]))
    } else {
        # select the primary tile from the gene that most already annotated copies point to
        if (length(primary_genes) > 1) {
	    # for primary, select the gene that has the most tiles already point to
	    t = table(na.omit(lib[i_aa,"aa_genecopyof"]))
	    primary_genes = primary_genes[order(t)]
	    # print(sprintf("Tile %s with %d occur. assigned primary gene %s with most tiles already pointing to this: %s",
	    #               aa, length(i_aa), primary_genes[1], paste0(names(t),": ",t, collapse=", ")))
	}
        i_primary = i_aa[which(lib[i_aa,"gene"] == primary_genes[1])]
    }
    stopifnot(! is.na(i_primary))
    stopifnot(is.na(lib[i_primary,"aa_copyof"]))
    i_copies = setdiff(i_aa, i_primary)
    lib[i_copies,"aa_copyof"] = i_primary
}

# tiles should have the same primary/main gene/protein as in the cytosol data frame
# lib index of copies for all tiles
il_cp = lib[il_cpgenes,"aa_copyof"]

# cytosol index of copies for all tiles
ic_cp = cytosol[match(lib[il_cpgenes,"gene"],cytosol$gene),"aa_copyof"]
print(sprintf("Number of tiles that has a aa_copy different from the aa_copy of the gene is %d",
              sum(lib[il_cp, "gene"] != cytosol[ic_cp, "gene"])))

# check that all aa-redundant tiles are annotated
t = table(lib$aa)
redun_tile_aa = names(t[t>1])
print(sprintf("Checking %d aa-redundant tiles (%d unique)", sum(t[t>1]), sum(t>1)))
stopifnot(all( sapply(redun_tile_aa, function(aa){ sum(is.na( lib[which(lib$aa==aa),"aa_copyof"] )) }) == 1 ))

# check residue numbers by reconstructing tile sequences
il_tile = which(grepl("tile",lib$tile))
lib_aa_check = cytosol[match(lib$gene,cytosol$gene),"seq_aa"]
lib_aa_check[il_tile] = substr(lib_aa_check[il_tile], lib[il_tile,"resi_first"], lib[il_tile,"resi_last"])
il_check = il_tile[ which(lib[il_tile,"aa"] == lib_aa_check[il_tile]) ]
il_bad = setdiff(il_tile, il_check)
lib[il_bad,c("resi_first","resi_last")] = NA
print(sprintf("Of %d cytosolic library tiles (excl. controls), %d sequences does not check (%.1f%%)",
              length(il_tile), length(il_bad), length(il_bad)/length(il_tile)*100))
print(head(lib[il_bad,]))


# make unique sub-libraries

subpool_indices = function(pool, sublib, add_controls=T) {
    i_subpool = which(lib$pool==pool & lib$sublib==sublib)
    if (add_controls) { i_subpool = c(i_subpool, which(lib$pool=="control")) }
    
    # tiles in the subpool that are copies
    i_subpool_copy = i_subpool[which(! is.na(lib[i_subpool,"dna_copyof"]))]

    # for tiles that are copies, add the primary tile (if not already present)
    i_primary = lib[i_subpool_copy,"dna_copyof"]
    i_subpool_nr = union(setdiff(i_subpool, i_subpool_copy), i_primary)
    # i_subpool_nr = i_subpool_nr[order(lib[i_subpool_nr,"gene"],lib[i_subpool_nr,"tile"])]
    i_subpool_nr = i_subpool_nr[order(lib[i_subpool_nr,"name"])]

    # check that all tiles are present and that all DNA sequencea are unique in the library
    stopifnot( setequal(lib[i_subpool,"dna"], lib[i_subpool_nr,"dna"]) )
    stopifnot( length(unique(lib[i_subpool_nr,"dna"])) == length(i_subpool_nr) )
    # check that all sequences are primary
    stopifnot( all(is.na(lib[i_subpool_nr,"dna_copyof"])) )

    # report
    n_discard = length(i_subpool) - length(i_subpool_nr)
    print(sprintf("Pool %s %s has %d of %d (%.2f%%) redundant tiles",
                  pool, sublib, n_discard, length(i_subpool), n_discard/length(i_subpool)*100))

    # look at copies of anything in the selected sub library
    i_copies_offlib = union(i_primary, which(lib$dna_copyof %in% i_subpool_nr))
    other_pools = setdiff(unique(lib[i_copies_offlib,"pool"]), pool)
    for (op in other_pools) { print(sprintf("    Found %d tiles also in pool %s", sum(lib[i_copies_offlib,"pool"]==op), op)) }
    other_sublibs = setdiff(unique(lib[i_copies_offlib,"sublib"]), sublib)
    for (os in other_sublibs) { print(sprintf("    Found %d tiles also in %s", sum(lib[i_copies_offlib,"sublib"]==os), os)) }
    
    return(i_subpool_nr)
}

print("")
i_even1 = subpool_indices("even", "sublib1")
print("")
i_even2 = subpool_indices("even", "sublib2")
print("")
i_odd3  = subpool_indices("odd", "sublib3")
print("")
i_odd4  = subpool_indices("odd", "sublib4")
print("")
i_ct1   = subpool_indices("ct", "sublib1")
print("")
	
lib$even1 = FALSE
lib[i_even1,"even1"] = TRUE

lib$even2 = FALSE
lib[i_even2,"even2"] = TRUE

lib$odd3 = FALSE
lib[i_odd3,"odd3"] = TRUE

lib$odd4 = FALSE
lib[i_odd4,"odd4"] = TRUE

lib$ct1 = FALSE
lib[i_ct1,"ct1"] = TRUE

# number of sub-libraries that each tile is found in
nlib = apply(lib[,c("even1","even2","odd3","odd4","ct1")], MARGIN=1, sum)
print(sprintf("Summary of sub-library population of %d tiles (%d controls are in all 5 sublibraries; zero is copied tiles)", nrow(lib), nrow(ctrl)))
print(table(nlib))

# all tiles that are not in a librarty should be DNA copies
i = which(nlib==0)
stopifnot( sum(is.na(lib[i,"dna_copyof"])) == 0 )


# === Split data into test and train
# for hold-out/test data, select random 10% of sequences that are different from others
solo_genes = cytosol[which(cytosol$cluster_size == 1),"gene"]
i_uniq = which(is.na(cytosol$aa_copyof))
# holdout_genes = sample(solo_genes, length(i_uniq) %/% 10)
# # use a static selection
# write.table(holdout_genes, row.names=F, quote=F, col.names=F, file="holdout_selection.txt")
holdout_genes = read.table("holdout_selection.txt")[,1]
train_genes = cytosol[which(! cytosol$gene %in% holdout_genes & is.na(cytosol$aa_copyof)),"gene"]
stopifnot(length(holdout_genes)+length(train_genes)==length(i_uniq))

# Checks at tile (training data) level
il_holdout = which(lib$gene %in% holdout_genes)
il_use = which(lib$pool != "control" & is.na(lib$aa_copyof))
print(sprintf("Selected %d of %d (%.2f%%) used tiles from %d of %d (%.2f%%) solo genes (cluster size 1), or %.2f%% of AA unique genes",
              length(il_holdout), length(il_use), length(il_holdout)/length(il_use)*100,
              length(holdout_genes), length(solo_genes), length(holdout_genes)/length(solo_genes)*100, length(holdout_genes)/length(i_uniq)*100))


# If a holdout protein has a tile which is an amino acid copy of a training protein, move that protein from holdout to training
print("Redundancy filter of exact tile copies")
il_redun = il_holdout[which(! is.na(lib[il_holdout,"aa_copyof"]))]
redun_genes = unique(lib[il_redun,"gene"])
stopifnot(all(redun_genes %in% holdout_genes))
print(sprintf("    Of %d holdout tiles, %d (%.2f%%) tiles from %d genes (%.2f%% of hold-out) are exact copies of a training tile",
              length(il_holdout), length(il_redun), length(il_redun)/length(il_holdout)*100,
	      length(redun_genes), length(redun_genes)/length(holdout_genes)*100))

holdout_genes = setdiff(holdout_genes, redun_genes)
stopifnot(! any(redun_genes %in% holdout_genes))
stopifnot(! any(redun_genes %in% train_genes))
train_genes = c(train_genes, redun_genes)
stopifnot(length(holdout_genes)+length(train_genes)==length(i_uniq))
print(sprintf("    Moved %d proteins from hold-out set to training set", length(redun_genes)))


il_holdout = which(lib$gene %in% holdout_genes)
write.table(lib[il_holdout,c("name","aa")], file="holdout_tiles.seq", row.names=F, col.names=F, quote=F)
il_train = which(lib$gene %in% train_genes) # has copies
write.table(lib[il_train,c("name","aa")], file="train_tiles.seq", row.names=F, col.names=F, quote=F)
# original run from Jun 20 2023 found 621 tiles from 89 genes using -c 0.5 and rerun from 22 Aug 2024 579 tiles from 95 genes using -c 0.9 and many more with -c 0.5 (needed with new tile names)
# /projects/prism/people/skr526/mmseqs/bin/mmseqs easy-search -e 0.001 -c 0.9 holdout_tiles.fasta train_tiles.fasta leek_tiles.txt tmp


print("Redundancy filter of similar tiles")
leek_tiles = read.table("leek_tiles.txt")
redun_tiles = unique(leek_tiles$V1)
splitlist = strsplit(redun_tiles,"_")
redun_genes = unique(sapply(splitlist, "[", 1))
stopifnot(length(redun_genes)==length(unique(redun_genes)))
stopifnot(all(redun_genes %in% holdout_genes))
print(sprintf("    Of %d holdout tiles, %d (%.2f%%) tiles from %d genes (%.2f%% of hold-out) are similar to a training tile",
              length(il_holdout), length(redun_tiles), length(redun_tiles)/length(il_holdout)*100,
              length(redun_genes), length(redun_genes)/length(holdout_genes)*100))

holdout_genes = setdiff(holdout_genes, redun_genes)
stopifnot(! any(redun_genes %in% holdout_genes))
stopifnot(! any(redun_genes %in% train_genes))
train_genes = c(train_genes, redun_genes)
stopifnot(length(holdout_genes)+length(train_genes)==length(i_uniq))
print(sprintf("    Moved %d proteins from hold-out set to training set", length(redun_genes)))
stopifnot( setequal(cytosol[which(is.na(cytosol$aa_copyof)),"gene"], c(holdout_genes,train_genes)) )

il_holdout = which(lib$gene %in% holdout_genes)
# since these are cluster centers, there's a good chance that hold-out tiles wil have unique amino acid sequences
stopifnot(all(is.na(lib[il_holdout,"aa_copyof"])))
il_train = which(lib$gene %in% train_genes & is.na(lib$aa_copyof))
print(sprintf("Selected %d (%.2f%%) tiles from %d (%.2f%%) proteins", length(il_holdout), length(il_holdout)/sum(lib$pool!="contole")*100,
              length(holdout_genes), length(holdout_genes)/length(i_uniq)*100))

# assign each tile to a data set for model training
lib$set = NA
lib[il_holdout,"set"] = "test"
lib[il_train,"set"] = "train"
il_ctrl = which(lib$pool=="control")
stopifnot(all(is.na( lib[il_ctrl,"set"] )))
lib[il_ctrl,"set"] = "control"
il_cp = which(! is.na(lib$aa_copy) & lib$pool!="control")
stopifnot(all(is.na( lib[il_cp,"set"] )))
lib[il_cp,"set"] = "copy"
stopifnot(all( ! is.na(lib$set) ))
# few tile have stops which is not good for training
il_nons = which(grepl("*", lib$aa, fixed=T))
lib[il_nons,"set"] = paste(lib[il_nons,"set"], "nonsense", sep="-")  # not tested!
print("Tile sets")
print(table(lib$set))

cytosol$set = NA
cytosol[which(cytosol$gene %in% holdout_genes),"set"] = "test"
cytosol[which(cytosol$gene %in% train_genes),"set"] = "train"
ic_cp = which(! is.na(cytosol$aa_copyof))
stopifnot(all(is.na(cytosol[ic_cp,"set"])))
cytosol[ic_cp,"set"] = "copy"
stopifnot(all( ! is.na(cytosol$set) ))

print("Protein sets")
print(table(cytosol$set))

# # mapped uniprot entries
# unip_map_fn = "unique_proteins_240809_uniprot_mapped.csv"
# unip = read.csv(unip_map_fn)
# stopifnot(length(unip$seq_aa)==length(unique(unip$seq_aa)))
# # from the list of uniprot entries, find the first with the fewest characters
# unip$first_short = sapply(unip$uniprot, function(s){ l = strsplit(gsub("[]'[]","",s), ", ")[[1]]; l[order(nchar(l))[1]] })
# cytosol_unip_map = unip[match(cytosol$seq_aa,unip$seq_aa),"uniprot"]
# cytosol_first_short = unip[match(cytosol$seq_aa,unip$seq_aa),"first_short"]
# # if the original uniprot is in the list, use this. Else use the first&shortest unless this is nan, then use the original.
# cytosol$uniprot2 = apply(cbind(cytosol$uniprot,cytosol_unip_map,cytosol_first_short), MARGIN=1,
#                          function(v){ if(grepl(v[1],v[2])) v[1] else if(v[3]=="nan") v[1] else v[3] })
# n_diff = sum(cytosol$uniprot!=cytosol$uniprot2)
# print(sprintf("Made uniprot2 annotation (differ for %d proteins (%.1f%%)) with %d 'NoUnip' and %d 'nan' using file: %s", n_diff, n_diff/nrow(cytosol)*100,
#               sum(grepl("NoUnip", cytosol$uniprot2)), sum(grepl("nan", cytosol$uniprot2)), unip_map_fn))


# === dump ===

# # do not write out DNA sequences - these do not have restriction sites removed
# i_uniq = which(is.na(cytosol$aa_copyof))
# i_uniq = i_uniq[order(nchar(cytosol[i_uniq,"seq_aa"]), cytosol[i_uniq,"seq_aa"])]
# write.csv(cytosol[i_uniq,c("gene","transcript","seq_aa","copies","aa_similarto")], file="unique_proteins.csv", row.names=F, quote=F)

write.csv(lib, file="cytosol_tiles.csv", row.names=F, quote=F)
write.csv(cytosol, file="cytosol_genes.csv", row.names=F, quote=F)

# Dump R file
print("Dump R file cytosol.rda")
save(cytosol, tile_lib, lib, ctrl, adapt, file="cytosol.rda")

