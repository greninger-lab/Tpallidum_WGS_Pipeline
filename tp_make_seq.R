# Adapted for Nextflow Aug 2020 for T. pallidum

# HSV script but works with any viral sequence: This script makes a new reference sequence from de novo assembled scaffolds
# Pavitra Roychoudhury
# Sep 2017

# Built to be called from hsv_wgs_pipeline.sh with input arguments specifying input filename
# Requires wgs_functions.R which contains several utility scripts plus multiple R packages listed below

rm(list = ls())
sessionInfo()
library(Rsamtools)
library(GenomicAlignments)
library(Biostrings)
library(RCurl)
library(parallel)
# Get args from command line
args <- (commandArgs(TRUE))
if (length(args) == 0) {
    print("No arguments supplied.")
} else {
    sampname <- args[[1]]
    ref <- args[[2]]
    # for(i in 1:length(args)){)
    # 	# eval(parse(text=args[[i]]))
    # 	# print(args[[i]])
    # }
}

# For testing
# sampname='GH120016_CGATCCAC-TCGCGCAT_L002'
# ref='NC_016842'

# First import scaffolds and filter by length (>200) and coverage (>10x)
# scaffname=paste0('./contigs/',sampname,'/scaffolds.fasta')
scaffname <- paste0("assembly.fasta")
# with unicycler
contigs <- readDNAStringSet(scaffname, format = "fasta")
contigs <- contigs[width(contigs) > 200]
# cov<-unlist(lapply(names(contigs),function(x)
#   as.numeric(strsplit(x,'_cov_')[[1]][2])))
# contigs<-contigs[cov>10]
scaffname_filtered <- gsub("assembly.fasta", "scaffolds_filtered.fasta", scaffname)
writeXStringSet(contigs, scaffname_filtered)

reffname <- paste0("", ref, ".fasta")
# Call bwa mem to align contigs against reference
# samfname<-gsub('scaffolds.fasta',
#                paste0(sampname,'_aligned_scaffolds_',strsplit(basename(reffname),'.fasta')[[1]][1],'.sam'),
#                scaffname)
samfname <- gsub(
    "assembly.fasta",
    paste0(sampname, "_aligned_scaffolds_", strsplit(basename(reffname), ".fasta")[[1]][1], ".sam"),
    scaffname
)
samfname <- gsub("./", "", samfname)
print(paste("Aligning contigs to reference", basename(reffname), "..."))
system(paste("/usr/local/miniconda/bin/bwa mem", reffname, scaffname_filtered, ">", samfname))
# Compress and clean up
print("Converting sam to bam ...")
bamfname <- gsub(".sam", ".bam", samfname)
system(paste(
    "/usr/src/samtools-1.9/samtools view -bh -o",
    bamfname, samfname,
    "-T", reffname
))
print("done ...")
# file.remove(samfname);
# rm(samfname)
print("Sort bam ...")
system(paste(
    "/usr/src/samtools-1.9/samtools sort -o",
    gsub(".bam", ".sorted.bam", bamfname),
    bamfname
))
# file.remove(bamfname);
# file.rename(gsub('.bam','.sorted.bam',bamfname),bamfname)
print("done ...")

bamfname <- gsub(".bam", ".sorted.bam", bamfname)
print("bamfname")
print(bamfname)


# Make a new reference from scaffolds
make_ref_from_assembly <- function(bamfname, reffname) {
    require(Rsamtools)
    require(GenomicAlignments)
    require(parallel)
    ncores <- detectCores()

    # Read reference sequence
    ref_seq <- readDNAStringSet(reffname)

    if (!is.na(bamfname) & class(try(scanBamHeader(bamfname), silent = T)) != "try-error") {
        # Index bam if required
        if (!file.exists(paste(bamfname, ".bai", sep = ""))) {
            baifname <- indexBam(bamfname)
        } else {
            baifname <- paste(bamfname, ".bai", sep = "")
        }

        # Import bam file
        params <- ScanBamParam(
            flag = scanBamFlag(isUnmappedQuery = FALSE),
            what = c("qname", "rname", "strand", "pos", "qwidth", "mapq", "cigar", "seq")
        )
        gal <- readGAlignments(bamfname, index = baifname, param = params)

        # Remove any contigs with width <200 bp
        gal <- gal[width(gal) > 200]

        # First lay contigs on reference space--this removes insertions and produces a seq of the same length as ref
        qseq_on_ref <- sequenceLayer(mcols(gal)$seq, cigar(gal), from = "query", to = "reference")
        qseq_on_ref_aligned <- stackStrings(qseq_on_ref, 1, max(mcols(gal)$pos + qwidth(gal) - 1, width(ref_seq)),
            shift = mcols(gal)$pos - 1, Lpadding.letter = "N", Rpadding.letter = "N"
        )

        # Make a consensus matrix and get a consensus sequence from the aligned scaffolds
        cm <- consensusMatrix(qseq_on_ref_aligned, as.prob = T, shift = 0)[c("A", "C", "G", "T", "N", "-"), ]
        # cm[c('N','-'),]<-0;
        cm["N", ] <- 0
        cm <- apply(cm, 2, function(x) if (all(x == 0)) {
            return(x)
        } else {
            return(x / sum(x))
        })
        cm["N", colSums(cm) == 0] <- 1
        con_seq <- DNAStringSet(gsub("\\?", "N", consensusString(cm, threshold = 0.25)))
        con_seq <- DNAStringSet(gsub("\\+", "N", con_seq))


        # Now fill in the Ns with the reference
        temp <- as.matrix(con_seq)
        temp[temp == "N"] <- as.matrix(ref_seq)[temp == "N"]
        con_seq <- DNAStringSet(paste0(temp, collapse = ""))
        names(con_seq) <- sub(".bam", "_consensus", basename(bamfname))

        # Look for insertions in bam cigar string
        cigs_ref <- cigarRangesAlongReferenceSpace(cigar(gal),
            with.ops = F, ops = "I",
            reduce.ranges = T, drop.empty.ranges = F,
            pos = mcols(gal)$pos
        )
        cigs_query <- cigarRangesAlongQuerySpace(cigar(gal),
            ops = "I", with.ops = F,
            reduce.ranges = T, drop.empty.ranges = F
        )
        all_ins <- mclapply(c(1:length(cigs_query)), function(i) {
            extractAt(mcols(gal)$seq[i], cigs_query[[i]])[[1]]
        })

        # Merge all insertions
        all_ins_merged <- do.call("rbind", mclapply(c(1:length(cigs_ref)), function(i) {
            return(data.frame(
                start_ref = start(cigs_ref[[i]]), end_ref = end(cigs_ref[[i]]),
                start_query = start(cigs_query[[i]]), end_query = end(cigs_query[[i]]),
                ins_seq = all_ins[[i]], width_ins = width(all_ins[[i]])
            ))
        },
        mc.cores = ncores
        ))
        all_ins_merged <- all_ins_merged[order(all_ins_merged$end_ref), ]

        # write.csv(all_ins_merged,'./testing/all_ins.csv',row.names=F);

        # TO DO: Check for overlaps--should be minimal since scaffolds don't usually overlap that much
        if (any(table(all_ins_merged$start_ref) > 1)) {
            print("Overlapping insertions")
            # not the best way, but just pick the first for now
            all_ins_merged <- all_ins_merged[!duplicated(all_ins_merged[, c("start_ref", "end_ref")]), ]
        }

        # Now the beauty part of inserting the strings back in
        # Split ref seq by the insert positions
        if (nrow(all_ins_merged) != 0) {
            new_strs <- DNAStringSet(rep("", nrow(all_ins_merged) + 1))
            for (i in 1:nrow(all_ins_merged)) {
                if (i == 1) {
                    new_strs[i] <- paste0(
                        extractAt(con_seq, IRanges(start = 1, end = all_ins_merged$end_ref[i]))[[1]],
                        all_ins_merged$ins_seq[i]
                    )
                } else {
                    new_strs[i] <- paste0(
                        extractAt(con_seq, IRanges(
                            start = all_ins_merged$start_ref[i - 1],
                            end = all_ins_merged$end_ref[i]
                        ))[[1]],
                        all_ins_merged$ins_seq[i]
                    )
                }
            }

            # Last bit
            new_strs[i + 1] <- paste0(extractAt(con_seq, IRanges(
                start = all_ins_merged$start_ref[i],
                end = width(con_seq)
            ))[[1]])
            temp_str <- paste0(as.character(new_strs), collapse = "")

            # Remove gaps to get final sequence
            con_seq_final <- DNAStringSet(gsub("-", "", temp_str))

            # No insertions
        } else {
            con_seq_final <- con_seq
        }
        names(con_seq_final) <- sub(".bam", "_consensus", basename(bamfname))

        writeXStringSet(
            con_seq_final,
            paste0("./", sampname, "_consensus.fasta")
        )

        # Delete bai file
        file.remove(baifname)
    } else {
        print("Bam file could not be opened.")
        return(NA)
    }
}


# Make a new reference scaffold
newref <- make_ref_from_assembly(bamfname, reffname)

if (newref == FALSE) print("Failed to generate consensus from scaffolds")