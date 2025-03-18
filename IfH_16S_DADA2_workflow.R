##### Analysis script for Pacbio CCS full-length 16S sequences #####
### Start script by running Rscript IfH_16S_DADA2_workflow.R [input path][Run ID]
### Input folder must contain CCS files in the following format: SampleID-16S-PB1.fastq.gz
### Output will contain ASV- & taxa-table as both .tsv and .Rdata files, a sequence list in .fna format & several QC files

##### Setup environment #####

args <- commandArgs(trailingOnly = TRUE)

path <- args[1]
projectname <- args[2]

output <- "/mnt/Netzlaufwerke/ukmsrv1488NGSequencing/16S/DADA2_pipeline_outputs"
SILVA <- "/mnt/Netzlaufwerke/ukmsrv1488NGSequencing/16S/silva_nr99_v138.2_toSpecies_trainset.fa.gz" 

#.libPaths(c("/usr/local/lib/R/site-library", "/usr/lib/R/site-library", "/usr/lib/R/library"))

# load libraries
library(dada2)
library(ggplot2)
library(logr)
library(stringr)

# specify primers
F27 <- "AGRGTTYGATYMTGGCTCAG" #these are the PacBio default primers for FL 16S (protocol version Feb 2020)
R1492 <- "RGYTACCTTGTTACGACTT"

# create output folders
t <- format(Sys.time(), "%F_%H-%M-%S")
Out <- file.path(output, paste("Out", projectname, t, sep = "_"))
QC <- file.path(Out, paste("QC", projectname, t, sep = "_"))
dir.create(Out)
dir.create(QC)

#set theme to apply to all plots (ggplot2)
theme_set(theme_bw()) 

### Import data ###

# get sample names from file names
fn <- list.files(path, full.names=TRUE)
samples <- basename(fn)
samples <- str_split_fixed(samples, "-", n=2)[,1]

# Check if any samples contain 0 Reads
empties <- character()
for (file in fn){
  n_reads <- length(ShortRead::readFastq(file))
  if (n_reads == 0){
    empties <- append(empties, file)
  }
}

empty_names <- str_split_fixed(basename(empties), "-", n=2)[,1]
fn <- fn[! fn %in% empties]
nz_samples <- samples[! samples %in% empty_names]

### Open a log-file and start time tracking
time.start <- Sys.time()
log_open(file_name = file.path(QC, "Logfile.log"), logdir = FALSE, show_notes = FALSE) #open log-file to track time

##### Run DADA2-Pipeline #####

# Remove primers and orient reads, discards reads that do not contain primers
nop <- file.path("dada_tmp", "noprimers", basename(nz_samples)) #creates a temporary subfolder in working directory
prim <- removePrimers(fn, nop, primer.fwd=F27, primer.rev=dada2:::rc(R1492), orient=TRUE, verbose=TRUE)
noprimers <- nop[file.exists(nop) == FALSE] #checks if every sample has reads and discards samples without reads
if (length(noprimers) > 0){
  np_names <- str_split_fixed(basename(noprimers), "-", n=2)[,1]
}
nop <- nop[file.exists(nop) == TRUE] 

t_prim <- Sys.time()

# QC step: Create a histogram of length distribution of oriented reads, should peak around ~1450
lens.fn <- lapply(nop, function(fn) nchar(getSequences(fn)))
lens <- do.call(c, lens.fn)
png(file = file.path(QC, "Length_distribution.png"))
hist(lens, 100)
dev.off()

# Filter reads by length & quality score
filts <- file.path("dada_tmp", "noprimers", "filtered", basename(nop))
reads.filt <- filterAndTrim(nop, filts, minQ=3, minLen=1000, maxLen=1600, maxN=0, rm.phix=FALSE, maxEE=2)
t_filt <- Sys.time()

# Check if any samples contain 0 Reads after filtering
empty_names_filt <- nz_samples[! nz_samples %in% basename(filts)]
nz_samples_filt <- nz_samples[! nz_samples %in% empty_names_filt]

### Run DADA2
# Dereplicate
drp <- derepFastq(filts, verbose=TRUE)
t_drp <- Sys.time()

# Learn error rates
err <- learnErrors(drp, errorEstimationFunction=PacBioErrfun, BAND_SIZE=32, multithread=TRUE)
t_err <- Sys.time()

# QC step: Plot error rates
png(file = file.path(QC, "Error_rates_plot.png"))
plotErrors(err)
dev.off()

# Denoise
dd <- dada(drp, err=err, BAND_SIZE=32, multithread=TRUE)
t_dd <- Sys.time()

# Make sequence table
seqtab <- makeSequenceTable(dd)

# Remove chimeras
st.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dimensions <- dim(st.nochim)

# QC step: Track reads through each processing step
prim1 <- prim[prim[,2] != 0,] #removes 0-Read samples
read.track <- cbind(ccs=prim1[,1], primers=prim1[,2], filtered=reads.filt[,2], denoised=sapply(dd, function(x) sum(x$denoised)), nonchim=rowSums(st.nochim))
read.track <- cbind(read.track, retained = round((read.track[,5]/read.track[,1]), digits = 2))
rownames(read.track) <- nz_samples_filt

### Save outputs
save(st.nochim, file = file.path(Out, "DADA2_out_seqtab.RData"))
write.table(t(st.nochim), file.path(Out, "ASV_table.tsv"), sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
uniquesToFasta(st.nochim, fout = file.path(Out, "Sequence_list.fna"), ids=colnames(st.nochim))

# remove temporary directory
unlink("dada_tmp", recursive = TRUE)
t_dada <- Sys.time()

### Assign taxonomy
taxa <- assignTaxonomy(st.nochim, SILVA, multithread=TRUE)
t_taxo <- Sys.time()
write.table(t(taxa), file.path(Out, "Taxa_table.tsv"), sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
save(taxa, file = file.path(Out, "taxa.RData"))

time.end <- Sys.time()
runtime <- round(difftime(time.end, time.start, units = "auto"), digits = 2)

### Write Log-file for QC
sep(paste("Run ID: ", projectname))

sep("Package versions")

for (package_name in sort(.packages())) {
  put(paste(package_name, packageVersion(package_name)), blank_after = FALSE )
}

sep("Input data")
put("Forward primer:", blank_after = FALSE)
put(F27)
put("Reverse primer:", blank_after = FALSE)
put(R1492)
put("Database for taxonomic classification:", blank_after = FALSE)
put(basename(SILVA))
put("Filtering parameters:", blank_after = FALSE)
put("minQ=3, minLen=1000, maxLen=1600, maxN=0, rm.phix=FALSE, maxEE=2")
put("Samples:", blank_after = FALSE)
put(samples)
put("The following samples were removed from the analysis because they contained 0 Reads:", blank_after = FALSE)
if(exists("empty_names")){
  put(empty_names, blank_after = FALSE)
}
if(exists("np_names")){
  put(np_names)
}
if(! exists("empty_names") && ! exists("np_names")){
  put("none")
}

sep("Quality control")
put(paste("The median read length after removing primers was:", median(lens)))
put("Read tracking through pipeline:", blank_after = FALSE)
put(read.track)
put(paste("On average,", round((mean(read.track[,6]))*100, digits = 2), "% of the reads were retained."))

put(paste(dimensions[2], "unique ASVs were found in", dimensions[1], "samples."))

sep("Time tracking")
put(paste("Start:", time.start), blank_after = FALSE)
put(paste("Primers removed:", t_prim), blank_after = FALSE)
put(paste("Filtered and trimmed:", t_filt), blank_after = FALSE)
put(paste("Dereplicated:", t_drp), blank_after = FALSE)
put(paste("Error rates learned:", t_err), blank_after = FALSE)
put(paste("Denoised:", t_dd), blank_after = FALSE)
put(paste("Dada completed:", t_dada), blank_after = FALSE)
put(paste("Taxonomy assigned:", t_taxo), blank_after = FALSE)
put(paste("Finished:", time.end), blank_after = TRUE)
put(runtime)

log_close()
unlink(file.path(QC, "Logfile.msg"), recursive = TRUE)
