#!/usr/bin/env Rscript

## setup -----------------------------------------------------------------------
rm(list=ls(all.names=TRUE))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(argparse))
options(stringsAsFactors=FALSE)

# arguments --------------------------------------------------------------------
args <- ArgumentParser()
args$add_argument("-i", "--ids", required=TRUE, 
                  help="file with sample IDs")
args$add_argument("-f", "--frontclip", required=TRUE, type="integer",
                  help="bases to clip from the front")
args$add_argument("-b", "--backclip", required=TRUE, type="integer",
                  help="bases to clip from the back")
args$add_argument("-o", "--odir", default=".", 
                  help="output directory")
args$add_argument("-t", "--threshold", default=0, type="integer",
                  help="SNV difference threshold for clustering")
args <- args$parse_args()

# helper functions -------------------------------------------------------------
mesg <- function(...) {
  x <- paste("echo -e '[MSG]", ..., "'")
  system(x)
}
err <- function(...) {
  x <- paste("echo -e '[ERR]", ..., "'")
  system(x)
  q(status=1)
}
read.vcf <- function(fname, front=args$frontclip, back=args$backclip) {
  fname %>%
    read.delim(comment.char="#",
               col.names=c("Segment", "Position", "SNV", "Ref", "Alt", 
                           "Qual", "Filter", "Info"),
               colClasses=c(Ref="character", Alt="character")) %>%
    filter(Position > front,
           Position < back) %>%
    mutate(Frequency=as.numeric(str_extract(Info, "(?<=AF=)[0-9\\.]+")),
           SNV=paste0(Segment, "-", Position, "-", Ref, "-", Alt),
           Filename=fname) %>%
    select(Filename, SNV, Frequency)
}

## check inputs ----------------------------------------------------------------
# check output directory
if(dir.exists(args$odir)) {
  mesg("Valid output directory:", args$odir)
} else {
  mesg("Creating output directory:", args$odir)
  dir.create(args$odir, recursive=TRUE)
}

# check and load sample info
if(file.exists(args$ids)) {
  mesg("Valid sample ID file:", args$ids)
  # first column is samples
  ids <- read.csv(args$ids, col.names=c("ID", "Filename"))
  mesg("Sample file contains", dim(ids)[1], "samples")
} else {
  err("Invalid sample ID file:", args$ids)
}

# check and load VCF files
vcf.exists <- file.exists(ids$Filename)
if(all(vcf.exists)) {
  mesg("All VCF files found; loading data...")
  snvs <- lapply(ids$Filename, read.vcf) %>%
          do.call(rbind, .) %>%
          # remove SNVs with frequency < 50%
          filter(Frequency >= 0.5) %>%
          # change to Sample ID
          left_join(ids, by="Filename") %>%
          select(SNV, ID, Frequency)
} else {
  missing <- ids$Filename[!which(vcf.exists)] %>%
             paste(collapse=", ")
  err("Missing VCF files detected:", missing)
}
rm(vcf.exists)

# done checking inputs
mesg("Done loading inputs!\n")

## pairwise comparison ---------------------------------------------------------
# fill in zeros and transform to frequency matrix
snvs <- expand.grid(SNV=unique(snvs$SNV),
                    ID=unique(ids$ID)) %>%
        left_join(snvs, by=c("SNV", "ID")) %>%
        # transform to wide format 
        reshape2::dcast(SNV ~ ID, fill=0, value.var="Frequency") %>%
        column_to_rownames("SNV") %>%
        as.matrix()

# change to binary; only had SNVs with >50%
snvs[snvs >= 0.5] <- 1
mode(snvs) <- "integer"

mesg("Calculating pairwise differences. This may take a while...")

# init difference matrix
nsamp <- dim(snvs)[2]
difs <- matrix(nrow=nsamp, ncol=nsamp, 
               dimnames=list(colnames(snvs), 
                             colnames(snvs)))
mode(difs) <- "integer"

# loop over samples (index=1, vector=v1)
for(i in 1:nsamp) {
  # get first vector
  v1 <- snvs[, i]
  
  # loop over top diagonal (index=2, vector=v2)
  for(j in i:nsamp) {
    # get second vector
    v2 <- snvs[, j]
    # calculate differences
    d <- sum(v1 != v2)
    
    # add value to both top and bottom diagonals of the matrix
    difs[i, j] <- d
    difs[j, i] <- d
  }
}
mode(difs) <- "integer"

diff.file <- paste0(args$odir, "/differences.csv")
mesg("Done calculating differences. Saving to", diff.file, "\n")
write.csv(difs, diff.file)

# clean up
rm(diff.file, d, i, j, nsamp, v1, v2, snvs, ids)

## assign clusters -------------------------------------------------------------
mesg("Assigning clusters using", args$threshold, "difference threshold")

# init cluster data
nsamp <- dim(difs)[1]
x <- data.frame(Rows=rownames(difs),
                RowIndex=1:nsamp)
y <- data.frame(Columns=colnames(difs),
                ColumnIndex=1:nsamp)
rm(nsamp)

# basing the cluster ID on the column
clusters <- difs %>%
            as.data.frame() %>%
            rownames_to_column("Rows") %>%
            reshape2::melt(id.vars="Rows",
                           variable.name="Columns",
                           value.name="Differences") %>%
            # add row indices
            full_join(x, by="Rows") %>% 
            # add column indices 
            full_join(y, by="Columns") %>%
            # remove matches to self (indices are equal) 
            # and where row index > column index
            filter(ColumnIndex > RowIndex) %>% 
            # get minimum difference 
            slice_min(by=Columns, n=1, order_by=Differences, with_ties=FALSE) 
rm(x, y)

# samples that start a new cluster will have minimum differences > threshold
# the first sample (column index 1) is an exception; it will always have ID=1
parent <- clusters %>%
          filter(Differences > args$threshold) %>% 
          # subset to column and index only
          select(Columns, ColumnIndex) %>%
          # add column name 1 (edge case)
          rbind(data.frame(Columns=colnames(difs)[1],
                           ColumnIndex=1)) %>%
          # arrange by column index
          dplyr::arrange(ColumnIndex) %>%
          # update sample ID column to "ID"
          rename(ID=Columns,
                 Index=ColumnIndex) %>%
          # set column order
          select(ID, Index)
# add cluster ID column
parent$Cluster <- 1:dim(parent)[1]

# subset clusters to those with matches <= threshold
clusters <- clusters %>%
            filter(Differences <= args$threshold) %>%
            # here the "Row" column indicates what matches to clusters$ID
            rename(ID=Rows) %>%
            select(ID, Columns, ColumnIndex) %>%
            left_join(parent, by="ID") %>%
            # format so that we only have Column, ColumnIndex, and Cluster
            select(-ID, -Index) %>% 
            rename(ID=Columns,
                   Index=ColumnIndex) %>%
            # set column order
            select(ID, Index, Cluster)

# combine parent clusters and child clusters
clusters <- rbind(parent, clusters) %>%
            # order by index
            dplyr::arrange(Index) %>%
            # remove index column 
            select(-Index)
rm(difs, parent)

# save cluster assignment
clust.file <- paste0(args$odir, "/clusters.csv")
mesg("Done assigning clusters; saving to", clust.file, "\n")
write.csv(clusters, clust.file, row.names=FALSE)

## all done! -------------------------------------------------------------------
mesg("Clustering complete! Outputs can be found in", args$odir)
sessionInfo()
