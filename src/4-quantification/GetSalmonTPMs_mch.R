require("tidyverse")
require("tximport")
require("rtracklayer")
require("dplyr")

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Provide a metric to extract TPM counts", call.=FALSE)
}

salmon_out <- args[2]
design_tab <- args[3]
gtf <- args[4]
outname <- args[5]

# Prepare annotation ------------------------------------------------------

# Get transcripts
human_gtf <- rtracklayer::import(gtf)
transcripts <- human_gtf %>% as.data.frame(.) %>%
  dplyr::select(c(gene_id, transcript_id)) %>%
  dplyr::rename(
    GENEID = gene_id, 
    TXNAME = transcript_id) %>%
  dplyr::filter(!is.na(TXNAME))

gene_names <-  human_gtf %>% as.data.frame(.) %>%
  dplyr::select(c(gene_id, gene_name)) %>% 
  dplyr::rename(
    GENEID = gene_id,
    GENENAME = gene_name) %>%
  dplyr::mutate(GENEID = str_remove_all(
    GENEID, ".[0-9]*$")) %>%
  dplyr::distinct(GENEID, GENENAME)

# Final annotation
# human_annotation <- transcripts %>% 
#   dplyr::left_join(., gene_names, by="GENEID")

rm(human_gtf)

# Prepare data ------------------------------------------------------------

design <- read_delim(
    design_tab, 
    delim = "\t",
    escape_double = FALSE, 
    trim_ws = TRUE
)

files <- file.path(salmon_out, design$SampleName , "quant.sf")
names(files) <- paste0(design$SampleName)

# Prepare counts ----------------------------------------------------------

# Raw TPM counts
txi.inf <- tximport(files, type = "salmon", tx2gene = transcripts, dropInfReps=TRUE, countsFromAbundance=args[1], ignoreAfterBar=TRUE)
counts_df <- txi.inf$counts %>%
  as.data.frame() %>%
  rownames_to_column("GENEID") %>%
  dplyr::left_join(., gene_names, by="GENEID")

write.table(
  counts_df, 
  str_interp("${outname}_COUNTS_AGGREGATED-${args[1]}.tab"),
  sep = "\t", quote = FALSE, 
  row.names = FALSE
)
rm(txi.inf)