##### Demutiplexing of mouse tumor sample barcodes #####

library(deMULTIplex2)

##### Demutiplex starting with raw sequencing data #####

# Read filtered cell barcodes
cell_barcodes <- read.csv("/c4/home/liw3/LWHDRB2/processing_output/2_cellranger_output_cat/outs/filtered_feature_bc_matrix/barcodes.tsv.gz", sep = "\t", header = FALSE)

# Extract the nucleotides of the barcodes. 10X output for barcodes has the format of Nucleotides-1
cell_barcodes_extract <- stringr::str_split(cell_barcodes$V1, pattern = "-", simplify = TRUE)
cell_barcodes_input <- cell_barcodes_extract[,1]


# Read barcode library
read_table <- readTags(dir = "/c4/home/liw3/LWHDRB2/reads/BC_cat",
                       name = "LWHDRB2_Bar",
                       barcode.type = "MULTIseq",
                       assay = "RNA",
                       filter.cells = cell_barcodes_input)
# Align tags
data(multiseq_oligos)
tag.ref <- multiseq_oligos[21:32]
tag_mtx <- alignTags(read_table, tag.ref)

# Obtain demultiplex results
res <- demultiplexTags(tag_mtx,
                       plot.path = "/c4/home/liw3/LWHDRB2/processing_output/3_demultiplex_output_cat",
                       plot.name = "demux2_",
                       plot.diagnostics = T)

# Write result files
bar_assignment <- res$assign_table
write.csv(bar_assignment, file = "/c4/home/liw3/LWHDRB2/processing_output/3_demultiplex_output_cat/barcode_assignment.csv")
