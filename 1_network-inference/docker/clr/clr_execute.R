library("BSplineMI")
library("optparse")
library("parmigene")
library("tidyr")

#
# Define, parse, and validate arguments
#

option_list <- list(
	make_option(c("--expression_sep"), action="store", type="character", default="\t",
		help="The separator character delimiting the `expression_file`. Default is \"\\t\"."),
	make_option(c("--n_cores"), type="integer", default=6,
		help="The number of CPU cores that should be used when inferring the regulatory network. Default is %default."),
	make_option(c("--rows"), action="store", type="character", default="genes",
		help="What dimension the row-names should be considered to be: \"genes\" or \"samples\".")
)

parser <- OptionParser(
	usage="%prog [options] expression_file regulators_file out_file",
	option_list=option_list,
	description=paste(
		"Run CLR, a mutual-information-based regulatory network inference tool, ",
		"on supplied gene expression data.\n",
		"\n",
		"See the following for more information:\n",
		"- https://search.r-project.org/CRAN/refmans/parmigene/html/clr.html\n",
		"- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1764438/\n",
		"- https://web.archive.org/web/20100727172659/http://gardnerlab.bu.edu/data/PLoS_2007/CLR.html\n",
		sep=""
	)
)

args <- parse_args(parser, positional_arguments=3)
option <- args$options
files <- args$args

if (option$rows != "genes"  && option$rows != "samples") {
	stop("Please specify the orientation of the data by passing `rows={\"genes\",\"samples\"}`.")
}

#
# Read and massage data
#

expression_matrix <- read.table(files[1], header=TRUE, sep=option$expression_sep, as.is=TRUE)
regulators_list <- readLines(files[2], warn=FALSE)

# coerce matrix to numeric

col_names <- colnames(expression_matrix)
row_names <- rownames(expression_matrix)
expression_matrix <- as.matrix(expression_matrix)
expression_matrix <- apply(expression_matrix, 2, function(x) { as.numeric(x) })
colnames(expression_matrix) <- col_names
rownames(expression_matrix) <- row_names

if (option$rows == "samples") {
	expression_matrix <- t(expression_matrix)
}

#
# Run CLR to infer regulatory network
#

Sys.setenv(OMP_NUM_THREADS=option$n_cores)

# mutual_info_mat  <- knnmi.all(expression_matrix) # uses K-nearest neighbor distances

# "All mutual information values were computed using 10 bins and third order B-splines."
# - from CLR paper, Faith 2007 (PMID 17214507)
mutual_info_mat <- calcSplineMI(expression_matrix, 10, 3, threads=option$n_cores)

full_wt_adjacency_mat <- clr(mutual_info_mat)

full_wt_adjacency_df <- as.data.frame(full_wt_adjacency_mat)

tfs_wt_adjacency_df <- full_wt_adjacency_df[regulators_list]
tfs_wt_adjacency_df$Gene <- rownames(tfs_wt_adjacency_df)
rownames(tfs_wt_adjacency_df) <- NULL

#
# Extract edges, label, subset and write
#

edges <- reshape(tfs_wt_adjacency_df,
	varying=regulators_list, v.names="Score", idvar="Gene", direction="long")

edges$Combo <- rownames(edges)
rownames(edges) <- NULL
edges <- edges |>
	separate_wider_delim(Combo, delim=".", names=c("Gene_dup", "RegulatorIdx"))

edges$Regulator <- regulators_list[as.numeric(unlist(edges$RegulatorIdx))]

# subset
edges <- edges[edges$Score > 0, c("Regulator", "Gene", "Score")]

write.table(edges[order(-edges$Score),],
	file=files[3], quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE, fileEncoding='UTF-8')
