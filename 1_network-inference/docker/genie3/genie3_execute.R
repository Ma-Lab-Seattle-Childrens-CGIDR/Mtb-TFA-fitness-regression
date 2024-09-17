library("GENIE3")
library("optparse")

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
		"Run GENIE3, a random-forest-based regulatory network inference tool, ",
		"on supplied gene expression data. Running time is ~14h for an E. coli regulatory network.\n",
		"\n",
		"See the following for more information:\n",
		"- https://bioconductor.org/packages/release/bioc/html/GENIE3.html\n",
		"- https://doi.org/10.1038/nmeth.4463\n",
		"- https://github.com/aertslab/GENIE3\n",
		"- https://doi.org/10.1371/journal.pone.0012776\n",
		"- https://github.com/vahuynh/GENIE3/tree/master/GENIE3_R_C_wrapper\n",
		"- https://raw.githubusercontent.com/vahuynh/GENIE3/master/GENIE3_running_times_DREAM5_nov16.png",
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

# Rules for expression input file:
# (as per https://github.com/vahuynh/GENIE3/blob/856416b/GENIE3_R_C_wrapper/GENIE3.R#L293)
# colnames and rownames must be directly inferable by read.table
# colnames and rownames must all be unique
# matrix data must be coercible to numeric values
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
# Run GENIE3 to infer regulatory network
#

weight_matrix <- GENIE3(expression_matrix, nCores=option$n_cores, regulators=regulators_list, verbose=TRUE)
link_list <- getLinkList(weight_matrix, reportMax=100000)

write.table(link_list, file=files[3], quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE, fileEncoding='UTF-8')
