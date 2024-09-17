#!/usr/bin/env python

import argparse
import numpy as np
import os
from multiprocessing import Pool
import sys
import warnings

warnings.filterwarnings("ignore", "\nPyarrow", DeprecationWarning)
import pandas as pd

import statsmodels.api as sm
from statsmodels.tools import add_constant

def regress(tf_expr_df, gene, gene_expr, n_cores=4, verbosity=0):

	x = tf_expr_df.to_numpy(dtype=np.float64)
	y = gene_expr.to_numpy(dtype=np.float64)

	# regression = skllm.LinearRegression(copy_X=False, n_jobs=n_cores).fit(x, y)
	model = sm.OLS(y, add_constant(x), hasconst=True)
	fit = model.fit()
	summary = fit.summary()

	coefs = fit.resid
	r_squared = fit.rsquared
	p_vals = fit.pvalues

	predictor_columns = tf_expr_df.columns.to_series().to_numpy()

	if verbosity >= 3:
		print(f'Coefficients for {gene} (R^2={r_squared:.2f}):',
			', '.join((
				f'{col} ({coef})' for col, coef in zip(predictor_columns, coefs)
			)),
		)
	elif verbosity >= 2:
		print(f'Model for {gene} completed (R^2={r_squared:.2f}).')

	return predictor_columns, coefs, p_vals, r_squared

#
# main
#

def main(expression_file, regulators_file, out_file, expression_sep='\t', n_cores=4,
		rows='samples', verbosity=0):

	if rows == 'genes':
		header = None
		index_col = 0
		transpose = True
	elif rows == 'samples':
		header = 0
		index_col = False
		transpose = False

	expression_df = pd.read_csv(expression_file,
		sep=expression_sep, header=header, index_col=index_col)

	if transpose:
		expression_df = expression_df.T

	if verbosity >= 1:
		print(f'Using expression data with {expression_df.shape[0]} samples and',
			f'{expression_df.shape[1]} genes.')

	regulators = pd.read_csv(regulators_file, header=None, index_col=False)[0]

	if verbosity >= 1:
		print(f'Generating predictions for {len(regulators)} regulators,',
			f'starting with {regulators.iloc[0]} and ending with {regulators.iloc[-1]},',
			f'using {n_cores} CPU cores.')

	regulators_df = expression_df[regulators.to_list()]

	# re-scale predictors to avoid large variance leading to over-representation in coefficients
	regulators_df = (regulators_df - regulators_df.mean(axis=0)) / regulators_df.std(axis=0)

	associations = []

	for gene in expression_df.columns:
		assoc_tfs, coefs, p_vals, r_squared = regress(regulators_df.drop(columns=gene, errors='ignore'),
			gene, expression_df[gene], n_cores=4, verbosity=verbosity)

		for tf, coef, p_val in zip(assoc_tfs, coefs, p_vals):
			associations.append((tf, gene, coef, p_val, r_squared))

	results_df = pd.DataFrame(
		data=associations, columns=['Regulator', 'Gene', 'Coef', 'p-val', 'R_squared'])

	results_df.to_csv(out_file, sep='\t', index=False)

	if verbosity >= 1:
		print(f'Results written to file {out_file}.')

if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		description=('Utility for inferring transcriptional regulatory networks from gene '
			'expression data using a simple linear model.'))

	parser.add_argument('expression_file',
		help=('The absolute or relative filename where the relevant expression data can be found. '
			'The expression file should contain a character-delimited matrix of expression values, '
			'with orientation corresponding to that indicated in `--rows`. '
			'Gene labels are required; sample labels should not be included.'
			'Delimiter can be configured with `--expression_sep`.'))
	parser.add_argument('regulators_file',
		help=('The absolute or relative filename where the relevant regulators can be found. '
			'The regulators file should contain a newline-delimited list of regulator '
			'(e.g. transcription factor) labels in the same format as the column headers '
			'in the expression file.'))
	parser.add_argument('out_file',
		help='The absolute or relative filename where the output should be written.')

	parser.add_argument('-s', '--expression_sep',
		default='\t',
		help='The delimiter of the data in `expression_file`.')
	parser.add_argument('-c', '--n_cores',
		default=4,
		type=int,
		help=('Number of CPU cores to use for processor-intensive operations.'))
	parser.add_argument('-r', '--rows',
		choices=('genes', 'samples'),
		default='samples',
		help='The dimension that rows should be considered to represent: \"genes\" or \"samples\".')
	parser.add_argument('-v', '--verbosity',
		action='count',
		default=0,
		help=('The level of verbosity, i.e., the amount of processing information that should be '
			'output in the course of inferring a network.'))

	args = parser.parse_args(sys.argv[1:])
	args_dict = vars(args)
	main(**args_dict)
