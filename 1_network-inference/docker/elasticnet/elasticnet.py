#!/usr/bin/env python

import argparse
import numpy as np
import os
from multiprocessing import Pool
import sys
import warnings

warnings.filterwarnings("ignore", "\nPyarrow", DeprecationWarning)
import pandas as pd

from sklearn.exceptions import ConvergenceWarning
import sklearn.linear_model as skllm
import sklearn.model_selection as sklms

def regress(tf_expr_df, gene, gene_expr, l1_ratio=[0.5, 0.8, 0.9, 0.95, 0.98, 0.99],
		penstrength_subdivisions=1000, n_cores=4, cv_folds=4, verbosity=0):

	x = tf_expr_df.to_numpy(dtype=np.float64)

	# re-scale phenotypes to make coefficients comparable across phenotypes: we will use to score
	gene_expr = (gene_expr - gene_expr.mean(axis=0)) / gene_expr.std(axis=0)
	y = gene_expr.to_numpy(dtype=np.float64)

	elasticnet_cv = skllm.ElasticNetCV(l1_ratio=l1_ratio, n_alphas=penstrength_subdivisions,
		cv=cv_folds, n_jobs=n_cores, verbose=verbosity-2)

	with warnings.catch_warnings():
		warnings.simplefilter('ignore', ConvergenceWarning)
		elasticnet_cv.fit(x, y)

	coef = elasticnet_cv.coef_
	r_squared = elasticnet_cv.score(x, y)

	if np.count_nonzero(coef) == 0:
		if verbosity >= 2:
			print(f'No non-zero coefficients for {gene}; moving on...')
		return [], [], r_squared

	nonzero_idxs = np.not_equal(coef, 0)
	predictor_columns = tf_expr_df.columns.to_series().to_numpy()

	if verbosity >= 2:
		print(f'Non-zero coefficients ({np.count_nonzero(coef)}) for {gene}:',
			', '.join((
				f'{col} ({coef})' for col, coef in zip(
					predictor_columns[nonzero_idxs], coef[nonzero_idxs])
			)),
			f'| R^2={r_squared:.2f}, l1_ratio={elasticnet_cv.l1_ratio_},',
			f'penalty_strength={elasticnet_cv.alpha_}')

	return predictor_columns[nonzero_idxs], coef[nonzero_idxs], r_squared


#
# main
#

def main(expression_file, regulators_file, out_file, expression_sep='\t', n_cores=4, cv_folds=4,
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

	#
	# run regression
	#

	associations = []

	for gene in expression_df.columns:
		assoc_tfs, coefs, r_squared = regress(regulators_df.drop(columns=gene, errors='ignore'),
			gene, expression_df[gene], n_cores=n_cores, cv_folds=cv_folds, verbosity=verbosity)

		for tf, coef in zip(assoc_tfs, coefs):
			associations.append((tf, gene, coef, r_squared))

	#
	# process and save results
	#

	results_df = pd.DataFrame(data=associations, columns=['Regulator', 'Gene', 'Coef', 'R_squared'])
	# results_df['Score'] = results_df['Coef'].mul(results_df['R_squared']).abs()
	results_df['Score'] = results_df['Coef'].abs()

	output_df = results_df[['Regulator', 'Gene', 'Score']].sort_values(by='Score', ascending=False)

	output_df.to_csv(out_file, sep='\t', header=False, index=False)
	results_df.to_csv('_details.'.join(out_file.rsplit('.', maxsplit=1)))

	if verbosity >= 1:
		print(f'Results written to file {out_file} and more info written to sibling _details file.')

if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		description=('Utility for inferring transcriptional regulatory networks from gene '
			'expression data using elastic-net regularization.'))

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
	parser.add_argument('-f', '--cv_folds',
		default=4,
		type=int,
		help=('Number of folds to use for x-fold cross-validation when evaluating regulation '
			'models for each gene.'))
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
