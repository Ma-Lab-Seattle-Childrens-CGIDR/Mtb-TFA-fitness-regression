#!/usr/bin/env python

import argparse
import sys

import pandas as pd


def main(in_file, tf_file, out_file):
	cmonkey_db_df = pd.read_csv(in_file,
		header=None, names=['Gene', 'Cluster', 'Iteration', 'Residual'])
	cmonkey_db_df['ClusterScore'] = 1 - cmonkey_db_df['Residual']

	tf_df = pd.read_csv(tf_file, header=None, names=['Gene'], index_col=False)

	iteration_cluster_tfs_df = cmonkey_db_df[cmonkey_db_df['Gene'].isin(tf_df['Gene'])]

	network_df = cmonkey_db_df.merge(
		iteration_cluster_tfs_df[['Gene', 'Cluster', 'Iteration']].rename(
			columns={'Gene': 'Regulator'}),
		how='outer', on=['Cluster', 'Iteration']).dropna()

	score_df = network_df.groupby(['Regulator', 'Gene'])[['ClusterScore']].sum().rename(
		columns={'ClusterScore': 'Sum'}).reset_index()
	score_df['PossibleCount'] = score_df['Regulator'].map(
		iteration_cluster_tfs_df.groupby('Gene')['Residual'].count())
	score_df['Score'] = score_df['Sum'].div(score_df['PossibleCount'])
	score_df.drop(index=score_df[score_df['Regulator'].eq(score_df['Gene'])].index, inplace=True)

	score_df['RelScore'] = score_df['Score'].sub(score_df['Score'].min()).div(
		score_df['Score'].max() - score_df['Score'].min())

	score_df[['Regulator', 'Gene', 'RelScore']].sort_values('RelScore', ascending=False).to_csv(
		out_file, sep='\t', header=False, index=False)


if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		description=('Utility for transforming the results of a cMonkey2 execution into an '
			'inferred transcriptional regulatory network in the style of the DREAM5 challenge. '
			'Execution will normally go through `cmonkey_process.sh`.'))

	parser.add_argument('in_file',
		help=('The absolute or relative filename where the relevant cMonkey data can be found. See '
			'`cmonkey_process.sh` for details on what this file should contain.'))
	parser.add_argument('tf_file',
		help=('The absolute or relative filename where the relevant regulators can be found. '
			'The regulators file should contain a newline-delimited list of regulator '
			'(e.g. transcription factor) labels in the same format as the column headers '
			'in the expression file.'))
	parser.add_argument('out_file',
		help='The absolute or relative filename where the output should be written.')

	args = parser.parse_args(sys.argv[1:])
	args_dict = vars(args)
	main(**args_dict)
