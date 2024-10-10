#!/usr/bin/env python

import argparse
import numpy as np
import os
import sys
import warnings

warnings.filterwarnings("ignore", "\nPyarrow", DeprecationWarning)
import pandas as pd

class LengthMatchError(ValueError):
	pass

def aggregate(network_dfs, methods=[], length_cap=100_000, threshold=2,
		verbosity=0):

	if not methods:
		methods = ['rank_average'] * len(network_dfs)

	if len(methods) != len(network_dfs):
		raise LengthMatchError

	truncate_by_fraction = False
	if length_cap <= 0:
		length_cap = sys.maxsize
	elif length_cap < 1:
		truncate_by_fraction = True

	aggregate_df = pd.DataFrame(columns=['Regulator', 'Gene', 'transform'])

	for i, network_df in enumerate(network_dfs):
		if truncate_by_fraction:
			length_cap = floor(len(network_df) * length_cap)

		network_df = network_df.iloc[:min(len(network_df), length_cap), :].reset_index()

		network_df['transform'] = transform(network_df, methods[i])

		aggregate_df = pd.concat((aggregate_df, network_df[['Regulator', 'Gene', 'transform']]))

	aggregate_df2 = aggregate_df.groupby(['Regulator', 'Gene']).agg(['sum', 'count'])
	aggregate_df2['Score'] = aggregate_df2[('transform', 'sum')]

	aggregate_df2 = aggregate_df2[aggregate_df2[('transform', 'count')].ge(threshold)]

	aggregate_df2 = aggregate_df2[['Score']]

	return aggregate_df2.sort_values('Score', ascending=False).reset_index()

def direct(network_df, direction_df, direction_score_col='Score'):
	if 'Direction' not in direction_df.columns:
		direction_df['Direction'] = direction_df[direction_score_col].map(
			lambda v: 'up' if v > 0 else 'down')

	direction_df = direction_df[
		['Regulator', 'Gene', 'Direction']
	].groupby(['Regulator', 'Gene']).max()

	network_df = network_df.set_index(['Regulator', 'Gene'])

	shared_idx = network_df.index.intersection(direction_df.index)

	network_df.loc[shared_idx, 'Direction'] = direction_df.loc[shared_idx, 'Direction']
	network_df = network_df[['Score', 'Direction']].dropna()

	return network_df.sort_values('Score', ascending=False).reset_index()

def transform(network_df, method):
	if method == 'rank_average':
		if len(network_df) > 10_000_000:
			raise ValueError((
				'Current rank average method will not work for networks larger than 10M edges ('
				f'attempting with a network of size {len(network_df)})'
			))
		return 10_000_000 - network_df.index.to_series()

	return None

#
# main
#

def main(network_files, direction_file='', out_file='network.txt', methods=[], length_cap=100_000,
		threshold=2, verbosity=0):

	network_dfs = []

	for network_file in network_files:
		network_dfs.append(
			pd.read_csv(network_file,
				sep='\t', names=['Regulator', 'Gene', 'Score'], index_col=False)
		)

	try:
		aggregate_df = aggregate(network_dfs,
			methods=methods, length_cap=length_cap, threshold=threshold, verbosity=verbosity)
	except LengthMatchError:
		print('Supplied methods must be the same in number as supplied network files.')
		sys.exit(1)

	if direction_file:
		direction_df = pd.read_csv(direction_file,
			sep='\t', names=['Regulator', 'Gene', 'value'], index_col=False)

		aggregate_df = direct(aggregate_df, direction_df, direction_score_col='value')

	aggregate_df.to_csv(out_file, sep='\t', header=False, index=False)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		description=('Utility for aggregating transcriptional regulatory networks generated from '
			'different sources into a single network.'))

	parser.add_argument('network_files',
		nargs='+',
		help=('The absolute or relative filenames where the relevant network data can be found. '
			'The network files should conform to the format expected by the DREAM5 challenge, '
			'a tab-delimited text file with 3 columns: regulator, regulated gene, score. Scores '
			'should fall between 0 and 1, with higher values indicating more confidence that the '
			'given regulatory relationship is valid.'))

	parser.add_argument('-d', '--direction_file',
		help=('The absolute or relative filename where a directed network can be found, to assign '
			'directionality to the resulting network. Contents should conform to the format '
			'expected by the DREAM5 challenge, but with the numeric third field containing signed '
			'values indicating the directionality of the network. If no direction_file is '
			'supplied, the resulting network will be undirected.'))

	parser.add_argument('-o', '--out_file',
		default='network.txt',
		help=('The absolute or relative filename where the aggregated network should be written. '
			'Contents will conform to the DREAM5 challenge format, as described under '
			'`network_files`.'))

	parser.add_argument('-l', '--length_cap',
		default=100_000,
		type=float,
		help=('The maximum length of an input network. Any additional edges beyond this count will '
			'be truncated. A value less than or equal to 0 will leave each input network '
			'untruncated. A value between 0 and 1 will be intepreted as a fraction of the original '
			'network. Defaults to 100,000.'))
	parser.add_argument('-t', '--threshold',
		default=2,
		type=int,
		help=('The number of individual networks a gene must be a part of in order to be included '
			'in the final aggregated network. Default is 2.'))

	# parser.add_argument('-m', '--methods',
	# 	nargs='+',
	# 	help='The methods that should be used to aggregate the networks together. WIP.')
	parser.add_argument('-v', '--verbosity',
		action='count',
		default=0,
		help=('The level of verbosity, i.e., the amount of processing information that should be '
			'output in the course of aggregating networks.'))

	args = parser.parse_args(sys.argv[1:])
	args_dict = vars(args)
	main(**args_dict)
