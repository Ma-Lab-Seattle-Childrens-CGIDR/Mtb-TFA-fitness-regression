#!/usr/bin/env python

import argparse
from datetime import datetime
import itertools
from sklearn.metrics import auc
import numpy as np
import os
import random
import sys
import warnings

warnings.filterwarnings("ignore", "\nPyarrow", DeprecationWarning)
import pandas as pd

class Logger:
	levels = (
		'RESULT:',
		'WARNING:',
		'INFO:',
		'DEBUG:'
	)

	def debug(self, message, verbosity):
		self.log(message, level=3, verbosity=verbosity)

	def info(self, message, verbosity):
		self.log(message, level=2, verbosity=verbosity)

	def log(self, message, level, verbosity):
		if verbosity >= level:
			print(f'{datetime.now().isoformat(sep="_")} {Logger.levels[level]:<8} {message}')

	def result(self, message, verbosity):
		self.log(message, level=0, verbosity=verbosity)

	def warn(self, message, verbosity):
		self.log(message, level=1, verbosity=verbosity)

log = Logger()

def calc_aupr(actual_df, full_prediction_df, space_len, step_size=-1, verbosity=0):
	if step_size <= 0:
		step_size = space_len // 100
		log.info(f'Using default step size of {space_len} / 100 = {step_size}.', verbosity)

	precisions, recalls = [], []
	for i in range(1, space_len, step_size):
		true_pos, false_pos, false_neg, true_neg = confusion_matrix(
			actual_df['Edge'], full_prediction_df.iloc[0:i,3], space_len)

		precision = true_pos / (true_pos + false_pos)
		precisions.append(precision)
		recall = true_pos / (true_pos + false_neg)
		recalls.append(recall)

		log.debug(f'Iteration={i: 7d}; precision={precision:.2%}; recall={recall:.2%}', verbosity)

	return auc(recalls, precisions)

def calc_false_neg(actual_df, full_prediction_df, space_len, step_size=-1, verbosity=0):
	if step_size <= 0:
		step_size = space_len // 100
		log.info(f'Using default step size of {space_len} / 100 = {step_size}.', verbosity)

	recalls = []
	for i in range(1, space_len, step_size):
		true_pos, false_pos, false_neg, true_neg = confusion_matrix(
			actual_df['Edge'], full_prediction_df.iloc[0:i,3], space_len)

		recall = true_pos / (true_pos + false_neg)
		recalls.append(recall)

		log.debug(f'Iteration={i: 7d}; size={i}; recall={recall:.2%}', verbosity)

	return auc(list(range(1, space_len, step_size)), recalls)

# positive likelihood ratio
# def calc_plr(true_pos, false_pos, false_neg, true_neg):
# 	true_pos_rate = true_pos / (true_pos + false_neg)
# 	false_pos_rate = false_pos / (false_pos + true_neg)
# 	return true_pos_rate / false_pos_rate

def confusion_matrix(actual, predicted, space_len):
	true_predictions_bool = predicted.isin(actual)
	true_pos = true_predictions_bool.sum()
	false_pos = true_predictions_bool.eq(False).sum()

	false_neg = actual.isin(predicted).eq(False).sum()
	true_neg = space_len - true_pos - false_pos - false_neg

	return true_pos, false_pos, false_neg, true_neg

def fill_network(network_df, tfs, genes, score=0, verbosity=0):
	all_edges = [
		(edge[0], edge[1], f'{edge[0]}_{edge[1]}') for edge in itertools.product(tfs, genes)
	]

	missing_edges = [edge for edge in all_edges if edge[2] not in network_df['Edge']]
	random.shuffle(missing_edges)

	missing_df = pd.DataFrame(missing_edges, columns=['Regulator', 'Gene', 'Edge'])
	missing_df.insert(2, 'Score', score)

	return pd.concat((network_df, missing_df)).drop_duplicates('Edge'), len(missing_df)

def matrix_pretty_print(confsn_matrix, label, prediction_label='Predicted'):
	print(label + ':')
	print(f'           | {prediction_label} positive | {prediction_label} negative |')
	print(f'Actual pos | True pos  = {confsn_matrix[0]:6g} | False neg = {confsn_matrix[2]:6g} |')
	print(f'Actual neg | False pos = {confsn_matrix[1]:6g} | True neg  = {confsn_matrix[3]:6g} |')
	print(f'Sensitivity = {confsn_matrix[0] / (confsn_matrix[0] + confsn_matrix[2]):.3g}')
	print(f'Specificity = {confsn_matrix[3] / (confsn_matrix[1] + confsn_matrix[3]):.3g}')

def random_network(tfs, genes, size):
	edges = []
	for i in range(size):
		tf = tfs[random.randrange(len(tfs))]
		gene = genes[random.randrange(len(genes))]
		edges.append((tf, gene, tf + '_' + gene))

	return pd.DataFrame(data=edges, columns=['Regulator', 'Gene', 'Edge'])

# def validate(pr(chip_df2, best_rand_net, space_len)
# 	# print('Random AUC:', random_auc)

# 	return best_matrix

#
# main
#

def main(chip_file, network_file, method='aupr', step_size=-1, verbosity=-1):
	log.info(f'Reading ChIP data from file {chip_file}.', verbosity)

	chip_df = pd.read_csv(chip_file, sep='\t', header=None,
		names=['Regulator', 'Gene', 'Score'])
	chip_df['Edge'] = chip_df['Regulator'].str.cat(chip_df['Gene'], sep='_')
	chip_df.drop_duplicates('Edge', inplace=True)

	log.info(f'Using {len(chip_df)} unique rows from file {chip_file}.', verbosity)
	log.info(f'Reading network data from file {network_file}.', verbosity)

	network_df = pd.read_csv(network_file, sep='\t', header=None,
		names=['Regulator', 'Gene', 'Score'])
	network_df['Edge'] = network_df['Regulator'].str.cat(network_df['Gene'], sep='_')
	network_df.drop_duplicates('Edge', inplace=True)

	log.info(f'Using {len(network_df)} unique rows from file {network_file}.', verbosity)

	mtb_tf_df = pd.read_csv(
		os.path.join(
			os.path.dirname(os.path.realpath(__file__)),
			'mount/mtb_tfs_20240220.txt'
		), header=None, names=['Gene'])
	mtb_genome_df = pd.read_excel(
		os.path.join(
			os.path.dirname(os.path.realpath(__file__)),
			'../src/datasets/MtbGeneAnnotation2015.xlsx'
		))

	network_filled_df, num_new_edges = fill_network(
		network_df, mtb_tf_df['Gene'], mtb_genome_df['Locus Tag'], score=-1, verbosity=verbosity)

	log.info(f'Added {num_new_edges} edges, randomly sorted, to complete the network.', verbosity)

	if method == 'aupr':
		log.info(f'Calculating AUPR of the supplied network.', verbosity)

		result = calc_aupr(chip_df, network_filled_df, len(network_filled_df), step_size=step_size,
			verbosity=verbosity)

		log.result(
			f'Supplied network in {network_file} achieves AUPR of {result} against ChIP data in ' + \
				f'{chip_file}.',
			verbosity)
	elif method == 'false-neg':
		log.info(f'Calculating false negative rate of the supplied network.', verbosity)

		result = calc_false_neg(chip_df, network_filled_df, len(network_filled_df),
			step_size=step_size, verbosity=verbosity)

		log.result(
			f'Supplied network in {network_file} achieves AUPR of {result} against ChIP data in ' + \
				f'{chip_file}.',
			verbosity)

	# generate many random networks (cache them!), and get a p-value for enrichment of this network
	# 	2k? 5k? 10k? 25k?
	# 	maybe only cache a pdf of this distribution?
	# 	sklearn.feature_selection import f_regression?

	# output network AUPR, average random network AUPR, and p-value
	return result

if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		description=('Utility for validating a provisional Mtb transcriptional regulatory network '
			'against ChIP-seq data.'))

	parser.add_argument('chip_file',
		help=('The absolute or relative filenames where the relevant ChIP-seq data can be found. '
			'The chip file should conform to the format expected by the DREAM5 challenge, a '
			'tab-delimited text file with 3 columns: regulator, regulated gene, score. Scores '
			'should fall between 0 and 1, and the file should be sorted with better scores first.'))
	parser.add_argument('network_file',
		help=('The absolute or relative filename where the relevant network data can be found. '
			'The network file should conform to the format expected by the DREAM5 challenge, '
			'a tab-delimited text file with 3 columns: regulator, regulated gene, score. Scores '
			'should fall between 0 and 1, with higher values indicating more confidence that the '
			'given regulatory relationship is valid.'))

	parser.add_argument('-s', '--step-size',
		default=-1,
		type=int,
		help=('The number of edges to add in each iteration when evaluating the AUPR of the '
			'supplied network. Configures the granularity of the resulting precision-recall curve, '
			'which has a slight effect on the AUPR. If no value is supplied, the system will '
			'default to 1%% of the network space, i.e. len(tfs) * len(genes) / 100, '
			'resulting in 100 points being measured for the PR curve.'))

	parser.add_argument('-v', '--verbosity',
		action='count',
		default=0,
		help=('The level of verbosity, i.e., the amount of processing that should be output in the '
			'course of processing. Default is 0, or only the result; 1 will yield warnings; 2, '
			'info messages; 3, debug messages.'))

	args = parser.parse_args(sys.argv[1:])
	args_dict = vars(args)
	main(**args_dict)
