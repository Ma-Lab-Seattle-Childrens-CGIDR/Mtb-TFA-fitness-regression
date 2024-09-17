#!/usr/bin/env python

import argparse
import sys

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans


# Adapted from https://github.com/SBRG/pymodulon/blob/master/src/pymodulon/core.py
def _kmeans_cluster(m_df, imodulon):
	data = m_df[imodulon]

	df = pd.DataFrame(data.abs())

	model = KMeans(n_clusters=3, random_state=1)
	model.fit(df.values.reshape(-1, 1))

	df['cluster'] = model.labels_

	# Get top two clusters
	counts = df['cluster'].value_counts().sort_values(ascending=True)
	idx1 = counts.index[0]
	idx2 = counts.index[1]
	clust1 = df[df.cluster == idx1]
	clust2 = df[df.cluster == idx2]

	# Get midpoint between lowest iModulon gene and highest insignificant
	# gene
	threshold = np.mean([clust1[imodulon].min(), clust2[imodulon].max()])
	return threshold


def main(in_file, tf_file, out_file):
	m_df = pd.read_csv(in_file, index_col=0)

	m_norm_df = m_df.reset_index(names='Gene').melt( # make df with cols Gene, iModulon, Weight
		id_vars='Gene', var_name='iModulon', value_name='Weight')
	thresholds = m_df.columns.to_series().map(lambda imod_name: _kmeans_cluster(m_df, imod_name))

	m_norm_df['Threshold'] = m_norm_df['iModulon'].map(thresholds) # set thresholds as calc'd above

	m_norm_df = m_norm_df.loc[
		m_norm_df['Weight'].gt(m_norm_df['Threshold']), # only keep members passing weight threshold
		~m_norm_df.columns.isin(['Threshold']) # get rid of Threshold col since we're done with it
	]

	tf_df = pd.read_csv(tf_file, header=None, names=['Gene'], index_col=False)

	m_tf_df = m_norm_df[m_norm_df['Gene'].isin(tf_df['Gene'])] # separate df for only TFs

	edge_df = m_norm_df.merge( # merge them together to produce a df of TF-gene edges
		m_tf_df.rename(columns={
			'Gene': 'Regulator', # since these are the TFs, rename the label to Regulator
			'Weight': 'RegWt',   # and the Weight to RegWt, too
		}),
		how='right', on='iModulon' # map all TFs to co-iModulon genes (via right outer merge)
	)

	edge_df = edge_df[edge_df['Gene'].ne(edge_df['Regulator'])] # drop self-regulatory edges

	edge_df['Score'] = edge_df['Weight'].mul(edge_df['RegWt']) # acct for both wts in final score

	# output_df = edge_df[['Regulator', 'Gene', 'Score']].sort_values('Score', ascending=False)
	output_df = edge_df.groupby(
		['Regulator', 'Gene'] # if same regulator and gene show up in multiple modulons
	)[[
		'Score'
	]].sum( # add scores if same edge occurs multiple times
	).sort_values(
		'Score', ascending=False # then reset the sort
	).reset_index() # and get a fresh index

	output_df.to_csv(out_file, sep='\t', header=False, index=False)


if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		description=('Utility for transforming the results of an imodulon execution into an '
			'inferred transcriptional regulatory network in the style of the DREAM5 challenge.'))

	parser.add_argument('in_file',
		help=('The absolute or relative filename where the relevant imodulon data can be found. '
			'Should be an `M.csv` file.'))
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
