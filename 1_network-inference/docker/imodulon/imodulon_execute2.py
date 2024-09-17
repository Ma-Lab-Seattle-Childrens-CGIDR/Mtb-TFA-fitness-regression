# Adapted from https://github.com/avsastry/modulome-workflow/blob/main/4_optICA

import argparse
import os
import shutil
import sys
import time
import warnings

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy import sparse, stats
from sklearn.cluster import DBSCAN
from sklearn.decomposition import PCA, FastICA
from sklearn.exceptions import EfficiencyWarning

def main(expression_file, iterations=1, tolerance=1e-7, out_dir='', distance=0.1,
		minimum_fraction=0.5, dimensionality='scan', dim_begin=20, dim_end=0, dim_step=0,
		finalize=True):

	# Set output files
	if not out_dir:
		out_dir = os.getcwd()
	else:
		if not os.path.exists(out_dir):
			os.makedirs(out_dir)

	out_subdir1 = os.path.join(out_dir, 'ica_runs')
	if not os.path.exists(out_subdir1):
		os.makedirs(out_subdir1)

	if not finalize == 'only':
		data_df = pd.read_csv(expression_file, sep='\t').transpose()
		n_genes, m_samples = data_df.shape

		if dimensionality == 'scan':
			if not dim_end:
				dim_end = m_samples
			if not dim_step:
				dim_step = (dim_end - dim_begin) // 25

			dimension_counts = range(dim_begin, dim_end, dim_step)
		elif dimensionality == 'infer':
			# Reduce dimensionality using PCA
			pca = PCA().fit(data_df.transpose())
			pca_var = np.cumsum(pca.explained_variance_ratio_)
			k_comp = np.where(pca_var > 0.99)[0][0] + 1
			print('Data: {} genes x {} samples'.format(n_genes, m_samples))
			print('Found {} dimensions from PCA'.format(k_comp))

			dimension_counts = [k_comp]
		else:
			raise ValueError(
				f'Dimensionality of "{dimensionality}" is invalid. See help text for more details.')

		for n_dimensions in dimension_counts:
			# Create temporary directory for files
			out_subdir2 = os.path.join(out_subdir1, str(n_dimensions))
			if not os.path.exists(out_subdir2):
				os.makedirs(out_subdir2)

			tmp_dir = os.path.join(out_subdir2, 'tmp')
			if not os.path.exists(tmp_dir):
				os.makedirs(tmp_dir)

			for i in range(iterations):
				run_ica(data_df, tolerance, n_dimensions, tmp_dir, i)

			compute_distances(tmp_dir, iterations)

			cluster_components(out_subdir2, tmp_dir, iterations, distance, minimum_fraction)

	if not finalize == False: # truthy check here would suffice, but check explicitly for clarity
		get_dimension(out_dir)


def run_ica(data_df, tolerance, n_dimensions, tmp_dir, i):
	# Run ICA

	S = []
	A = []

	ica = FastICA(whiten=True, max_iter=int(1e10), tol=tolerance, n_components=n_dimensions)
	S = pd.DataFrame(ica.fit_transform(data_df), index=data_df.index)
	A = pd.DataFrame(ica.mixing_, index=data_df.columns)

	S.to_csv(os.path.join(tmp_dir, 'proc_{}_S.csv'.format(i)))
	A.to_csv(os.path.join(tmp_dir, 'proc_{}_A.csv'.format(i)))

	print(f'Finished round {i+1} of ICA')


def compute_distances(tmp_dir, iterations):
	for i in range(iterations):
		for j in range(i, iterations):
			S1 = pd.read_csv(os.path.join(tmp_dir, 'proc_{}_S.csv'.format(i)), index_col=0)
			S2 = pd.read_csv(os.path.join(tmp_dir, 'proc_{}_S.csv'.format(j)), index_col=0)
			dist = abs(np.dot(S1.T, S2))
			dist[dist < 0.5] = 0
			dist_file = os.path.join(tmp_dir, 'dist_{}_{}.npz'.format(i, j))
			sparse_dist = sparse.coo_matrix(np.clip(dist, 0, 1))
			sparse.save_npz(dist_file, sparse_dist)

	print('Finished computing distances between ICA runs')


def cluster_components(out_dir, tmp_dir, iterations, distance, minimum_fraction):
	min_samples = int(np.ceil(args.minimum_fraction * iterations))

	# -----------------------------------------------------------
	# Combine distance matrix

	print('Combining D matrix')

	block = []
	block_size = {}
	for i in range(iterations):
		col = []
		for j in range(iterations):
			if i <= j:
				mat = sparse.load_npz(os.path.join(tmp_dir, 'dist_{}_{}.npz'.format(i, j)))
				col.append(mat)
				block_size[i] = mat.shape[0]
				block_size[j] = mat.shape[1]
			else:
				mat = sparse.load_npz(os.path.join(tmp_dir, 'dist_{}_{}.npz'.format(j, i)))
				col.append(mat.T)
		block.append(col)
	D = sparse.bmat(block, 'csr')

	# ----------------------------------------------------------
	# Convert to true distance matrix (originally dissimilarity)

	data = D.data
	indices = D.indices
	indptr = D.indptr
	D = sparse.csr_matrix((1 - D.data, D.indices, D.indptr))

	print('Clustering')

	# Run DBSCAN
	with warnings.catch_warnings():
		warnings.filterwarnings('ignore', category=EfficiencyWarning)
		dbscan = DBSCAN(eps=distance, min_samples=min_samples, metric='precomputed')
		labels = dbscan.fit_predict(D)
		n_clusters = max(labels) + 1

	print('Identified', n_clusters, 'clusters')

	# -----------------------------------------------------------

	# Place clustered components into correct bins

	print('Loading individual S and A matrices')
	start = 0
	end = 0
	S_bins = {i: [] for i in range(n_clusters)}
	A_bins = {i: [] for i in range(n_clusters)}

	for i in range(iterations):
		# Get labels for each partial matrix
		start = end
		end += block_size[i]
		proc_labels = labels[start:end]

		# Add parts of matrix to full table
		S_part = pd.read_csv(
			os.path.join(tmp_dir, 'proc_{}_S.csv'.format(i)), index_col=0
		)
		A_part = pd.read_csv(
			os.path.join(tmp_dir, 'proc_{}_A.csv'.format(i)), index_col=0
		)
		S_part.columns = range(S_part.shape[1])
		A_part.columns = range(A_part.shape[1])
		for i, label in enumerate(proc_labels):
			if label != -1:
				S_bins[label].append(S_part[i])
				A_bins[label].append(A_part[i])

	# ----------------------------------------------------------

	# Gather final S and A matrices

	print('Gathering final S and A matrices')

	S_final = pd.DataFrame(columns=range(n_clusters), index=S_part.index)
	A_final = pd.DataFrame(columns=range(n_clusters), index=A_part.index)
	df_stats = pd.DataFrame(
		columns=['S_mean_std', 'A_mean_std', 'count'], index=range(n_clusters)
	)

	for label in range(n_clusters):
		S_clust = S_bins[label]
		A_clust = A_bins[label]

		# First item is base component
		Svec0 = S_clust[0]
		Avec0 = A_clust[0]

		# Make sure base component is facing positive
		if abs(min(Svec0)) > max(Svec0):
			Svec0 = -Svec0
			Avec0 = -Avec0

		S_single = [Svec0]
		A_single = [Avec0]

		# Add in rest of components
		for j in range(1, len(S_clust)):
			Svec = S_clust[j]
			Avec = A_clust[j]
			if stats.pearsonr(Svec, Svec0)[0] > 0:
				S_single.append(Svec)
				A_single.append(Avec)
			else:
				S_single.append(-Svec)
				A_single.append(-Avec)

		# Add centroid of cluster to final S matrix
		S_final[label] = np.array(S_single).T.mean(axis=1)
		A_final[label] = np.array(A_single).T.mean(axis=1)

		# Get component stats
		df_stats.loc[label, 'S_mean_std'] = np.array(S_single).std(axis=1).mean()
		df_stats.loc[label, 'A_mean_std'] = np.array(A_single).std(axis=1).mean()
		df_stats.loc[label, 'count'] = len(S_single)

	print('Final components created')

	# Write to file
	# Remove components that exist in under 50% of runs
	good_comps = df_stats[df_stats['count'] > (iterations * 0.5)].index

	S_final.columns = range(len(S_final.columns))
	A_final.columns = range(len(A_final.columns))

	print('Writing files to ' + out_dir)
	S_final.to_csv(os.path.join(out_dir, 'M.csv'))
	A_final.T.to_csv(os.path.join(out_dir, 'A.csv'))

	# Clean up tmp directories
	try:
		shutil.rmtree(tmp_dir)
	except PermissionError:
		print(f'(failed to remove temp directory {tmp_dir})')

	print('Clustering complete!')


def get_dimension(out_dir):
	print('Computing optimal set of independent components')

	def load_mat(dim, mat):
		df = pd.read_csv(
			os.path.join(out_dir, 'ica_runs', str(dim), mat + '.csv'), index_col=0
		)
		df.columns = range(len(df.columns))
		return df.astype(float)

	dims = sorted([int(x) for x in os.listdir(os.path.join(out_dir, 'ica_runs'))])

	M_data = [load_mat(dim, 'M') for dim in dims]
	A_data = [load_mat(dim, 'A') for dim in dims]

	# -----------------------------------
	# Check large iModulon dimensions

	final_a = A_data[-1]
	while np.allclose(final_a, 0, atol=0.01):
		A_data = A_data[:-1]
		M_data = M_data[:-1]
		dims = dims[:-1]
		final_a = A_data[-1]

	final_m = M_data[-1]

	n_components = [m.shape[1] for m in M_data]

	# -----------------------------------
	# Get iModulon statistics

	thresh = 0.7

	n_final_mods = []
	n_single_genes = []
	for m in M_data:
		# Find iModulons similar to the highest dimension
		l2_final = np.sqrt(np.power(final_m, 2).sum(axis=0))
		l2_m = np.sqrt(np.power(m, 2).sum(axis=0))
		dist = (
			pd.DataFrame(abs(np.dot(final_m.T, m)))
			.divide(l2_final, axis=0)
			.divide(l2_m, axis=1)
		)
		n_final_mods.append(len(np.where(dist > thresh)[0]))

		# Find iModulons with single gene outliers
		counter = 0
		for col in m.columns:
			sorted_genes = abs(m[col]).sort_values(ascending=False)
			if sorted_genes.iloc[0] > 2 * sorted_genes.iloc[1]:
				counter += 1
		n_single_genes.append(counter)

	non_single_components = np.array(n_components) - np.array(n_single_genes)

	DF_stats = pd.DataFrame(
		[n_components, n_final_mods, non_single_components, n_single_genes],
		index=[
			'Robust Components',
			'Final Components',
			'Multi-gene Components',
			'Single Gene Components',
		],
		columns=dims,
	).T
	DF_stats.sort_index(inplace=True)

	dimensionality = (
		DF_stats[DF_stats['Final Components'] >= DF_stats['Multi-gene Components']]
		.iloc[0]
		.name
	)

	print('Optimal Dimensionality:', dimensionality)

	# Plot dimensions
	fig, ax = plt.subplots()
	ax.plot(dims, n_components, label='Robust Components')
	ax.plot(dims, n_final_mods, label='Final Components')
	ax.plot(dims, non_single_components, label='Non-single-gene Components')
	ax.plot(dims, n_single_genes, label='Single Gene Components')

	ax.vlines(dimensionality, 0, max(n_components), linestyle='dashed')

	ax.set_xlabel('Dimensionality')
	ax.set_ylabel('# Components')
	ax.legend(bbox_to_anchor=(1, 1))
	plt.savefig(
		os.path.join(out_dir, 'dimension_analysis.pdf'),
		bbox_inches='tight',
		transparent=True,
	)

	# Save final matrices
	final_M_file = os.path.join(out_dir, 'ica_runs', str(dimensionality), 'M.csv')
	final_A_file = os.path.join(out_dir, 'ica_runs', str(dimensionality), 'A.csv')

	final_M_dest = os.path.join(out_dir, 'M.csv')
	final_A_dest = os.path.join(out_dir, 'A.csv')
	shutil.copyfile(final_M_file, final_M_dest)
	shutil.copyfile(final_A_file, final_A_dest)

	# Clean up tmp directories
	try:
		shutil.rmtree(os.path.join(out_dir, 'ica_runs'))
	except PermissionError:
		print(f'(failed to remove ica_runs directory)')


if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		description=('Performs ICA multiple times on a dataset to identify robust independent '
			'components. The output files are stored in a temporary directory between '
			'processing steps.'))

	parser.add_argument('-d', '--distance',
		default=0.1,
		type=float,
		help='Maximum distance between points in a cluster for DBSCAN.')
	parser.add_argument('-f', '--expression_file',
		required=True,
		help='Path to expression data file.')
	parser.add_argument('-i', '--iterations',
		default=1,
		type=int,
		help='Number of ICA runs.')
	parser.add_argument('-m', '--minimum_fraction',
		default=0.5,
		type=float,
		help='Minimum fraction of samples in a cluster for DBSCAN.')
	parser.add_argument('-o', '--out_dir',
		default='',
		help='Path to output file directory. Defaults to current working directory.')
	parser.add_argument('-t', '--tolerance',
		default=1e-7,
		type=float,
		help=('ICA convergence tolerance. Larger values run faster but provide less accurate '
			'components.'))

	parser.add_argument('-n', '--dimensionality',
		choices=['infer', 'scan'],
		default='scan',
		help=('Method of determining dimensionality for search. "infer" will determine the optimal '
			'number of dimensions by running PCA. "scan" will test a large variety of dimensions '
			'iteratively and aggregate the results together.'))
	parser.add_argument('--dim_begin',
		default=20,
		type=int,
		help='Minimum dimensionality for ICA search (default: 20).')
	parser.add_argument('--dim_end',
		default=0,
		type=int,
		help='Maximum dimensionality for ICA search (default: n_samples).')
	parser.add_argument('--dim_step',
		default=0,
		type=int,
		help='Dimensionality step size (default: n_samples / 25).')

	parser.add_argument('--finalize',
		choices=[True, False, 'only'],
		default=True,
		type=lambda f: True if f == 'True' else False if f == 'False' else f,
		help=('Whether or not the results should be aggregated together (when '
			'dimensionality=scan). "true" will aggregate all results generated in the given '
			'execution; this is the default. "false" will cause no aggregation to be performed. '
			'"only" will skip all other processing steps and only aggregate any results present in '
			'the appropriate location: this allows for parallelizing many executions with '
			'finalize=false and aggregating all the results with finalize=only.'))

	args = parser.parse_args(sys.argv[1:])
	args_dict = vars(args)
	main(**args_dict)
