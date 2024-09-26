import argparse
import os
import time
from pathlib import Path
from collections import defaultdict

import numpy as np

import nca


BASE_DIR = Path.cwd()

# Network related
NETWORK_DIR = BASE_DIR / 'in' / 'network'
NETWORK_FILE = NETWORK_DIR / 'Table4_aggregate.txt'

# Sequencing related
SEQ_DIR = BASE_DIR / 'in' / 'seq'
SEQ_LABELS_DIR = SEQ_DIR / 'labels'

# Output related
OUTPUT_DIR = BASE_DIR / 'out'
CACHE_DIR = OUTPUT_DIR / 'cache'
FILTERED_DIR = CACHE_DIR / 'filtered'
NET_CACHE_DIR = CACHE_DIR / 'net'
NET_CACHE_FILE = NET_CACHE_DIR / 'network.npy'
NET_LABELS_DIR = NET_CACHE_DIR / 'labels'
NET_LABELS_FILE = NET_LABELS_DIR / 'netLabels.txt'


def txt2npy(netFile: Path, randomize = False) -> tuple[dict[str, int], np.ndarray]:
    """
    Loads a network .txt file and converts it into a numpy array.

    Args:
        netFile: path to the network .txt file
        randomize: if True, randomizes the TF-gene pairings for each TF while preserving the total number of pairings

    Returns:
        rows: mapping of gene labels to row numbers of the network matrix
        net: the network matrix
    """
    
    genes = set()
    pairs: dict[str, dict[str, int]] = defaultdict(dict) # dict[TF, dict[gene, value]]
    rows: dict[str, int] # genes
    cols: dict[str, int] # TFs

    with open(netFile, 'r') as f:
        for line in f.readlines():
            tf, gene = line.split()[:2]
            genes.add(gene)
            pairs[tf].update({gene:1})
            # tf, gene, value, sign = line.split()
            # genes.add(gene)
            # value = float(value)
            # if sign == 'down':
            #     value *= -1
            # pairs[tf].update({gene:value})
            
    tfs = pairs.keys()

    # Sort using casefold to match MS Excel sorting
    rows = {gene:row for (row, gene) in enumerate(sorted(list(genes), key=str.casefold))}
    cols = {tf:col for (col, tf) in enumerate(sorted(list(tfs), key=str.casefold))}

    net = np.zeros((len(rows), len(cols))) # genes x TFs
    if randomize:
        rng = np.random.default_rng()

    for tf in pairs:
        col = cols[tf]
        if randomize:
            rngRows = rng.choice(len(rows), len(pairs[tf]), replace=False)
            for row in rngRows:
                net[row, col] = 1
        else:
            for gene, value in pairs[tf].items():
                row = rows[gene]
                net[row, col] = 1
    
    return genes, net

def filter_net(netLabels: set[str], sharedGenes: list[str], net: np.ndarray) -> np.ndarray:
    """
    Filters out any network matrix rows that are not shared with the seq matrix.

    Args:
        netLabels: gene labels for the network matrix rows
        sharedGenes: ordered list of genes common to both matrices
        net: the network matrix to be filtered
        
    Returns:
        filteredNet: the filtered network matrix
    """

    rows = {gene:row for (row, gene) in enumerate(sorted(list(netLabels), key=str.casefold))}
    rowIdxs = [rows[gene] for gene in sharedGenes]
    return net[rowIdxs, :]

def filter_seq(seqLabels: list[str], sharedGenes: list[str], seq_data: np.ndarray) -> np.ndarray:
    """
    Filters out any seq matrix rows that are not shared with the seq matrix.

    Args:
        seqLabels: gene labels for the seq matrix rows
        sharedGenes: ordered list of genes common to both matrices
        seq_data: the seq matrix to be filtered
        
    Returns:
        filteredSeq: the filtered seq matrix
    """

    rows = {gene:row for (row, gene) in enumerate(seqLabels)}
    rowIdxs = [rows[gene] for gene in sharedGenes]
    return seq_data[rowIdxs, :]

def parse_args() -> argparse.Namespace:
    """Parse command line args for options to run."""

    parser = argparse.ArgumentParser()

    default_nca = 'robust_nca'
    default_robust_iterations = 1000

    # General options
    # parser.add_argument('-o', '--output',
    #     default=OUTPUT_DIR,
    #     help=f'Absolute path of output directory to save results to (default: {OUTPUT_DIR.relative_to(BASE_DIR)}).')
    parser.add_argument('-v', '--verbose',
        action='store_true',
        help='If set, prints status updates while running.')

    # Required arguments
    network = parser.add_mutually_exclusive_group(required=True)
    network.add_argument('-n', '--netfile', nargs='?', const=NETWORK_FILE.name,
        help=f"Absolute path of the .txt network file, or file name if located in the '{NETWORK_DIR.relative_to(BASE_DIR)}' directory (default: {NETWORK_FILE.name}).")
    network.add_argument('-c', '--cachenet', nargs='?', const=NET_CACHE_FILE.name,
        help=f"Absolute path of the .npy network matrix file, or file name if located in the '{NET_CACHE_DIR.relative_to(BASE_DIR)}' directory (default: {NET_CACHE_FILE.name}).")
    
    parser.add_argument('-s', '--seqfile',
        required=True,
        help=f"Absolute path of the .tsv sequencing data file, or file name if located in the '{SEQ_DIR.relative_to(BASE_DIR)}' directory.")

    # Filtering options
    parser.add_argument('-q', '--seqlabel',
        help=f"Absolute path of the .txt file containing row/gene labels for the seq data, or file name if located in the '{SEQ_LABELS_DIR.relative_to(BASE_DIR)}' directory."
        ' Required if the rows of the network and seq matrix do not match.')
    parser.add_argument('-t', '--netlabel', nargs='?', const=NET_LABELS_FILE.name,
        help=f"Absolute path of the .txt file containing row/gene labels for the network, or file name if located in the '{NET_LABELS_DIR.relative_to(BASE_DIR)}' directory (default: {NET_LABELS_FILE.name})."
        ' Required if the rows of the network and seq matrix do not match while also using a cached network matrix.')
    
    # Data options
    parser.add_argument('-l', '--linear',
        action='store_true',
        help='If set, use linear counts from sequencing data, otherwise keep log2 counts.')
    parser.add_argument('-r', '--randomize',
        action='store_true',
        help='If set, will randomize the TF-gene pairings when generating the network matrix.')

    # NCA options
    parser.add_argument('-m', '--method',
        choices=nca.METHODS,
        default=default_nca,
        help=f'NCA method to use, defined in nca.py (default: {default_nca}).')

    ## ROBNCA specific options
    parser.add_argument('-i', '--robust-iterations',
        type=int,
        default=default_robust_iterations,
        help=f'Maximum iterations for robust_nca (default: {default_robust_iterations}).')

    return parser.parse_args()


if __name__ == '__main__':
    start = time.time()
    args = parse_args()

    os.makedirs(FILTERED_DIR, exist_ok=True)
    os.makedirs(NET_LABELS_DIR, exist_ok=True)

    if args.netfile:
        if args.verbose:
            print('Generating network matrix from network .txt file...')
        networkFile = Path(args.netfile)
        if not networkFile.is_absolute():
            networkFile = NETWORK_DIR / networkFile
        netLabels, net = txt2npy(networkFile, randomize=args.randomize)
        with open(NET_LABELS_DIR / 'networkLabels.txt', 'w') as f:
            for gene in netLabels:
                f.write(gene + '\n')
        np.save(NET_CACHE_FILE, net)
        np.savetxt(NET_CACHE_DIR / 'network.tsv', net, delimiter=',')
    elif args.cachenet:
        if args.verbose:
            print('Loading network matrix from cache...')
        networkFile = Path(args.cachenet)
        if not networkFile.is_absolute():
            networkFile = NET_CACHE_DIR / networkFile
        net = np.load(networkFile)

    if args.verbose:
        print('Loading sequencing data...')
    seqFile = Path(args.seqfile)
    if not seqFile.is_absolute():
        seqFile = SEQ_DIR / seqFile
    seq_data = np.loadtxt(seqFile)
    if args.linear:
        seq_data = 2**seq_data

    if args.seqlabel:
        if args.verbose:
            print('Filtering network and seq matrices to shared genes...')
        labelFile = Path(args.seqlabel)
        if not labelFile.is_absolute():
            labelFile = SEQ_LABELS_DIR / labelFile
        with open(labelFile, 'r') as f:
            seqLabels = {line.strip() for line in f.readlines()}
        
        if args.netlabel:
            labelFile = Path(args.netlabel)
            if not labelFile.is_absolute():
                labelFile = NET_LABELS_DIR / labelFile
            with open(labelFile, 'r') as f:
                netLabels = {line.strip() for line in f.readlines()}

        sharedGenes = sorted(list(netLabels & seqLabels), key=str.casefold)
        with open(FILTERED_DIR / 'rowLabels.txt', 'w') as f:
            for gene in sharedGenes:
                f.write(gene + '\n')

        net = filter_net(netLabels, sharedGenes, net)
        np.save(FILTERED_DIR / 'filteredNetwork.npy', net)
        np.savetxt(FILTERED_DIR / 'filteredNetwork.tsv', net, delimiter=',')
        
        seq_data = filter_seq(seqLabels, sharedGenes, seq_data)
        np.save(FILTERED_DIR / 'filteredSeqData.npy', seq_data)
        np.savetxt(FILTERED_DIR / 'filteredSeqData.tsv', seq_data, delimiter=',')
    
    print(f'Network matrix (A) shape: {net.shape}')
    print(f'Expression matrix (E) shape: {seq_data.shape}')

    # Solve NCA problem
    nca_method = getattr(nca, args.method)
    A, P = nca_method(seq_data, net, n_iters=args.robust_iterations)

    # Save results
    np.savetxt(OUTPUT_DIR / 'A.csv', A, delimiter=',')
    np.savetxt(OUTPUT_DIR / 'P.csv', P, delimiter=',')
    print(f'A matrix shape: {A.shape}')
    print(f'P matrix shape: {P.shape}')

    print(f'Completed in {(time.time() - start) / 60:.1f} min')
