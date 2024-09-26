# Model Prediction

The code in this folder corresponds to section **2.9 Calculating transcription factor activity profiles from network component analysis** under Methods.

The input data found under `in\` are taken from Supplementary Tables 1, 2, and 4. The RNG networks under `out\cache\net` are the 10 randomized networks referred to in **3.3 Network component analyses reveal per-sample Mtb TF activities under different conditions** under Results and in Supplementary Figure 6.
Other networks and sequence data can be used by placing them in their respective folders under `in\`.
The NCA code in `nca.py` was taken from Covert Lab's [Whole Cell Model](https://github.com/CovertLab/WholeCellEcoliRelease/blob/master/reconstruction/ecoli/scripts/nca/nca.py)

## Prerequisites:
- Python 3.10 or later

The following libraries can be installed simply by running:
```sh
pip install libraryname
```
- numpy
- scipy
- swiglpk

## Usage
Running `python main.py -h` will return the following help message:
```
usage: main.py [-h] [-v] (-n [NETFILE] | -c [CACHENET]) -s SEQFILE [-q SEQLABEL] [-t [NETLABEL]] [-l] [-r]
               [-m {robust_nca,constrained_nca,fast_nca,random_nca}] [-i ROBUST_ITERATIONS]

options:
  -h, --help            show this help message and exit
  -v, --verbose         If set, prints status updates while running.
  -n [NETFILE], --netfile [NETFILE]
                        Absolute path of the .txt network file, or file name if located in the 'in\network' directory (default: Table4_aggregate.txt).
  -c [CACHENET], --cachenet [CACHENET]
                        Absolute path of the .npy network matrix file, or file name if located in the 'out\cache\net' directory (default:
                        network.npy).
  -s SEQFILE, --seqfile SEQFILE
                        Absolute path of the .tsv sequencing data file, or file name if located in the 'in\seq' directory.
  -q SEQLABEL, --seqlabel SEQLABEL
                        Absolute path of the .txt file containing row/gene labels for the seq data, or file name if located in the 'in\seq\labels'    
                        directory. Required if the rows of the network and seq matrix do not match.
  -t [NETLABEL], --netlabel [NETLABEL]
                        Absolute path of the .txt file containing row/gene labels for the network, or file name if located in the
                        'out\cache\net\labels' directory (default: netLabels.txt). Required if the rows of the network and seq matrix do not match      
                        while also using a cached network matrix.
  -l, --linear          If set, use linear counts from sequencing data, otherwise keep log2 counts.
  -r, --randomize       If set, will randomize the TF-gene pairings when generating the network matrix.
  -m {robust_nca,constrained_nca,fast_nca,random_nca}, --method {robust_nca,constrained_nca,fast_nca,random_nca}
                        NCA method to use, defined in nca.py (default: robust_nca).
  -i ROBUST_ITERATIONS, --robust-iterations ROBUST_ITERATIONS
                        Maximum iterations for robust_nca (default: 1000).
```
Example: `python main.py -n -s Table2_RNAseq.tsv -q Table2_RNAseq.txt`

The above command applies a network (which defaults to `Table4_aggregate.txt`) to the sequence data in `Table2_RNAseq.tsv`. The `-q` seqlabel option is required since the network and sequence data matrices do not have a matching number of rows.


## Contact
Oliver Gu - NJMS Yang Lab: [og127@njms.rutgers.edu](og127@njms.rutgers.edu)