# inPhynite

A C++ tool for Bayesian phylogenetic inference under the infinite-sites mutation model. Created by Ivan Specht (ispecht@stanford.edu).

## Installation

### Prerequisites
- C++ compiler with C++17 support (e.g., clang++ or g++)
- Make

### Building from source

1. Clone the repository:
```bash
git clone https://github.com/yourusername/inPhynite.git
cd inPhynite
```

2. Compile the code:
```bash
make
```

This will create the `inPhynite` executable in the current directory.

## Usage

inPhynite can operate in two modes: analyzing real sequence data from a FASTA file, or analyzing simulated data.

### Basic command structure
```bash
./inPhynite [options]
```

### Required parameters (both modes)

- `--rooted <0|1>` - Whether the root is fixed (0=false, 1=true). If true and using real data, the first sequence in the FASTA is taken to be the root.
- `--n_global_iters <value>` - Number of global iterations
- `--n_tree_iters <value>` - Number of tree iterations per global iteration
- `--sample_every <value>` - Sample every N iterations
- `--infer_nt <0|1>` - Whether to infer effective population size (0=false, 1=true)
- `--out_dir <path>` - Output directory for results

### Mode 1: Real data (FASTA input)

When using real sequence data, provide:

- `--in_file <path>` - Input FASTA file with aligned sequences

**Example:**
```bash
./inPhynite \
  --in_file sequences.fasta \
  --rooted 0 \
  --n_global_iters 10000 \
  --n_tree_iters 10000 \
  --sample_every 10 \
  --infer_nt 1 \
  --out_dir ./output \
  --seed 12345
```

### Mode 2: Simulated data

When simulating data, provide:

- `--n <value>` - Number of leaves
- `--nt_init <value>` - Initial effective population size
- `--nt_log_slope <value>` - Log slope for effective population size over time

**Example:**
```bash
./inPhynite \
  --n 40 \
  --nt_init 1.0 \
  --nt_log_slope 0.0 \
  --rooted 0 \
  --n_global_iters 10000 \
  --n_tree_iters 10000 \
  --sample_every 10 \
  --infer_nt 0 \
  --out_dir ./output \
  --seed 12345
```

### Optional parameters

- `--seed <value>` - Random seed (default: 123)

### Input file format

For real data mode, the input FASTA file should contain aligned sequences. The tool will then:
1. Filter to biallelic sites (positions with exactly two nucleotides across all sequences)
2. Remove sites with missing data (N's or gaps)
3. Apply compatibility filtering to ensure data conform to the infinite-sites model

### Output files

The tool generates the following files in the output directory:

**For both modes:**
- `trees.nexus` - Sampled phylogenetic trees in NEXUS format
- `nt.csv` - Sampled effective population sizes over time

**For simulated data only:**
- `true_tree.nwk` - The true simulated tree in Newick format
- `true_nt.csv` - The true effective population sizes used in simulation
- `sequences.fasta` - Generated sequence data (100,000 bp)

## Cleaning build files

To remove compiled objects and the executable:
```bash
make clean
```

## License

Copyright 2025 Ivan Specht

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## Citation

Specht, I. & Palacios, J.A. Efficient Bayesian Phylogenetics under the Infinite-Sites Model with inPhynite. [arXiv preprint coming soon].