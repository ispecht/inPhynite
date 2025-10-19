# inPhynite

A C++ tool for Bayesian phylogenetic inference using coalescent theory.

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

### Basic command structure
```bash
./inPhynite [options]
```

### Required parameters

- `--nt_init <value>` - Initial effective population size
- `--nt_log_slope <value>` - Log slope for nt
- `--n <value>` - Number of leaves
- `--rooted <0|1>` - Whether the root is fixed (0=false, 1=true)
- `--n_global_iters <value>` - Number of global iterations
- `--n_tree_iters <value>` - Number of tree iterations
- `--sample_every <value>` - Sample every N iterations
- `--infer_nt <0|1>` - Whether to infer nt (0=false, 1=true)
- `--out_dir <path>` - Output directory for results

### Optional parameters

- `--seed <value>` - Random seed (default: 123)

### Example

```bash
./inPhynite \
  --nt_init 10000 \
  --nt_log_slope 0.01 \
  --n 50 \
  --rooted 1 \
  --n_global_iters 10000 \
  --n_tree_iters 100 \
  --sample_every 10 \
  --infer_nt 1 \
  --out_dir ./output \
  --seed 42
```

## Cleaning build files

To remove compiled objects and the executable:
```bash
make clean
```

## License

[Add your license information here]

## Citation

[Add citation information here]