# inPhynite

A C++ tool for Bayesian phylogenetic inference under the infinite-sites mutation model.

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
  --nt_init 1.0 \
  --nt_log_slope 0.0 \
  --n 40 \
  --rooted 0 \
  --n_global_iters 10000 \
  --n_tree_iters 10000 \
  --sample_every 10 \
  --infer_nt 0 \
  --out_dir ./output \
  --seed 12345
```

## Cleaning build files

To remove compiled objects and the executable:
```bash
make clean
```

## License

Copyright 2025 Ivan Specht

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## Citation

Specht, I. & Palacios, J.A. Efficient Bayesian Phylogenetics under the Infinite-Sites Model with inPhynite. arXiv, 2025.