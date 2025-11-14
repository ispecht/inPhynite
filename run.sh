# Fixed nt
./inPhynite --nt_init 4.0 --nt_log_slope 0.0 --n 32 --rooted 0 --n_global_iters 10000 --n_tree_iters 10000 --sample_every 10 --infer_nt 0 --out_dir /Users/ispecht/Desktop/inPhynite-analyses/runs/inPhynite/fixed_N_eff/run_1

./inPhynite --nt_init 16.0 --nt_log_slope 0.0 --n 32 --rooted 0 --n_global_iters 10000 --n_tree_iters 10000 --sample_every 10 --infer_nt 0 --out_dir /Users/ispecht/Desktop/inPhynite-analyses/runs/inPhynite/fixed_N_eff/run_2

./inPhynite --nt_init 64.0 --nt_log_slope 0.0 --n 32 --rooted 0 --n_global_iters 10000 --n_tree_iters 10000 --sample_every 10 --infer_nt 0 --out_dir /Users/ispecht/Desktop/inPhynite-analyses/runs/inPhynite/fixed_N_eff/run_3

./inPhynite --nt_init 4.0 --nt_log_slope 0.0 --n 128 --rooted 0 --n_global_iters 10000 --n_tree_iters 10000 --sample_every 10 --infer_nt 0 --out_dir /Users/ispecht/Desktop/inPhynite-analyses/runs/inPhynite/fixed_N_eff/run_4

./inPhynite --nt_init 16.0 --nt_log_slope 0.0 --n 128 --rooted 0 --n_global_iters 10000 --n_tree_iters 10000 --sample_every 10 --infer_nt 0 --out_dir /Users/ispecht/Desktop/inPhynite-analyses/runs/inPhynite/fixed_N_eff/run_5

./inPhynite --nt_init 64.0 --nt_log_slope 0.0 --n 128 --rooted 0 --n_global_iters 10000 --n_tree_iters 10000 --sample_every 10 --infer_nt 0 --out_dir /Users/ispecht/Desktop/inPhynite-analyses/runs/inPhynite/fixed_N_eff/run_6

./inPhynite --nt_init 4.0 --nt_log_slope 0.0 --n 512 --rooted 0 --n_global_iters 10000 --n_tree_iters 10000 --sample_every 10 --infer_nt 0 --out_dir /Users/ispecht/Desktop/inPhynite-analyses/runs/inPhynite/fixed_N_eff/run_7

./inPhynite --nt_init 16.0 --nt_log_slope 0.0 --n 512 --rooted 0 --n_global_iters 10000 --n_tree_iters 10000 --sample_every 10 --infer_nt 0 --out_dir /Users/ispecht/Desktop/inPhynite-analyses/runs/inPhynite/fixed_N_eff/run_8

./inPhynite --nt_init 64.0 --nt_log_slope 0.0 --n 512 --rooted 0 --n_global_iters 10000 --n_tree_iters 10000 --sample_every 10 --infer_nt 0 --out_dir /Users/ispecht/Desktop/inPhynite-analyses/runs/inPhynite/fixed_N_eff/run_9

# Inferred nt
./inPhynite --nt_init 128.0 --nt_log_slope -1.0 --n 40 --rooted 0 --n_global_iters 100000 --n_tree_iters 10000 --sample_every 100 --infer_nt 1 --out_dir /Users/ispecht/Desktop/inPhynite-analyses/runs/inPhynite/inferred_N_eff/drift_-1

./inPhynite --nt_init 32.0 --nt_log_slope 0.0 --n 40 --rooted 0 --n_global_iters 100000 --n_tree_iters 10000 --sample_every 100 --infer_nt 1 --out_dir /Users/ispecht/Desktop/inPhynite-analyses/runs/inPhynite/inferred_N_eff/drift_0

./inPhynite --nt_init 1.5 --nt_log_slope 1.0 --n 40 --rooted 0 --n_global_iters 100000 --n_tree_iters 10000 --sample_every 100 --infer_nt 1 --out_dir /Users/ispecht/Desktop/inPhynite-analyses/runs/inPhynite/inferred_N_eff/drift_1

# Real data
./inPhynite --rooted 0 --n_global_iters 100000 --n_tree_iters 10000 --sample_every 100 --infer_nt 1 --in_file /Users/ispecht/Desktop/inPhynite-analyses/data/CEU.fasta --out_dir /Users/ispecht/Desktop/inPhynite-analyses/runs/inPhynite/real_data/CEU

./inPhynite --rooted 0 --n_global_iters 100000 --n_tree_iters 10000 --sample_every 100 --infer_nt 1 --in_file /Users/ispecht/Desktop/inPhynite-analyses/data/FIN.fasta --out_dir /Users/ispecht/Desktop/inPhynite-analyses/runs/inPhynite/real_data/FIN

./inPhynite --rooted 0 --n_global_iters 100000 --n_tree_iters 10000 --sample_every 100 --infer_nt 1 --in_file /Users/ispecht/Desktop/inPhynite-analyses/data/PEL.fasta --out_dir /Users/ispecht/Desktop/inPhynite-analyses/runs/inPhynite/real_data/PEL

./inPhynite --rooted 0 --n_global_iters 100000 --n_tree_iters 10000 --sample_every 100 --infer_nt 1 --in_file /Users/ispecht/Desktop/inPhynite-analyses/data/YRI.fasta --out_dir /Users/ispecht/Desktop/inPhynite-analyses/runs/inPhynite/real_data/YRI
