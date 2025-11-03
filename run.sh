# Fixed nt at 1, 2, 4, 8
#./inPhynite --nt_init 1.0 --nt_log_slope 0.0 --n 40 --rooted 0 --n_global_iters 10000 --n_tree_iters 10000 --sample_every 10 --infer_nt 0 --out_dir /Users/ispecht/Desktop/PhyloSMC-analyses/runs/PhyloSMC/fixed_N_eff/Ne_1

#./inPhynite --nt_init 2.0 --nt_log_slope 0.0 --n 40 --rooted 0 --n_global_iters 10000 --n_tree_iters 10000 --sample_every 10 --infer_nt 0 --out_dir /Users/ispecht/Desktop/PhyloSMC-analyses/runs/PhyloSMC/fixed_N_eff/Ne_2

#./inPhynite --nt_init 4.0 --nt_log_slope 0.0 --n 40 --rooted 0 --n_global_iters 10000 --n_tree_iters 10000 --sample_every 10 --infer_nt 0 --out_dir /Users/ispecht/Desktop/PhyloSMC-analyses/runs/PhyloSMC/fixed_N_eff/Ne_4

#./inPhynite --nt_init 8.0 --nt_log_slope 0.0 --n 40 --rooted 0 --n_global_iters 10000 --n_tree_iters 10000 --sample_every 10 --infer_nt 0 --out_dir /Users/ispecht/Desktop/PhyloSMC-analyses/runs/PhyloSMC/fixed_N_eff/Ne_8

# Inferred nt, start 128, log slope -1
#./inPhynite --nt_init 128.0 --nt_log_slope -1.0 --n 40 --rooted 0 --n_global_iters 100000 --n_tree_iters 10000 --sample_every 100 --infer_nt 1 --out_dir /Users/ispecht/Desktop/PhyloSMC-analyses/runs/PhyloSMC/inferred_N_eff/drift_-1

# Inferred nt, constant 32
./inPhynite --nt_init 32.0 --nt_log_slope 0.0 --n 40 --rooted 0 --n_global_iters 100000 --n_tree_iters 10000 --sample_every 100 --infer_nt 1 --out_dir /Users/ispecht/Desktop/PhyloSMC-analyses/runs/PhyloSMC/inferred_N_eff/drift_0

# start 1.5, log slope 1
#./inPhynite --nt_init 1.5 --nt_log_slope 1.0 --n 40 --rooted 0 --n_global_iters 100000 --n_tree_iters 10000 --sample_every 100 --infer_nt 1 --out_dir /Users/ispecht/Desktop/PhyloSMC-analyses/runs/PhyloSMC/inferred_N_eff/drift_1

# Real data
#./inPhynite --rooted 1 --n_global_iters 100000 --n_tree_iters 10000 --sample_every 100 --infer_nt 1 --in_file /Users/ispecht/Desktop/PhyloSMC-analyses/data/YRI.fasta --out_dir /Users/ispecht/Desktop/PhyloSMC-analyses/runs/PhyloSMC/real_data
