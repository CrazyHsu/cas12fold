#!/bin/bash
#SBATCH -o /share/home/xufeng/sbatch_submit_templates/sbatch_job_logs/job.%j.out
#SBATCH -e /share/home/xufeng/sbatch_submit_templates/sbatch_job_logs/job.%j.err
#SBATCH -J xufeng
##SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH -w gpu02

#fasta_file="/share/home/xufeng/Projects/af2_optimizm/cas12fold_test/batch_test_20240906/cas12fold_pred/cas12_1.dbtime_20990101/fastas/cas12_1.fa"
fasta_file="/share/home/xufeng/Projects/af2_optimizm/cas12fold_test/batch_test_20240906/cas12fold_pred/cas12_2.dbtime_20990101/fastas/cas12_2.fa"
#python run_cas12fold_controller.py --fasta_path $fasta_file --gpu_device 3 --option_file db_options --output_dir cas12_refine_test --run_af2_wt --af2_use_precomputed_msas
python run_cas12fold_controller.py --fasta_path $fasta_file --gpu_device 3 --option_file db_options --output_dir cas12_refine_test1 --run_cas12fold --af2_use_precomputed_msas
