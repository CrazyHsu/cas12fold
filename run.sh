module load oneapi24u1
conda activate /share/apps/miniconda3/envs/af23


python configure.py --template_option_file template_db_options \
  --conda_env_dir /share/apps/miniconda3/envs/af23 \
  --tools_dir /share/home/xufeng/Projects/test_cas12fold_refine/refine_cas12/cas12fold_20241014/tools \
  --afdb_dir /share/home/xufeng/Data/alphafold2_db \
  --af2_dir /share/home/xufeng/Projects/test_cas12fold_refine/refine_cas12/cas12fold_20241014/alphafold \
  --out_option_file db_options


python configure.py --template_option_file template_db_options --conda_env_dir /share/apps/miniconda3/envs/af23 --tools_dir /path/to/tools --afdb_dir /share/home/xufeng/Software/alphafold-2.3.2/ --af2_dir /share/home/xufeng/Software/alphafold-2.3.2/alphafold --out_option_file db_options

fasta_file="/share/home/xufeng/Projects/af2_optimizm/cas12fold_test/batch_test_20240906/cas12fold_pred/cas12_1.dbtime_20990101/fastas/cas12_1.fa"
msas_dir="/share/home/xufeng/Projects/af2_optimizm/cas12fold_test/batch_test_20240906/cas12fold_pred/cas12_1.dbtime_20990101/wt/cas12_1/msas/"
pdb_path="/share/home/xufeng/Projects/af2_optimizm/cas12fold_test/batch_test_20240906/cas12fold_pred/cas12_1.dbtime_20990101/wt/cas12_1/relaxed_model_3_ptm_pred_0.pdb"
pdb_name="cas12_1"
pkl_path="/share/home/xufeng/Projects/af2_optimizm/cas12fold_test/batch_test_20240906/cas12fold_pred/cas12_1.dbtime_20990101/wt/cas12_1/result_model_3_ptm_pred_0.pkl"
python run_cas12fold.py --fasta_path $fasta_file --msas_dir $msas_dir --option_file db_options --output_dir cas12_refine_test --pdb_name $pdb_name --pdb_path $pdb_path --pkl_path $pkl_path --run_refinement

python run_cas12fold.py --fasta_path $fasta_file --gpu_device 1 --option_file db_options --output_dir cas12_refine_test


fasta_file="/share/home/xufeng/Projects/af2_optimizm/cas12fold_test/batch_test_20240906/cas12fold_pred/cas12_1.dbtime_20990101/fastas/cas12_1.fa"
pdb_path="/share/home/xufeng/Projects/test_cas12fold_refine/refine_cas12/cas12fold_20241014/cas12_refine_test1/cas12_1/ranked_0.pdb"
msas_dir="/share/home/xufeng/Projects/test_cas12fold_refine/refine_cas12/cas12fold_20241014/cas12_refine_test1/cas12_1/msas"
pkl_path="/share/home/xufeng/Projects/test_cas12fold_refine/refine_cas12/cas12fold_20241014/cas12_refine_test1/cas12_1/result_model_1_pred_0.pkl"
pdb_name="cas12_1"
python run_cas12fold_refiner_controller.py --fasta_path $fasta_file --msas_dir $msas_dir --option_file db_options --output_dir cas12_refine_test2 --pdb_name $pdb_name --pdb_path $pdb_path --pkl_path $pkl_path


fasta_file="/share/home/xufeng/Projects/af2_optimizm/cas12fold_test/batch_test_20240906/cas12fold_pred/cas12_2.dbtime_20990101/fastas/cas12_2.fa"
pdb_path="/share/home/xufeng/Projects/test_cas12fold_refine/refine_cas12/cas12fold_20241014/cas12_refine_test1/cas12_2/ranked_0.pdb"
msas_dir="/share/home/xufeng/Projects/test_cas12fold_refine/refine_cas12/cas12fold_20241014/cas12_refine_test1/cas12_2/msas"
pkl_path="/share/home/xufeng/Projects/test_cas12fold_refine/refine_cas12/cas12fold_20241014/cas12_refine_test1/cas12_2/result_model_1_pred_0.pkl"
pdb_name="cas12_2"
python run_cas12fold_refiner_controller.py --fasta_path $fasta_file --msas_dir $msas_dir --option_file db_options --output_dir cas12_refine_test2 --pdb_name $pdb_name --pdb_path $pdb_path --pkl_path $pkl_path --gpu_device 3 --result_overwrite


fasta_file="/share/home/xufeng/Projects/af2_optimizm/cas12fold_test/batch_test_20240906/cas12fold_pred/cas12_2.dbtime_20990101/fastas/cas12_2.fa"
pdb_path="/share/home/xufeng/Projects/test_cas12fold_refine/refine_cas12/cas12fold_20241020/cas12_refine_test1/cas12_2/ranked_0.pdb"
msas_dir="/share/home/xufeng/Projects/test_cas12fold_refine/refine_cas12/cas12fold_20241020/cas12_refine_test1/cas12_2/msas"
pkl_path="/share/home/xufeng/Projects/test_cas12fold_refine/refine_cas12/cas12fold_20241020/cas12_refine_test1/cas12_2/result_model_1_pred_0.pkl"
pdb_name="cas12_2"
python run_cas12fold_refiner_controller.py --fasta_path $fasta_file --msas_dir $msas_dir --option_file db_options --output_dir cas12_refine_test2 --pdb_name $pdb_name --pdb_path $pdb_path --pkl_path $pkl_path --gpu_device 3

for i in cas12_74 cas12_103 cas12_104; do
  fasta_file="/share/home/xufeng/Projects/af2_optimizm/cas12fold_test/batch_test_20240906/cas12fold_pred/$i.dbtime_20990101/fastas/$i.fa"
  msas_dir="/share/home/xufeng/Projects/af2_optimizm/cas12fold_test/batch_test_20240906/cas12fold_pred/$i.dbtime_20990101/merged_wt_cas12folddb/$i/msas/"
  pdb_path="/share/home/xufeng/Projects/af2_optimizm/cas12fold_test/batch_test_20240906/cas12fold_pred/$i.dbtime_20990101/merged_wt_cas12folddb/$i/ranked_0.pdb"
  pdb_name=$i
  pkl_path="/share/home/xufeng/Projects/af2_optimizm/cas12fold_test/batch_test_20240906/cas12fold_pred/$i.dbtime_20990101/merged_wt_cas12folddb/$i/result_model_1_ptm_pred_0.pkl"
  echo -e "python cas12_refiner.py --fasta_path $fasta_file --msas_dir $msas_dir --option_file db_options --output_dir cas12_refine_test --pdb_name $pdb_name --pdb_path $pdb_path --pkl_path $pkl_path --run_refinement 1>test_logs/$i.log 2>&1 &"
done > test_cas12refiner.sh

for i in cas12_95 cas12_9 cas12_13; do
  fasta_file="/share/home/xufeng/Projects/af2_optimizm/cas12fold_test/batch_test_20240906/cas12fold_pred/$i.dbtime_20990101/fastas/$i.fa"
  msas_dir="/share/home/xufeng/Projects/af2_optimizm/cas12fold_test/batch_test_20240906/cas12fold_pred/$i.dbtime_20990101/merged_wt_cas12folddb/$i/msas/"
  pdb_path="/share/home/xufeng/Projects/af2_optimizm/cas12fold_test/batch_test_20240906/cas12fold_pred/$i.dbtime_20990101/merged_wt_cas12folddb/$i/ranked_0.pdb"
  pdb_name=$i
  pkl_path="/share/home/xufeng/Projects/af2_optimizm/cas12fold_test/batch_test_20240906/cas12fold_pred/$i.dbtime_20990101/merged_wt_cas12folddb/$i/result_model_1_ptm_pred_0.pkl"
  echo -e "python cas12_refiner.py --fasta_path $fasta_file --msas_dir $msas_dir --option_file db_options --output_dir cas12_refine_test --pdb_name $pdb_name --pdb_path $pdb_path --pkl_path $pkl_path --run_refinement 1>test_logs/$i.log 2>&1 &"
done > test_cas12refiner_2.sh

ls /share/home/xufeng/Projects/af2_optimizm/cas12fold_test/batch_test_20240906/low_plddt_cas12/*.fa | rush "basename {}" | sed 's/.fa//g' | grep -v -f <(ls /share/home/xufeng/Projects/test_cas12fold_refine/refine_cas12/pythonProject/cas12_refine_test/) - | head -8 | while read i;do
  fasta_file="/share/home/xufeng/Projects/af2_optimizm/cas12fold_test/batch_test_20240906/cas12fold_pred/$i.dbtime_20990101/fastas/$i.fa"
  msas_dir="/share/home/xufeng/Projects/af2_optimizm/cas12fold_test/batch_test_20240906/cas12fold_pred/$i.dbtime_20990101/merged_wt_cas12folddb/$i/msas/"
  pdb_path="/share/home/xufeng/Projects/af2_optimizm/cas12fold_test/batch_test_20240906/cas12fold_pred/$i.dbtime_20990101/merged_wt_cas12folddb/$i/ranked_0.pdb"
  pdb_name=$i
  pkl_path="/share/home/xufeng/Projects/af2_optimizm/cas12fold_test/batch_test_20240906/cas12fold_pred/$i.dbtime_20990101/merged_wt_cas12folddb/$i/result_model_1_ptm_pred_0.pkl"
  echo -e "python cas12_refiner.py --fasta_path $fasta_file --msas_dir $msas_dir --option_file db_options --output_dir cas12_refine_test --pdb_name $pdb_name --pdb_path $pdb_path --pkl_path $pkl_path --run_refinement 1>test_logs/$i.log 2>&1 &"
done > test_cas12refiner_3.sh

ls /share/home/xufeng/Projects/af2_optimizm/cas12fold_test/batch_test_20240906/low_plddt_cas12/*.fa | rush "basename {}" | sed 's/.fa//g' | grep -v -f <(ls /share/home/xufeng/Projects/test_cas12fold_refine/refine_cas12/pythonProject/cas12_refine_test/) - | tail -8 | while read i;do
  fasta_file="/share/home/xufeng/Projects/af2_optimizm/cas12fold_test/batch_test_20240906/cas12fold_pred/$i.dbtime_20990101/fastas/$i.fa"
  msas_dir="/share/home/xufeng/Projects/af2_optimizm/cas12fold_test/batch_test_20240906/cas12fold_pred/$i.dbtime_20990101/merged_wt_cas12folddb/$i/msas/"
  pdb_path="/share/home/xufeng/Projects/af2_optimizm/cas12fold_test/batch_test_20240906/cas12fold_pred/$i.dbtime_20990101/merged_wt_cas12folddb/$i/ranked_0.pdb"
  pdb_name=$i
  pkl_path="/share/home/xufeng/Projects/af2_optimizm/cas12fold_test/batch_test_20240906/cas12fold_pred/$i.dbtime_20990101/merged_wt_cas12folddb/$i/result_model_1_ptm_pred_0.pkl"
  echo -e "python cas12_refiner.py --fasta_path $fasta_file --msas_dir $msas_dir --option_file db_options --output_dir cas12_refine_test --pdb_name $pdb_name --pdb_path $pdb_path --pkl_path $pkl_path --run_refinement 1>test_logs/$i.log 2>&1 &"
done > test_cas12refiner_4.sh


ls /share/home/xufeng/Projects/af2_optimizm/cas12fold_test/batch_test_20240903/rcsb_cas12_fa/*.fa | rush "basename {}" | sed 's/.fa//g' | head | while read i;do
  fasta_file="/share/home/xufeng/Projects/af2_optimizm/cas12fold_test/batch_test_20240903/cas12fold_pred/$i.dbtime_20990101/fastas/${i}.single.fa"
  msas_dir="/share/home/xufeng/Projects/af2_optimizm/cas12fold_test/batch_test_20240903/cas12fold_pred/$i.dbtime_20990101/merged_wt_cas12folddb/${i}.single/msas/"
  pdb_path="/share/home/xufeng/Projects/af2_optimizm/cas12fold_test/batch_test_20240903/cas12fold_pred/$i.dbtime_20990101/merged_wt_cas12folddb/${i}.single/ranked_0.pdb"
  pdb_name=$i
  pkl_path="/share/home/xufeng/Projects/af2_optimizm/cas12fold_test/batch_test_20240903/cas12fold_pred/$i.dbtime_20990101/merged_wt_cas12folddb/${i}.single/result_model_1_ptm_pred_0.pkl"
  echo -e "python cas12_refiner.py --fasta_path $fasta_file --msas_dir $msas_dir --option_file db_options --output_dir cas12_refine_test --pdb_name $pdb_name --pdb_path $pdb_path --pkl_path $pkl_path --run_refinement 1>test_logs/$i.log 2>&1 &"
done > test_cas12refiner_5.sh

j=0
ls /share/home/xufeng/Projects/af2_optimizm/cas12fold_test/batch_test_20240906/low_plddt_cas12/*.fa | rush "basename {}" | grep -v -f <(ls /share/home/xufeng/Projects/test_cas12fold_refine/refine_cas12/pythonProject/cas12_refine_test/ | awk '{print $1".fa"}') - | sed 's/.fa//g' | head -42 | while read i;do
  gpu=$(($j%4))
  if [[ $gpu == 1 ]]; then
    j=$(($j+1))
    gpu=$(($j%4))
  fi
  j=$(($j+1))
  fasta_file="/share/home/xufeng/Projects/af2_optimizm/cas12fold_test/batch_test_20240906/cas12fold_pred/$i.dbtime_20990101/fastas/$i.fa"
  msas_dir="/share/home/xufeng/Projects/af2_optimizm/cas12fold_test/batch_test_20240906/cas12fold_pred/$i.dbtime_20990101/merged_wt_cas12folddb/$i/msas/"
  pdb_path="/share/home/xufeng/Projects/af2_optimizm/cas12fold_test/batch_test_20240906/cas12fold_pred/$i.dbtime_20990101/merged_wt_cas12folddb/$i/ranked_0.pdb"
  pdb_name=$i
  pkl_path="/share/home/xufeng/Projects/af2_optimizm/cas12fold_test/batch_test_20240906/cas12fold_pred/$i.dbtime_20990101/merged_wt_cas12folddb/$i/result_model_1_ptm_pred_0.pkl"
  echo -e "python cas12_refiner.py --fasta_path $fasta_file --msas_dir $msas_dir --option_file db_options --output_dir cas12_refine_test --pdb_name $pdb_name --pdb_path $pdb_path --pkl_path $pkl_path --gpu_device $gpu --run_refinement 1>test_logs/$i.log 2>&1"
done > run_105_low_plddt_cas12.head42.sh

j=0
ls /share/home/xufeng/Projects/af2_optimizm/cas12fold_test/batch_test_20240906/low_plddt_cas12/*.fa | rush "basename {}" | grep -v -f <(ls /share/home/xufeng/Projects/test_cas12fold_refine/refine_cas12/pythonProject/cas12_refine_test/ | awk '{print $1".fa"}') - | sed 's/.fa//g' | tail -41 | while read i;do
  gpu=$(($j%4))
  if [[ $gpu == 0 || $gpu == 2 ]]; then
    j=$(($j+1))
    gpu=$(($j%4))
  fi
  j=$(($j+1))
  fasta_file="/share/home/xufeng/Projects/af2_optimizm/cas12fold_test/batch_test_20240906/cas12fold_pred/$i.dbtime_20990101/fastas/$i.fa"
  msas_dir="/share/home/xufeng/Projects/af2_optimizm/cas12fold_test/batch_test_20240906/cas12fold_pred/$i.dbtime_20990101/merged_wt_cas12folddb/$i/msas/"
  pdb_path="/share/home/xufeng/Projects/af2_optimizm/cas12fold_test/batch_test_20240906/cas12fold_pred/$i.dbtime_20990101/merged_wt_cas12folddb/$i/ranked_0.pdb"
  pdb_name=$i
  pkl_path="/share/home/xufeng/Projects/af2_optimizm/cas12fold_test/batch_test_20240906/cas12fold_pred/$i.dbtime_20990101/merged_wt_cas12folddb/$i/result_model_1_ptm_pred_0.pkl"
  echo -e "python cas12_refiner.py --fasta_path $fasta_file --msas_dir $msas_dir --option_file db_options --output_dir cas12_refine_test --pdb_name $pdb_name --pdb_path $pdb_path --pkl_path $pkl_path --gpu_device $gpu --run_refinement 1>test_logs/$i.log 2>&1"
done > run_105_low_plddt_cas12.tail41.sh

sbatch -J cas12_refiner --ntasks-per-node=6 ~/sbatch_submit_templates/sbatch_submit_rush.head.gpu01.sh run_105_low_plddt_cas12.head42.sh 55 6
sbatch -J cas12_refiner --ntasks-per-node=6 ~/sbatch_submit_templates/sbatch_submit_rush.tail.gpu02.sh run_105_low_plddt_cas12.tail41.sh 55 6


j=0
ls /share/home/xufeng/Projects/af2_optimizm/cas12fold_test/batch_test_20240906/low_plddt_cas12/*.fa | rush "basename {}" | grep -v -f <(find /share/home/xufeng/Projects/test_cas12fold_refine/refine_cas12/pythonProject/cas12_refine_test/ -name "ranked_0.pdb" | cut -f 10 -d '/' | grep -v "bak" | awk '{print $1".fa"}') - | sed 's/.fa//g' | while read i;do
  gpu=$(($j%4))
  if [[ $gpu == 1 ]]; then
    j=$(($j+1))
    gpu=$(($j%4))
  fi
  j=$(($j+1))
  fasta_file="/share/home/xufeng/Projects/af2_optimizm/cas12fold_test/batch_test_20240906/cas12fold_pred/$i.dbtime_20990101/fastas/$i.fa"
  msas_dir="/share/home/xufeng/Projects/af2_optimizm/cas12fold_test/batch_test_20240906/cas12fold_pred/$i.dbtime_20990101/merged_wt_cas12folddb/$i/msas/"
  pdb_path="/share/home/xufeng/Projects/af2_optimizm/cas12fold_test/batch_test_20240906/cas12fold_pred/$i.dbtime_20990101/merged_wt_cas12folddb/$i/ranked_0.pdb"
  pdb_name=$i
  pkl_path="/share/home/xufeng/Projects/af2_optimizm/cas12fold_test/batch_test_20240906/cas12fold_pred/$i.dbtime_20990101/merged_wt_cas12folddb/$i/result_model_1_ptm_pred_0.pkl"
  echo -e "python cas12_refiner.py --fasta_path $fasta_file --msas_dir $msas_dir --option_file db_options --output_dir cas12_refine_test --pdb_name $pdb_name --pdb_path $pdb_path --pkl_path $pkl_path --gpu_device $gpu --run_refinement 1>test_logs/$i.log 2>&1"
done > run_105_low_plddt_cas12.remain.sh

sbatch -J cas12_refiner --ntasks-per-node=9 ~/sbatch_submit_templates/sbatch_submit_rush.head.gpu01.sh run_105_low_plddt_cas12.remain.sh 20 9