### 激活环境
```bash
module load oneapi24u1
conda activate /share/apps/miniconda3/envs/af23
```

### 到指定目录下面，然后进行参数配置
```bash
cd /share/share_data/xufeng/cas12fold

python configure.py \
  --template_option_file template_db_options.cas12fold \
  --conda_env_dir /share/apps/miniconda3/envs/af23 \
  --install_dir /share/share_data/xufeng/cas12fold/ \
  --af2db_dir /share/share_data/xufeng/cas12fold/alphafold2_db \
  --out_option_file db_options
```

### 运行AlphaFold2
```bash

fasta_file="examples/8DC2.fa"
python run_cas12fold_controller.py \
  --fasta_path $fasta_file \
  --gpu_device 0 \
  --option_file db_options \
  --output_dir af2_cas12_pred \
  --run_af2_wt \
  --af2_use_precomputed_msas
```

### 运行Cas12Fold
```bash

fasta_file="examples/8DC2.fa"
python run_cas12fold_controller.py \
  --fasta_path $fasta_file \
  --gpu_device 0 \
  --option_file db_options \
  --output_dir cas12fold_cas12_pred \
  --run_cas12fold \
  --af2_use_precomputed_msas
```

### 以Cas12Fold运行的结果进行结构微调
```bash
name="8DC2"
fasta_file="examples/$name.fa"
msas_dir=$(realpath "cas12fold_cas12_pred/$name/msas")
pdb_path=$(realpath "cas12fold_cas12_pred/$name/ranked_0.pdb")
pkl_path=$(realpath "cas12fold_cas12_pred/$name/result_model_1_pred_0.pkl")

python run_cas12fold_refiner_controller.py \
  --fasta_path $fasta_file \
  --msas_dir $msas_dir \
  --option_file db_options \
  --output_dir cas12fold_refine_cas12_pred \
  --pdb_name $name \
  --pdb_path $pdb_path \
  --pkl_path $pkl_path \
  --gpu_device 0 \
  --result_overwrite
```

### 打印预测pdb的平均pLDDT值
```bash
name="8DC2"
python cal_mean_plddt.py alphafold2_cas12_pred/$name/ranked_0.pdb
python cal_mean_plddt.py cas12fold_cas12_pred/$name/ranked_0.pdb
python cal_mean_plddt.py cas12fold_refine_cas12_pred/$name/iteration1/alphafold/ranked_0.pdb
# af2_cas12_pred/8DC2/ranked_0.pdb  57.29467788264259
# cas12fold_cas12_pred/8DC2/ranked_0.pdb    90.36121050471893
# cas12fold_refine_cas12_pred/8DC2/iteration1/alphafold/ranked_0.pdb        92.8979286007386
```


