# Cas12Fold: Sequence-Enhanced and Structure-Informed Tertiary Structure Prediction for Cas12 Proteins

---
## Setup and installation
### Download cas12fold
```bash
git clone https://github.com/CrazyHsu/cas12fold.git
cd cas12fold
```

### Install AlphaFold2 (non-docker version)
#### **Install miniconda**

``` bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && bash Miniconda3-latest-Linux-x86_64.sh
```

#### **Create a new conda environment and update**

``` bash
conda create --name af23 python==3.8
conda update -n base conda
```

#### **Activate conda environment**

``` bash
conda activate af23
```

#### **Install dependencies**

- Change `cudatoolkit==11.2.2` version if it is not supported in your system

``` bash
conda install -y -c conda-forge openmm==7.5.1 cudatoolkit==11.2.2 pdbfixer
conda install -y -c bioconda hmmer hhsuite==3.3.0 kalign2
```

- Change `jaxlib==0.3.25+cuda11.cudnn805` version if this is not supported in your system

``` bash
pip install absl-py==1.0.0 biopython==1.79 chex==0.0.7 dm-haiku==0.0.9 dm-tree==0.1.6 immutabledict==2.0.0 jax==0.3.25 ml-collections==0.1.0 numpy==1.21.6 pandas==1.3.4 protobuf==3.20.1 scipy==1.7.0 tensorflow-cpu==2.9.0

pip install --upgrade --no-cache-dir jax==0.3.25 jaxlib==0.3.25+cuda11.cudnn805 -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
```

### **Download chemical properties to the common folder**

``` bash
wget -q -P alphafold/common/ https://git.scicore.unibas.ch/schwede/openstructure/-/raw/7102c63615b64735c4941278d92b554ec94415f8/modules/mol/alg/src/stereo_chemical_props.txt
```

### **Apply OpenMM patch**

``` bash
# Replace $MULTICOM3_INSTALL_DIR with your MULTICOM3 installation directory

cd ~/anaconda3/envs/af23/lib/python3.8/site-packages/ && patch -p0 < $MULTICOM3_INSTALL_DIR/tools/alphafold-v2.3.2/docker/openmm.patch

# or

cd ~/miniconda3/envs/af23/lib/python3.8/site-packages/ && patch -p0 < $MULTICOM3_INSTALL_DIR/tools/alphafold-v2.3.2/docker/openmm.patch
```

### **Download required databases**
```bash
bash download_dbs.sh -d </home/crazyyhsu/alphafold_data>
```

## **Genetic databases used by Cas12Fold*

Assume the following databases have been installed as a part of the AlphaFold2/AlphaFold-Multimer installation
*   [BFD](https://bfd.mmseqs.com/),
*   [MGnify](https://www.ebi.ac.uk/metagenomics/),
*   [PDB70](http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/),
*   [PDB](https://www.rcsb.org/) (structures in the mmCIF format),
*   [PDB seqres](https://www.rcsb.org/)
*   [UniRef30](https://uniclust.mmseqs.com/),
*   [UniProt](https://www.uniprot.org/uniprot/),
*   [UniRef90](https://www.uniprot.org/help/uniref).

Additional databases will be used for the Cas12Fold prediction process which can be downloaded in [zenodo]():
*   **Cas12FoldDB sequence**: A candidate database that mines Cas12 sequences from the JGI metagenome database and integrates publicly reported Cas12 protein sequences.
*   **RCSB FoldSeek database**: A FoldSeek database built with RCSB PDBs (downloaded in 24/4/2024, 218500 sequences).
*   **Cas12FoldDB template hhdatabase**: A hhdatabase built with the result by aligning Cas12FoldDB to RCSB FoldSeek database using FoldSeek. 
*   **Cas12FoldDB-based AFDB50 FoldSeek database**: A FoldSeek database built with a subset of AFDB50 by aligning the structures of Cas12FoldDB to AFDB50 FoldSeek database using FoldSeek.
*   **Cas12FoldDB-based AFDB50 FoldComp database**: A FoldComp database built with PDBs which used to build Cas12FoldDB-based AFDB50 FoldSeek database



