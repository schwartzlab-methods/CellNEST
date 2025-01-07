# Python package installation

On Digital Alliance, you have to load the Python module first as follows (this can be skipped if you are working on your local machine):
```
module load python/3.10
```

Then, run the following commands to create a virtual environment for running CellNEST and activate it before installing Python libraries to setup the environment:
```
virtualenv --no-download /project/[group_name]/[user_name]/cellnest_venv
source /project/[group_name]/[user_name]/cellnest_venv/bin/activate
```
Please replace [group_name] and [user_name] with your assigned group and name. 

Upgrade pip3:
```
pip3 install --no-index --upgrade pip
```

After that, use the requirements.txt to install the required Python libraries as follows:
```
pip3 install -r requirements.txt
```

The above command should setup the environment properly. If the setup fails due to unmatched CUDA in your system, then please perform following steps to setup Torch with CUDA manually.

## Install Torch with CUDA support for GPU usage
See which Torch and CUDA versions are available:
```
avail_wheels "torch*"
module load cuda [press tab to see the available versions]
```

Next, based on the Pytorch and CUDA version, run with [wheels links](https://pytorch.org/get-started/previous-versions/
) for installing. Here we are using pytorch 1.13.1 with CUDA 11.7:

```
pip3 install torch-scatter -f https://data.pyg.org/whl/torch-1.13.1+cu117.html
pip3 install torch-sparse -f https://data.pyg.org/whl/torch-1.13.1+cu117.html
pip3 install torch-cluster -f https://data.pyg.org/whl/torch-1.13.1+cu117.html
pip3 install torch-geometric -f https://data.pyg.org/whl/torch-1.13.1+cu117.html
```

After the environment is setup, you can exit from your virtual environment using following command:
```
deactivate
```

Please remember to activate the 'cellnest_venv' using 'source' command everytime you run CellNEST. 

# Submit GPU jobs

There is a sample shell script written for submitting GPU job as follows:
```
sbatch gpu_job_submit_compute_canada.sh
```
Please see the 'gpu_job_submit_digital_alliance.sh' in the CellNEST repository to understand the requested GPU resources. This will execute just one run of CellNEST. If enough GPU is available, multiple CellNEST runs can be executed in parallel by replacing the contents of line 35 with the following:
```
nohup cellnest run --data_name='V1_Human_Lymph_Node_spatial' --num_epoch 80000 --seed=1 --model_name 'CellNEST_V1_Human_Lymph_Node_spatial' --run_id=1 > output_human_lymph_node_run1.log &
nohup cellnest run  --data_name='V1_Human_Lymph_Node_spatial' --num_epoch 80000 --seed=2 --model_name 'CellNEST_V1_Human_Lymph_Node_spatial' --run_id=2 > output_human_lymph_node_run2.log &
nohup cellnest run  --data_name='V1_Human_Lymph_Node_spatial' --num_epoch 80000 --seed=3 --model_name 'CellNEST_V1_Human_Lymph_Node_spatial' --run_id=3 > output_human_lymph_node_run3.log &
nohup cellnest run  --data_name='V1_Human_Lymph_Node_spatial' --num_epoch 80000 --seed=4 --model_name 'CellNEST_V1_Human_Lymph_Node_spatial' --run_id=4 > output_human_lymph_node_run4.log &
nohup cellnest run  --data_name='V1_Human_Lymph_Node_spatial' --num_epoch 80000 --seed=5 --model_name 'CellNEST_V1_Human_Lymph_Node_spatial' --run_id=5 > output_human_lymph_node_run5.log &
```




