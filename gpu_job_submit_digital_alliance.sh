#!/bin/bash
# ---------------------------------------------------------------------
# SLURM script for a gpu job on a Compute Canada cluster. 
# ---------------------------------------------------------------------
#SBATCH --account=dummy_group
#SBATCH --gres=gpu:1        # Request GPU "generic resources"
#SBATCH --cpus-per-task=16  # Cores proportional to GPUs: 6 on Cedar, 16 on Graham.
#SBATCH --mem=63500M        # Memory proportional to GPUs: 31500 Cedar, 63500 Graham.
#SBATCH --time=24:00:00
#SBATCH --job-name=test_job
#SBATCH --output=V1_human_lymph_log-%j.out # See this log to confirm if your model is running. 
# ---------------------------------------------------------------------
echo "Current working directory: `pwd`"
echo "Starting run at: `date`"
# ---------------------------------------------------------------------
echo ""
echo "Job Array ID / Job ID: $SLURM_ARRAY_JOB_ID / $SLURM_JOB_ID"
echo "This is job $SLURM_ARRAY_TASK_ID out of $SLURM_ARRAY_TASK_COUNT jobs."
echo ""
# ---------------------------------------------------------------------
# Activate your virtual environment
source /home/[user_name]/ENV/bin/activate

# Load necessary modules
module load python/3.10

echo 'Available gpu: '
nvidia-smi

# Change current working directory to the downloaded NEST repository
cd /project/[group_name]/[user_name]/NEST

# Run your python script with parameters
echo 'Running NEST'
cellnest run --data_name='V1_Human_Lymph_Node_spatial' --num_epoch 80000 --manual_seed='yes' --seed=1 --model_name 'NEST_V1_Human_Lymph_Node_spatial' --run_id=1 
# If the process starts successfully, you will see the progress log here: V1_human_lymph_log-%j.out
# Optionally, you can also redirect the output of the above statement to some other log file by appending this to 
# the above line: > output_human_lymph_node_run1.log 
#echo 'To track the training progress see the output_human_lymph_node_*.log under the working directory'

# If you want to start multiple training in parallel using the same GPU, then use following statements:
# nohup cellnest run --data_name='V1_Human_Lymph_Node_spatial' --num_epoch 80000 --manual_seed='yes' --seed=1 --model_name='NEST_V1_Human_Lymph_Node_spatial' --run_id=1 > output_human_lymph_node_run1.log &
# nohup cellnest run  --data_name='V1_Human_Lymph_Node_spatial' --num_epoch 80000 --manual_seed='yes' --seed=2 --model_name='NEST_V1_Human_Lymph_Node_spatial' --run_id=2 > output_human_lymph_node_run2.log &
# nohup cellnest run  --data_name='V1_Human_Lymph_Node_spatial' --num_epoch 80000 --manual_seed='yes' --seed=3 --model_name='NEST_V1_Human_Lymph_Node_spatial' --run_id=3 > output_human_lymph_node_run3.log &
# nohup cellnest run  --data_name='V1_Human_Lymph_Node_spatial' --num_epoch 80000 --manual_seed='yes' --seed=4 --model_name='NEST_V1_Human_Lymph_Node_spatial' --run_id=4 > output_human_lymph_node_run4.log &
# nohup cellnest run  --data_name='V1_Human_Lymph_Node_spatial' --num_epoch 80000 --manual_seed='yes' --seed=5 --model_name='NEST_V1_Human_Lymph_Node_spatial' --run_id=5 > output_human_lymph_node_run5.log &
# In this case, please note that the above five processes are running in the background. But you cannot exit the current shell 
# script as that will revoke those processes as well. Therefore, you have to put current shell script into sleep as below:
# sleep 86400
# The sleep time is an approximate time required by the model to finish training. It is encouraged to see the progress of 
# the running processes using the log files. If your log files are saying that all the instances have finished training, then you should
# kill this sleeping shell script using the corresponding JOBID. The JOBID can be found using this command: squeue $USR
# ---------------------------------------------------------------------


echo "Job finished with exit code $? at: `date`"
