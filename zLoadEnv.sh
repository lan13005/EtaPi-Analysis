conda deactivate
source ~/setUp_websubmit_configs.sh
module unload gcc
module load gcc/4.9.2

echo "gcc: " $(which gcc)
echo "python: " $(which python)