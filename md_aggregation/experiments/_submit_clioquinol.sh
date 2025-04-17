aslurm -cn haicore_4gpu \
    cmd "conda activate md_aggregation ; export CUDA_VISIBLE_DEVICES=0 ; python aggregation_simulation__clioquinol.py" \
    cmd "conda activate md_aggregation ; export CUDA_VISIBLE_DEVICES=1 ; python aggregation_simulation__clioquinol.py" \
    cmd "conda activate md_aggregation ; export CUDA_VISIBLE_DEVICES=2 ; python aggregation_simulation__clioquinol.py" \
    cmd "conda activate md_aggregation ; export CUDA_VISIBLE_DEVICES=3 ; python aggregation_simulation__clioquinol.py"

aslurm -cn haicore_4gpu \
    cmd "conda activate md_aggregation ; export CUDA_VISIBLE_DEVICES=0 ; python aggregation_simulation__clioquinol_methyl.py" \
    cmd "conda activate md_aggregation ; export CUDA_VISIBLE_DEVICES=1 ; python aggregation_simulation__clioquinol_methyl.py" \
    cmd "conda activate md_aggregation ; export CUDA_VISIBLE_DEVICES=2 ; python aggregation_simulation__clioquinol_methyl.py" \
    cmd "conda activate md_aggregation ; export CUDA_VISIBLE_DEVICES=3 ; python aggregation_simulation__clioquinol_methyl.py"

aslurm -cn haicore_4gpu \
    cmd "conda activate md_aggregation ; export CUDA_VISIBLE_DEVICES=0 ; python aggregation_simulation__clioquinol_without.py" \
    cmd "conda activate md_aggregation ; export CUDA_VISIBLE_DEVICES=1 ; python aggregation_simulation__clioquinol_without.py" \
    cmd "conda activate md_aggregation ; export CUDA_VISIBLE_DEVICES=2 ; python aggregation_simulation__clioquinol_without.py" \
    cmd "conda activate md_aggregation ; export CUDA_VISIBLE_DEVICES=3 ; python aggregation_simulation__clioquinol_without.py"