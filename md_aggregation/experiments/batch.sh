aslurm -cn haicore_4gpu \
    cmd "conda activate md_aggregation ; export CUDA_VISIBLE_DEVICES=0 ; python aggregation_simulation__nesabi28_org.py" \
    cmd "conda activate md_aggregation ; export CUDA_VISIBLE_DEVICES=1 ; python aggregation_simulation__nesabi28_org.py" \
    cmd "conda activate md_aggregation ; export CUDA_VISIBLE_DEVICES=2 ; python aggregation_simulation__nesabi28_org.py" \
    cmd "conda activate md_aggregation ; export CUDA_VISIBLE_DEVICES=3 ; python aggregation_simulation__nesabi28_cf1.py"
