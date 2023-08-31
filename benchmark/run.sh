# 1. Environment setup.
echo "Please make sure you have the dependencies installed :)"
echo "conda env create -f env.yml / pip install -r requirements.txt"

# 2. Example: train AttE.
CUDA_VISIBLE_DEVICES=0 python main.py --model AttE --dataset whole --valid 5 --neg_sample_size 150 --optimizer Adam --learning_rate 0.001

# 3. More examples below. Please checkout ./configs/ folder for full configurations files in JSON format.
# CUDA_VISIBLE_DEVICES=1 python main.py --model RotE --dataset instance --valid 5 --neg_sample_size 150 --optimizer Adam --learning_rate 0.001
# CUDA_VISIBLE_DEVICES=2 python main.py --model RefE --dataset whole --valid 5 --neg_sample_size 150 --optimizer Adam --learning_rate 0.001
# CUDA_VISIBLE_DEVICES=3 python main.py --model AttH --dataset instance --valid 5 --neg_sample_size 150 --optimizer Adam --learning_rate 0.001 --multi_c
