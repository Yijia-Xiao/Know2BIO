import os
import argparse
import know2bio
from know2bio.config import Trainer, Tester
from know2bio.module.model import *
from know2bio.module.loss import SoftplusLoss, MarginLoss, SigmoidLoss
from know2bio.module.strategy import NegativeSampling
from know2bio.data import TrainDataLoader, TestDataLoader


parser = argparse.ArgumentParser()
parser.add_argument('--model_name', type=str, default='TransE', help='Name of the model to be evaluated.')
parser.add_argument('--split', type=str, default='BMKG', help='The dataset/split to be used.')
parser.add_argument('--data_path', type=str, default='./data/', help='Path to the datasets: GO, Know2BIO\'s ontology/instance/aggregate view, etc.')
parser.add_argument('--ckpt_path', type=str, default='./checkpoint/', help='Path to the model checkpoints storage.')
parser.add_argument('--lr', type=float, default=5e-5, help='Learning rate in training process.')
parser.add_argument('--batch_size', type=int, default=1024, help='Batch size of for modeling training.')
parser.add_argument('--train_epochs', type=int, default=5000, help='Number of training epochs.')
parser.add_argument('--opt_method', type=str, default='adam', help='Optimization method for model training.')


args = parser.parse_args()
model_name = args.model_name
split = args.split
lr = args.lr
batch_size = args.batch_size
train_epochs = args.train_epochs
opt_method = args.opt_method
data_path = args.data_path
ckpt_path = args.ckpt_path


train_dataloader = TrainDataLoader(
	in_path = f"./{data_path}/{split}/", 
	batch_size=batch_size, 
	threads = 32, 
	sampling_mode = "normal", 
	bern_flag = 1, 
	filter_flag = 1, 
	neg_ent = 25,
	neg_rel = 0)
test_dataloader = TestDataLoader(f"./{data_path}/{split}/", "link")


if model_name == 'complex':
    func = ComplEx(
        ent_tot = train_dataloader.get_ent_tot(),
        rel_tot = train_dataloader.get_rel_tot(),
        dim = 200
    )
    model = NegativeSampling(
        model = func, 
        loss = SoftplusLoss(),
        batch_size = train_dataloader.get_batch_size(), 
        regul_rate = 1.0
    )

elif model_name == 'distmult':
    func = DistMult(
        ent_tot = train_dataloader.get_ent_tot(),
        rel_tot = train_dataloader.get_rel_tot(),
        dim = 1024
    )

    model = NegativeSampling(
        model = func, 
        loss = SoftplusLoss(),
        batch_size = train_dataloader.get_batch_size(), 
        regul_rate = 1.0
    )
    opt_method = "adagrad"

elif model_name == 'rotate':
    func = RotatE(
        ent_tot = train_dataloader.get_ent_tot(),
        rel_tot = train_dataloader.get_rel_tot(),
        dim = 1024,
        margin = 6.0,
        epsilon = 2.0,
    )
    model = NegativeSampling(
        model = func, 
        loss = SigmoidLoss(adv_temperature = 2),
        batch_size = train_dataloader.get_batch_size(), 
        regul_rate = 0.0
    )

elif model_name == 'transe':
    func = TransE(
        ent_tot = train_dataloader.get_ent_tot(),
        rel_tot = train_dataloader.get_rel_tot(),
        dim = 1024, 
        p_norm = 1, 
        norm_flag = True)
    model = NegativeSampling(
        model = func, 
        loss = MarginLoss(margin = 5.0),
        batch_size = train_dataloader.get_batch_size()
    )

else:
    raise NotImplementedError


trainer = Trainer(model = model, data_loader = train_dataloader, train_times = train_epochs, save_steps=200, checkpoint_dir=f'./{ckpt_path}/{model_name}-{split}-{lr}', alpha = lr, use_gpu = True, opt_method = opt_method)

trainer.run()
func.save_checkpoint(f'./{ckpt_path}/{model_name}-{split}-{lr}.ckpt')

func.load_checkpoint(f'./{ckpt_path}/{model_name}-{split}-{lr}.ckpt')
tester = Tester(model = func, data_loader = test_dataloader, use_gpu = True)
tester.run_link_prediction(type_constrain = False)
