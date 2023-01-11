import torch
import numpy as np
import time
import import_ipynb
import torch.nn as nn
import pandas as pd 
from torch.distributions import Normal
from matplotlib import pyplot as plt
from sklearn.preprocessing import scale 
from sklearn.model_selection import train_test_split
import scipy.stats as scst
import pathlib,os
from scipy.special import expit
from sklearn import random_projection
import argparse
from joblib import Parallel, delayed
from torch.autograd import Variable
from torchvision import datasets, transforms
import torch.nn.functional as F
from scipy.linalg import svd
from scipy.stats import multinomial
from sklearn.metrics import roc_curve, auc
from sklearn.utils.validation import check_random_state
from scipy.special import gammainc

#### tools file content
#os.getcwd()
#os.chdir('cancer_result')

def sigmoid(z):
    return 1. / (1 + torch.exp(-z))


def log_gaussian(x, mu, sigma):
    """
        log pdf of one-dimensional gaussian
    """
    if not torch.is_tensor(sigma):
        sigma = torch.tensor(sigma)
    return float(-0.5 * np.log(2 * np.pi)) - torch.log(sigma) - (x - mu)**2 / (2 * sigma**2)


class MLPLayer(nn.Module):
    """
        Layer of our BNN
    """
    def __init__(self, input_dim, output_dim, rho_prior, device):
        # initialize layers
        super(MLPLayer, self).__init__()
        # set input and output dimensions
        self.input_dim = input_dim
        self.output_dim = output_dim
        self.w_mu = nn.Parameter(torch.Tensor(input_dim, output_dim).uniform_(-0.6,0.6))
        self.w_rho = nn.Parameter(torch.Tensor(input_dim, output_dim).uniform_(-0.6,0.6))
        self.b_mu = nn.Parameter(torch.Tensor(output_dim).uniform_(-6.0,-6.0))
        self.b_rho = nn.Parameter(torch.Tensor(output_dim).uniform_(-6.0,-6.0))
        self.rho_prior = rho_prior
        self.device = device
        self.w = None
        self.b = None
        self.kl = 0

    def forward(self, X):
        """
            For one Monte Carlo sample
            :param X: [batch_size, input_dim]
            :return: output for one MC sample, size = [batch_size, output_dim]
        """
        # sample weights and biases
        sigma_w = torch.log(1 + torch.exp(self.w_rho))
        sigma_b = torch.log(1 + torch.exp(self.b_rho))
        sigma_prior = torch.log(1 + torch.exp(self.rho_prior))

        epsilon_w = Normal(0, 1).sample(self.w_mu.shape)
        epsilon_b = Normal(0, 1).sample(self.b_mu.shape)
        epsilon_w = epsilon_w.to(self.device)
        epsilon_b = epsilon_b.to(self.device)

        self.w = self.w_mu + sigma_w * epsilon_w
        self.b = self.b_mu + sigma_b * epsilon_b
        output = torch.mm(X, self.w) + self.b.expand(X.size()[0], self.output_dim)

        kl_w = (torch.log(sigma_prior) - torch.log(sigma_w) +
                        0.5 * (sigma_w ** 2 + self.w_mu ** 2) / sigma_prior ** 2 - 0.5)

        kl_b = (torch.log(sigma_prior) - torch.log(sigma_b) +
                        0.5 * (sigma_b ** 2 + self.b_mu ** 2) / sigma_prior ** 2 - 0.5)

        self.kl = torch.sum(kl_w) + torch.sum(kl_b)

        return output
    
#### sparsefunc file content
#### sparsefunc file content
class SFunc(nn.Module):
    """
        Our BNN
    """
    def __init__(self, data_dim, hidden_dim1, target_dim, device, rho_prior=1):

        # initialize the network using the MLP layer
        super(SFunc, self).__init__()
        self.rho_prior = torch.Tensor([rho_prior]).to(device)
        self.device = device
        self.l1 = MLPLayer(data_dim, hidden_dim1, self.rho_prior, self.device)
        self.l1_sigmoid = nn.Sigmoid()
        self.l4 = MLPLayer(hidden_dim1, target_dim, self.rho_prior, self.device)
        self.target_dim = target_dim

    def forward(self, X):
        """
            output of the BNN for one Monte Carlo sample
            :param X: [batch_size, data_dim]
            :return: [batch_size, target_dim]
        """
        output = self.l1_sigmoid(self.l1(X))
        output = self.l4(output)
        return output

    def kl(self):
        # calculate the kl over all the layers of the BNN
        kl = self.l1.kl +  self.l4.kl
        return kl

    def sample_elbo(self, X, y, n_samples, num_batches):
        """
            calculate the loss function - negative elbo
            :param X: [batch_size, data_dim]
            :param y: [batch_size]
            :param n_samples: number of MC samples
            :return:
        """
        outputs = torch.zeros(n_samples, y.shape[0], self.target_dim).to(self.device)
        kls = 0.
        log_likes = 0.
        # make predictions and calculate prior, posterior, and likelihood for a given number of MC samples
        for i in range(n_samples):  # ith mc sample
            outputs[i] = self.forward(X)
            sample_kl = self.kl()  # get kl (a number)
            kls += sample_kl      
            # The following does not work when number of KL samples =1
            log_likes += -F.cross_entropy(outputs[i].squeeze(), y, reduction='sum')
        # calculate MC estimates of log prior, vb and likelihood
        kl_MC = kls/float(n_samples)
        # calculate negative loglikelihood
        nll_MC = - log_likes/float(n_samples)
        loss = kl_MC / num_batches + nll_MC
        return loss, outputs.squeeze()

        
    
def read_data(split_state):
    # The data has been generated from the following equation
    # y=I[e^{x_1}+x_2^2+5sin(x_3x_4)-3>0]
    # for n=3000 and p=20
    df=pd.read_excel("non_linear_classification_small.xlsx", index_col=None)
    df=np.array(df)
    X_orig = df[:, 1:-1]
    Y_orig = df[:, -1]
    X_til=scale(X_orig,with_mean=True,with_std=True)
    x_train, x_test, y_train, y_test = train_test_split(X_til, Y_orig, test_size=0.3, random_state=split_state) 
    x_train = torch.Tensor(x_train)
    y_train =  torch.flatten(torch.Tensor(y_train).long())
    x_test = torch.Tensor(x_test)
    y_test = torch.flatten(torch.Tensor(y_test).long())
    return x_train, x_test, y_train, y_test
    

def train_model(x_train,x_test,y_train,y_test,data_dim,hidden_dim1,target_dim,epochs,learning_rate,batch_size,jump):
    torch.set_default_dtype(torch.float64)
    num_batches = x_train.shape[0]/batch_size
    net = SFunc(data_dim=data_dim, hidden_dim1 = hidden_dim1, target_dim = target_dim, device=device).to(device)
    optimizer = torch.optim.Adam(net.parameters(), lr=learning_rate)
    for epoch in range(epochs): 
        permutation = torch.randperm(x_train.shape[0])
        for i in range(0, x_train.shape[0], batch_size):
            indices = permutation[i:i + batch_size]
            batch_x, batch_y = x_train[indices], y_train[indices]
            batch_x = batch_x.to(device)
            batch_y = batch_y.to(device)        
            loss, _ = net.sample_elbo(batch_x, batch_y, 1, num_batches)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step() 
        if (epoch+1)%jump==0:
            with torch.no_grad():
                loss1, pred = net.sample_elbo(x_train, y_train,10,1)
                pred=pred.mean(dim=0)
                x_test = x_test.to(device)
                y_test = y_test.to(device)
                _, pred2 = net.sample_elbo(x_test, y_test,10,1)
                pred2=pred2.mean(dim=0)
                pred = torch.max(pred, 1)[1].to(device)
                train_accur = np.mean(pred.detach().numpy() == y_train.detach().numpy())
                pred2 = torch.max(pred2, 1)[1].to(device)
                test_accur =  np.mean(pred2.detach().numpy() == y_test.detach().numpy())
                print('epoch: {}, loss: {}, train acc: {} test acc: {}'.format(epoch,loss1,train_accur,test_accur))            
    return loss1,train_accur,test_accur




parser = argparse.ArgumentParser()
parser.add_argument('--epochs', type=int, default = 2000, help='What is the the number of epochs')
parser.add_argument('--batch_size', type=int, default = 512, help='What is the batch_size')
parser.add_argument('--nodes', type=int, default = 32, help='what is the number of nodes in hidden layer')
parser.add_argument('--lr', type=float, default = 5e-3, help='what is the learning rate')
parser.add_argument('--jump', type=int, default = 100, help='what is the jump size')
parser.add_argument('--split', type=int, default =1, help='split state')
args = parser.parse_args()

if __name__ == '__main__':
    
    device = torch.device('cpu')
    torch.set_default_dtype(torch.float64)
    x_train, x_test, y_train, y_test=read_data(args.split)
    data_dim=x_train.shape[1]
    hidden_dim1=args.nodes
    learning_rate=args.lr
    jump=args.jump
    epochs=args.epochs
    target_dim=2
    batch_size=args.batch_size
    loss1, train_accur, test_accur=train_model(x_train, x_test, y_train, y_test, data_dim, hidden_dim1, target_dim, epochs, learning_rate, batch_size, jump)
    
