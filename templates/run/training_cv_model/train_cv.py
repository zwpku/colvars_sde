#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
import torch
import math 
import torch.nn as nn
import random
from tqdm import tqdm
import os
import sys
import time
import yaml
import argparse

# this line may need to be changed 
sys.path.append('../../../colvars-finder')

from colvarsfinder.core import AutoEncoderTask, EigenFunctionTask, RegAutoEncoderTask
from colvarsfinder.nn import AutoEncoder, EigenFunctions, RegAutoEncoder, RegModel 
from colvarsfinder.utils import WeightedTrajectory, integrate_sde_overdamped, calc_weights

def set_all_seeds(seed):
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
        np.random.seed(seed)
    random.seed(seed)

def dict2namespace(config):
    namespace = argparse.Namespace()
    for key, value in config.items():
        if isinstance(value, dict):
            new_value = dict2namespace(value)
        else:
            new_value = value
        setattr(namespace, key, new_value)
    return namespace

def show_trajectory(trajectory, x_domain, y_domain, save_fig=True):
### Plot the trajectory 
    fig = plt.figure(figsize=(14,6))
    ax0 = fig.add_subplot(1, 2, 1)
    ax1 = fig.add_subplot(1, 2, 2)
    nx = ny = 200

    dx = (x_domain[1] - x_domain[0]) / nx
    dy = (y_domain[1] - y_domain[0]) / ny

    h, xedges, yedges = np.histogram2d(trajectory[:,0], trajectory[:,1], bins=[nx, ny], range=[[x_domain[0],x_domain[1]],[y_domain[0],y_domain[1]]], density=False)
    X, Y = np.meshgrid(xedges, yedges)
    im = ax0.pcolormesh(X, Y, h.T,  cmap='RdBu_r', shading='auto',norm=colors.LogNorm(1,1000) )
    cbar = fig.colorbar(im, ax=ax0, shrink=1.0)
    cbar.ax.tick_params(labelsize=15)

    ax0.set_title('histogram of data', fontsize=20)
    ax0.set_xlim(x_domain)
    ax0.set_ylim(y_domain)
    ax0.set_xlabel(r'$x_1$',fontsize=20)
    ax0.set_ylabel(r'$x_2$',fontsize=20, rotation=0)
    ax0.tick_params(axis='both', labelsize=20)

    ax1.scatter(trajectory[:,0], trajectory[:,1])
    ax1.set_title('scatter plot', fontsize=20)
    ax1.set_xlim(x_domain)
    ax1.set_ylim(y_domain)
    ax1.set_xlabel(r'$x_1$',fontsize=20)
    ax1.set_ylabel(r'$x_2$',fontsize=20, rotation=0)
    ax1.tick_params(axis='both', labelsize=20)

    if save_fig :
        filename = f"./traj.jpg" 
        fig.savefig(filename)
        print (" trajectory plot saved to file: %s" % filename)

def plot_cv(cv_model, x_domain, y_domain, model_path):
    gridx = np.linspace(x_domain[0], x_domain[1], 100)
    gridy = np.linspace(y_domain[0], y_domain[1], 100)
    x_plot = np.outer(gridx, np.ones(100)) 
    y_plot = np.outer(gridy, np.ones(100)).T 
    # prepare data
    x2d = torch.from_numpy(np.concatenate((x_plot.reshape(100 * 100, 1), y_plot.reshape(100 * 100, 1)), axis=1)).float()
    # evaluate model on data
    cv_on_grid = cv_model(x2d).detach().numpy()

    cv = cv_on_grid[:,0].reshape(100,100)
    #print ( "min and max values of %dth dimension of encoder: (%.4f, %.4f)" % (idx, encoder.min(), encoder.max()) )

    fig = plt.figure(figsize=(12,5))
    ax0 = fig.add_subplot(1, 2, 1, projection='3d')
    ax1 = fig.add_subplot(1, 2, 2)  
    ax0.plot_surface(x_plot, y_plot, cv, cmap='coolwarm', edgecolor='none')

    ax0.set_xlabel(r'$x_1$',fontsize=20)
    ax0.set_ylabel(r'$x_2$',fontsize=20, rotation=0)

    im = ax1.pcolormesh(x_plot, y_plot, cv, cmap='coolwarm',shading='auto')
    contours = ax1.contour(x_plot, y_plot, cv, 15, colors='black')
    ax1.clabel(contours, inline=True, fontsize=8)
    fig.colorbar(im, ax=ax1)

    ax1.set_xlabel(r'$x_1$',fontsize=20)
    ax1.set_ylabel(r'$x_2$',fontsize=20, rotation=0)

    fig_name = f"{model_path}/cv.jpg"
    fig.savefig(fig_name)
    #print ( "encoder profiles saved to file: %s" % fig_name )
    plt.close()

def read_parameter_config(config_file):

    with open(config_file, 'r') as f:
        args = yaml.safe_load(f)
    args = dict2namespace(args)

    return args

def load_traj(args):

    traj_fullname = os.path.join(args.data_path, args.traj_filename)

    traj_weight_fullname = None
    if args.weight_filename != None :
        traj_weight_fullname = os.path.join(args.data_path, args.weight_filename)
        if os.path.exists(traj_weight_fullname) == False : 
            traj_weight_fullname = None

# construct trajectory class
    traj = WeightedTrajectory(traj_filename=traj_fullname, weight_filename=traj_weight_fullname)

    return traj

if __name__ == '__main__':

    if len(sys.argv) < 2:
        print ("Error: config file not provided.")
        print ("Usage: ./train_cv.py config_file")
        exit(1)

    config_file = sys.argv[1]
        
    if os.path.exists(config_file) == False : 
        print (f"Error: config file ({config_file}) not found!")
        exit(1)

    args = read_parameter_config(config_file)

    print ("\n------- Parameters -------")

    for key, value in vars(args).items():
        print (f"{key} = {value}")

    set_all_seeds(args.seed)

    traj = load_traj(args)

    dim = traj.trajectory.shape[1]

    print (f"\nTrajectory:\n dim={dim}, len={traj.n_frames}")

    show_trajectory(traj.trajectory, args.x_domain, args.y_domain, save_fig=args.save_traj_fig)

    # use raw position data by setting preprocessing layer to identity
    pp_layer = torch.nn.Identity()

    e_dims = [dim] + args.e_dims + [args.k]
    d_dims = [args.k] + args.d_dims + [dim]

    model = AutoEncoder(e_dims, d_dims)

    print ("\nModel:", model)
    model_path = os.path.join(f'autoencoder-k={args.k}-' + time.strftime("%Y-%m-%d-%H:%M:%S", time.localtime()))

# define training task
    train_obj = AutoEncoderTask(traj, pp_layer, model, model_path, learning_rate=args.learning_rate, 
                                batch_size=args.batch_size, test_ratio=args.test_ratio, num_epochs=args.num_epochs,  verbose=False)

# train autoencoder
    train_obj.train()

# get cv model
    cv = train_obj.colvar_model()

# display results
    if dim == 2 :
        plot_cv(cv, args.x_domain, args.y_domain, model_path)

