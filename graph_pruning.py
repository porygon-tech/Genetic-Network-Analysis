#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 11:17:25 2023
This script is intended for educational purposes at the block course ECO 365 Ecological Networks from the University of Zurich.
@author: Miguel Roman
"""

#python3 add_noise.py --f=lopho_graph_adj.csv --out=noisy_adj.csv
#python3 graph_pruning.py --method=mincov --f=noisy_adj.csv --out=pruned_mincov.csv

import pandas as pd
import networkx as nx
import numpy as np
import argparse
#from os import chdir
from pathlib import Path
#chdir('/home/roman/LAB/Teaching') #this line is for Spyder IDE only
root = Path(".")

#%% ARGS
parser = argparse.ArgumentParser(description='''
Script developed by Miguel Roman.
''')

parser.add_argument('--method', help='algorithm for pruning. Can be \'threshold\' or \'singlecomp\' or \'mincov\'',
                    required=True)
parser.add_argument('--out', help='Name of the output file',
                    default='pruned.csv')
parser.add_argument('--f', help='Name of the input file',
                    required=True)
parser.add_argument('--param', help='params for algorithm.',
                    default=0)

#args = parser.parse_args(['-h'])
args = parser.parse_args()
METHOD = args.method
FILE_OUT = args.out
FILE_IN = args.f
PARAM = float(args.param)

#%% SETUP TEST
'''
nodeData =  pd.read_csv(root / 'lopho_nodeData.csv')
nodeData = nodeData.set_index('name')
nodeData['pos'] = list(nodeData[['Longitude','Latitude']].to_numpy())

#df = pd.read_csv(root / 'lopho_graph_adj.csv')
df = pd.read_csv(root / 'noisy_adj.csv')
df = pd.read_csv(root / 'pruned.csv')

df=df.rename(columns={'Unnamed: 0' : 'name'}) 
df=df.set_index('name')
'''
#%% SETUP
nodeData =  pd.read_csv(root / 'lopho_nodeData.csv')
nodeData = nodeData.set_index('name')
nodeData['pos'] = list(nodeData[['Longitude','Latitude']].to_numpy())

#df = pd.read_csv(root / 'lopho_graph_adj.csv')
df = pd.read_csv(root / FILE_IN)
df=df.rename(columns={'Unnamed: 0' : 'name'}) 
df=df.set_index('name')

#%% LOAD / RELOAD
G = nx.from_pandas_adjacency(df)
nx.set_node_attributes(G, nodeData.to_dict('index'))

print('INITIAL N_edges ' + str(len(G.edges)))
#%% TEST VIEW
'''
sizes = list( 5*np.array( list(nx.get_node_attributes(G, 'size').values())))
linewidths = list( 5/np.array( list(nx.get_edge_attributes(G,'weight').values())))

nx.draw(G, 
        nx.get_node_attributes(G, 'pos'), 
        with_labels=True, 
        node_size=sizes, 
        width=linewidths)
'''
#%% PRUNE 1 
#remove = [node for node,degree in dict(G.degree()).items() if degree > 8]
if METHOD == 'threshold':
    print("METHOD: threshold")
    threshold = PARAM
    remove = [edge for edge,weight in nx.get_edge_attributes(G,'weight').items() if weight > threshold]
    G.remove_edges_from(remove)
    adjacency = nx.adjacency_matrix(G).todense()

#%% PRUNE 2

if METHOD == 'singlecomp':
    print("METHOD: singlecomp")
    #threshold = PARAM
    while nx.number_connected_components(G) == 1:
        gsave = G.copy()
        d = nx.get_edge_attributes(G,'weight')
        remove = max(d, key=d.get)
        G.remove_edges_from([remove])
    
    G = gsave
    adjacency = nx.adjacency_matrix(G).todense()


#%% PRUNE 3
#D= pd.read_csv(root / 'gfake_adj.csv').rename(columns={'Unnamed: 0' : 'name'}).set_index('name').to_numpy()
from scipy.stats import chi2
if METHOD == 'mincov':
    print("METHOD: minimum covariance")
    print(PARAM)
    D = nx.adjacency_matrix(G).todense()
    D=1/D; np.fill_diagonal(D,0) # invert, to obtain the closeness instead of genetic distance matrix
    n = D.shape[0]
    avgD1 = D.mean(1)[:,np.newaxis] * np.ones((1,n)) #the axis for the mean calculation does not matter since the matrix is symmetric
    avgD2 = D.mean(0)[np.newaxis,:] * np.ones((n,1)) 
    avgDg=D.mean()                 * np.ones(D.shape)# global avg
    C = -0.5*(D - avgD1 - avgD2 + avgDg) # this is the covariance matrix as specified in Gower 1966 and https://doi.org/10.1111/j.1365-294X.2004.02177.x
    # P = np.linalg.pinv(C)
    # P=P/P.sum()
    
    P = np.linalg.pinv(C.T @ C) * C.T
    denom = np.sqrt(np.outer(np.diag(P),np.diag(P)))
    R = - P / denom
    np.fill_diagonal(R,1)
    eed = -n*np.log(1-R**2)#edge exclusion deviance
    
    #-------------
    if PARAM > 0.:
        print("Performing EED removal based on c.i.=" + str(100*(1-PARAM)) +"%...")
        #showdata(D*(eed>threshold))
        #showdata(eed)
        threshold = chi2.ppf(1-PARAM, df=1)
        adjacency = D*(eed>threshold)
        
    #-------------
    else:
        print("No parameter specified. Performing EED sequential removal...")
        dft = pd.DataFrame(eed)
        dft.index = dft.columns = nodeData.index
        Geed = nx.from_pandas_adjacency(dft)
        nx.set_node_attributes(Geed, nodeData.to_dict('index'))
        while nx.number_connected_components(Geed) == 1:
            #print(nx.number_connected_components(Geed))
            gsave = Geed.copy()
            d = nx.get_edge_attributes(Geed,'weight')
            remove = min(d, key=d.get)
            Geed.remove_edges_from([remove])
            #print('edges left: ' + str(len(Geed.edges)))
        
        Geed = gsave
        
        eedTopo = nx.adjacency_matrix(Geed).todense()
        adjacency = D*(eedTopo>0)
        adjacency[np.where(adjacency != 0)] **=-1

'''
showdata(eed)
showdata(adjacency)
nx.is_connected(nx.from_numpy_array(adjacency))
nx.draw(nx.from_numpy_array(adjacency))
'''
#%% PRUNE 4

if METHOD == 'min_spanning_tree':
    mint = nx.minimum_spanning_tree(G)
    adjacency = nx.adjacency_matrix(mint).todense()


                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            

'''
d = nx.get_edge_attributes(Geed,'weight')
remove = min(d, key=d.get)
Geed.remove_edges_from([remove])
showdata(nx.adjacency_matrix(Geed).todense())
len(Geed.edges)
nx.draw(Geed)
nx.number_connected_components(Geed)
showdata(eedTopo)
showdata(C)
showdata(R)
showdata(R_zeroDiag)
showdata(R_zeroDiag<1)
showdata(adjacency)
nx.number_connected_components(adjacency)

'''

#%% 
'''
A = pd.read_csv(root / 'lopho_graph_adj.csv').rename(columns={'Unnamed: 0' : 'name'}).set_index('name').to_numpy()
A = pd.read_csv(root / 'pruned.csv').rename(columns={'Unnamed: 0' : 'name'}).set_index('name').to_numpy()

maxv=np.max(A)
minv=np.min(A[np.nonzero(A)])
meanv=np.mean(A[np.nonzero(A)])
noise = np.random.uniform(low=maxv*1.5, high=2*maxv,size=A.shape)
#noise = np.random.uniform(low=minv, high=maxv,size=A.shape)
#noise = np.random.uniform(low=meanv-np.mean((maxv-meanv,meanv-minv)), high=meanv+np.mean((maxv-meanv,meanv-minv)),size=A.shape)

noise[np.where(A != 0)] = 0
noise = np.triu(noise,1)+np.triu(noise,1).T
showdata(noise)
D = A+noise


#D = nx.adjacency_matrix(G).todense()

showdata(A)
showdata(D)


n = D.shape[0]
avgD1 = D.mean(1)[:,np.newaxis] * np.ones((1,n)) #the axis for the mean calc. does not matter since the matrix is symmetric
avgD2 = D.mean(0)[np.newaxis,:] * np.ones((n,1)) #the axis for the mean calc. does not matter since the matrix is symmetric
avgDg=D.mean()                 * np.ones(D.shape)# global avg
C = -0.5*(D - avgD1 - avgD2 + avgDg) # this is the covariance matrix as specified in Gower 1966 and https://doi.org/10.1111/j.1365-294X.2004.02177.x
np.linalg.cond(C)


if np.linalg.cond(C) < 10e2: # if the condition number is very large (significantly higher than one), the matrix is ill-conditioned and it may be difficult to compute its inverse accurately.
    P = np.linalg.inv(C)
    diag_sqrt = np.sqrt(np.diag(P))
    diag_mat = np.diag(diag_sqrt)
    prec_sqrt_mat = np.dot(diag_mat, P)
    corr_mat = C / np.dot(diag_mat, P)

    np.fill_diagonal(R, 1)
else:
    stds = np.sqrt(np.diag(C))
    R = C / np.outer(stds, stds)
    #R = np.corrcoef(C, rowvar=False)

R_zeroDiag = R.copy() 
np.fill_diagonal(R_zeroDiag,0)
eed = -n*np.log(1-R_zeroDiag**2)#edge exclusion deviance


#%%
showdata(eed)


dft = pd.DataFrame(eed)
dft.index = dft.columns = nodeData.index
G = nx.from_pandas_adjacency(dft)
nx.set_node_attributes(G, nodeData.to_dict('index'))
while nx.number_connected_components(G) == 1:
    gsave = G.copy()
    d = nx.get_edge_attributes(G,'weight')
    remove = min(d, key=d.get)
    G.remove_edges_from([remove])

G = gsave

eedTopo = nx.adjacency_matrix(G).todense()
adjacency = D*(eedTopo>0)
showdata(adjacency)

#i need to get a correlation matrix from a covariance matrix through a precision matrix  in python
#how is the precision matrix standardized to a correlation matrix, R, using normal matrix routines?
#but what if the precision matrix has negative numbers? the sqrt will not work in that case, and would be very common 

lopho = pd.read_csv(root / 'lopho_graph_adj.csv').rename(columns={'Unnamed: 0' : 'name'}).set_index('name')
showdata(lopho.to_numpy())
showdata(D)
showdata(A)
showdata(D*(eed>7))
showdata(D*((eed/D)>.2))
showdata(C)

showdata(R)

df = pd.DataFrame(adjacency)
df.index = df.columns = nodeData.index
G = nx.from_pandas_adjacency(df)
nx.set_node_attributes(G, nodeData.to_dict('index'))

sizes = list( 5*np.array( list(nx.get_node_attributes(G, 'size').values())))
linewidths = list( 5/np.array( list(nx.get_edge_attributes(G,'weight').values())))

nx.draw(G, 
        nx.get_node_attributes(G, 'pos'), 
        with_labels=True, 
        node_size=sizes, 
        width=linewidths)
'''
#%%
'''
while nx.sigma(G) <= 1:
    gsave = G.copy()
    d = nx.get_edge_attributes(G,'weight')
    remove = max(d, key=d.get)
    G.remove_edges_from([remove])
    
G = gsave
'''
#%% SAVE ADJ
df = pd.DataFrame(adjacency)
df.index = df.columns = nodeData.index
df.to_csv(root / FILE_OUT)

print("Saved to " + FILE_OUT)
