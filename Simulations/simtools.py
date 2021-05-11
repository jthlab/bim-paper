import numpy as np

import sys

sys.path = ['/home/enes/.conda/envs/gene/lib/python3.8/site-packages'] + sys.path
sys.path.append('../')
import pandas as pd
import tsinfer
import msprime as msp
from functools import lru_cache

def recap(ts, Ne, r = 1e-8, SEED = 48):
    # recapitate slim outputs
    rts = ts.recapitate(recombination_rate = r, Ne=Ne, random_seed=SEED % 2**32 )
    return rts

def mutate(ts, mu = 1e-8, SEED = 1516):
    # adds mutation to not mutated trees
    mts = msp.mutate(ts, rate=mu, random_seed = SEED % 2**32 )
    return mts

def subsample(ts, N, SEED = 2342):
    # simplyfy trees to a subsample
    np.random.seed(SEED % 2**32)
    s1 = np.random.choice(ts.samples(), N, False)
    ts1 = ts.simplify(s1)
    return ts1

def infer_tree(ts):
    # takes a tree file and infer another tree by using genotype matrix
    GM = ts.genotype_matrix()
    muts = ts.mutations()
    positions = [next(muts).position for _ in range(GM.shape[0])] 
    L = GM.shape[0]
    seqlen = ts.get_sequence_length()
    with tsinfer.SampleData(sequence_length=seqlen) as sample_data:
        for i in range(L):
            sample_data.add_site(positions[i], GM[i])
    return tsinfer.infer(sample_data).simplify()

def get_betas(ts, btree):
    # calculates beta statistic and colless for true tree and inferred tree
        trees = ts.trees()
        ntrees = ts.num_trees
        betaest = []        
        for i in range(ntrees):
            Tree = next(trees) # traverse the tree sequence    
            start, end = Tree.interval
            res = {'start':start, 'end':end} # locaiton
            n, k = btree.tree_to_splits(Tree)['splits'] # tree -> splits
            mdl1 = btree.predict(n, k, n-2) # splits -> betahat weighted
            mdl2 = btree.predict(n, k) # splits -> betahat
            res['btree_w'] = mdl1['betahat']
            res['btree'] = mdl2['betahat']
            res['colless'] = mdl1['colless']            
            res.update() 
            betaest.append(res)
        df = pd.DataFrame(betaest)
        L = ts.get_sequence_length()
        df['w'] = (df['end'] - df['start'])/L
        btree_w = sum(df['w']*df['btree_w'])
        btree = sum(df['w']*df['btree'])
        colless = sum(df['w']*df['colless'])
        return btree_w, btree, colless

def calc_stats(ts, bsfs, btree):
    '''
    Statistics to be caclulated for the simulations. All simulations will share 
    the same statistics.

    Parameters
    ----------
    ts : tskit.trees.TreeSequence
    bsfs : BimSFS.bSFS
    btree : BimTree.bTree

    Returns
    -------
    ret : dict
        Dictionary of calcualted stats.

    '''
    
    if all([i is None for i in [ts, bsfs, btree]]):
        # test
        from BimTree import bTree
        from BimSFS import bSFS
        N = 200
        ts = msp.simulate(sample_size = N, Ne = 10000, length = 10000, 
                          recombination_rate = 1e-8, mutation_rate = 1e-8)
        btree = bTree(N, rho1 = 0., rho2 = 0.)
        bsfs = bSFS(N)
        
    
    # SFS based model
    sfs = ts.allele_frequency_spectrum(span_normalise=False, polarised=True)[1:-1]
    ret = bsfs.predict(sfs)
    ret.pop('h0')
    ret.pop('h1')
    ret['bsfs'] = ret['betahat']
    ret['nMuts'] = ts.num_mutations
    ret['nTrees'] = ts.num_trees
    del ret['betahat']
    
    # Tree based model with true geneology       
    beta, betawF, colless = get_betas(ts, btree)
    
    # Tree based model with inferred geneology     
    ibeta, ibetawF, icolless = get_betas(infer_tree(ts), btree)

    ret.update({'btreeF':betawF, 'ibetaF':ibetawF, 'btree':beta, 'ibtree':ibeta, 'Ic':colless, 'iIc': icolless})
    return ret

def ROC(ax, y0, y1, score_ascending = False, label = None):
    if not score_ascending:
        y0 = -y0
        y1 = -y1
    
    leny1 = len(y1)
    leny0 = len(y0)
    
    y_true = np.r_[np.zeros(leny0), np.ones(leny1)]
    y_score = np.r_[y0, y1]

    fpr, tpr, _ = roc_curve(y_true, y_score)
    
    ax.plot(fpr, tpr, label = label, linewidth = 3)

    return roc_auc_score(y_true, y_score)
    
    
def hist(y, color = 'blue', label = None):
    plt.hist(y, bins = 50, density = True, alpha = 0.5, color = color, label = label)