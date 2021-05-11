import sys, os
sys.path = ['/home/enes/.conda/envs/gene/lib/python3.8/site-packages'] + sys.path
import pyslim
import numpy as np
import msprime as msp

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

if __name__ == '__main__':
    print('recap is starting\n')
    simID = int(sys.argv[1])
    N = int(sys.argv[2])
    Ne = int(sys.argv[3])
    r = float(sys.argv[4])
    mu = float(sys.argv[5])
    
    recap_path = 'trees/r'+str(simID)+'.trees'
    origi_path = 'trees/'+str(simID)+'.trees'

    ts = pyslim.load(origi_path).simplify()
    ts = recap(ts, Ne = Ne, r = r, SEED = simID)
    print('recap is done\n')
    ts = mutate(ts, mu = mu, SEED = simID)
    print('mutate is done\n')
    ts = subsample(ts, N, SEED = simID)
    print('subsample is done\n')
    ts.dump(recap_path)
    print('it\'s done\n')
    import subprocess
    print(subprocess.check_output('mv '+origi_path + ' ' + '/scratch/stats_dept_root/stats_dept1/enes/'+origi_path,
                                  stderr=subprocess.STDOUT,
                                  shell=True,
                                  encoding = 'utf-8'))

    