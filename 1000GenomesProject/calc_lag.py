import sys, os
sys.path.append('/home/enes/pyslurm/')
import numpy as np
import pandas as pd

if __name__ == '__main__':
    
    key = sys.argv[1] 
    no_lag = int(sys.argv[2])
    folder = sys.argv[3]
    try:
        mode = sys.argv[4]
    except:
        mode = ''
    
        
    def lag(x, r):
        e = len(x)-r
        covmat = np.cov(x[r:], x[:e])
        if len(covmat.shape) == 2:
            return covmat[0, 1]
        else:
            return covmat
    
    poplag = {}
    for pop in range(26):
        path = os.path.join(folder, 's'+str(pop)+mode+'.csv')
        dc = pd.read_csv(path)
        dc = dc[dc[key].notna()]
        locvars = ['Chromosome', 'start', 'end']
        mu_g = dc[key].mean()
        vals = [i[1].to_numpy() for i in dc.groupby('Chromosome')[key]]
        no_arr = len(vals)
        lagmat = np.zeros((no_arr, no_lag))
        for i in range(no_arr):
            for j in range(no_lag):
                lagmat[i, j] = lag(vals[i], j)

        w = dc.groupby('Chromosome').size().to_numpy()
        w = w/w.sum()
        poplag[pop] = (lagmat*w[:, None]).sum(0)
    pd.DataFrame(poplag).to_csv('lags/'+key+'.csv', index = False)