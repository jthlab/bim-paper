import sys
sys.path.append('../')
import numpy as np
import pandas as pd
import ruptures as rpt
import os

def segment_to_p(df, key, mu_g, cov, weights = 'gaussian'):
    from scipy.stats import norm
    r = len(cov)
    segment = 's'+key
        
    def cF(k):
        if k >= r:
            return 0
        else:
            return cov[k]

    def Zscan(x):
        nX = len(x)

        if weights == 'linear':
            w = np.arange(nX)
            w = np.min((nX-w, w+1), 0)
        elif weights == 'uniform':
            w = np.ones(nX)
        elif weights == 'gaussian':
            w = np.arange(1, nX+1)
            m = (nX+1)/2
            w = norm.pdf(5*(w-m)/(nX+1))  
        w = w/w.sum()

        x_bar =  (x*w).sum()

        w2 = w**2

        x_bar_var = w2.sum()*cF(0)
        for i in range(1, nX):
            cW = (w[:-i]*w[i:]).sum()   
            x_bar_var += 2*cF(i)*cW
        sem = np.sqrt(x_bar_var)

        return norm.cdf((x_bar-mu_g)/sem)
    
    df['p'+key] = df[key].copy()
    
    dx = df.groupby(segment, sort = False).agg({'Chromosome':'first', 
                                                'start':'first', 
                                                'end':'last', 
                                                'p'+key:Zscan,
                                                key:'mean'}).reset_index()
    
   
    key = 'p'+key
    
    cps = np.sort(np.r_[dx['start'].to_numpy(), dx['end'].to_numpy()])

    ps = dx[key].to_numpy()
    return cps, np.repeat(ps,2)

def segment_to_Z(df, key, mu_g, cov, weights = 'gaussian'):
    from scipy.stats import norm
    r = len(cov)
    segment = 's'+key
        
    def cF(k):
        if k >= r:
            return 0
        else:
            return cov[k]

    def Zscan(x):
        nX = len(x)

        if weights == 'linear':
            w = np.arange(nX)
            w = np.min((nX-w, w+1), 0)
        elif weights == 'uniform':
            w = np.ones(nX)
        elif weights == 'gaussian':
            w = np.arange(1, nX+1)
            m = (nX+1)/2
            w = norm.pdf(5*(w-m)/(nX+1))  
        w = w/w.sum()

        x_bar =  (x*w).sum()

        w2 = w**2

        x_bar_var = w2.sum()*cF(0)
        for i in range(1, nX):
            cW = (w[:-i]*w[i:]).sum()   
            x_bar_var += 2*cF(i)*cW
        sem = np.sqrt(x_bar_var)

        return (x_bar-mu_g)/sem
    
    df['Z'+key] = df[key].copy()
    
    dx = df.groupby(segment, sort = False).agg({'Chromosome':'first', 
                                                'start':'first', 
                                                'end':'last', 
                                                'Z'+key:Zscan,
                                                key:'mean'}).reset_index()
    
   
    key = 'Z'+key
    
    cps = np.sort(np.r_[dx['start'].to_numpy(), dx['end'].to_numpy()])

    zs = dx[key].to_numpy()
    return cps, np.repeat(zs,2)


def read_and_load(folder, popid, chrno, mode = ''):
    if mode != '':
        mode = 'C'+mode
    df2 = pd.read_csv(os.path.join(folder,str(popid)+'_'+str(chrno)+mode+'.csv'), comment = '#')
    df2.insert(0, "Chromosome", chrno)
    return df2

# segmenter
def cpd_genomescan(dc, key, popid, n_bkps = 8):
    
    df = dc.copy()
    df[key] = df[key].interpolate(method='nearest').fillna(method='ffill').fillna(method='bfill')
    y = df[key].to_numpy() # signal
    
    algo_c = rpt.KernelCPD(kernel="linear").fit(y/y.std()) 
    result = algo_c.predict(n_bkps = n_bkps)

    result = np.r_[[0], result]

    df['s'+key] = ['p'+str(popid)+'c'+str(chrno)+'.'+str(i) for i in np.repeat(np.arange(len(result)-1), np.diff(result))]
    return df

if __name__ == '__main__':
    
    popid = int(sys.argv[1])
    chrno = int(sys.argv[2])
   
    try:
        avglen = int(sys.argv[3])
    except:
        avglen = float(sys.argv[3])
    
    folder = sys.argv[4]
    
    try:
        mode = sys.argv[5]
    except:
        mode = ''
    
    locvars = ['Chromosome', 'start', 'end']
    dc = read_and_load(folder, popid, chrno, mode)
    stats = list(set(dc.keys()).difference(set(locvars+['N', 'path'])))

    if mode in ['medi', 'mean']:
        if mode == 'mean':
            dx = pd.read_csv('meanstats.csv')
        else:
            dx = pd.read_csv('medianstats.csv')

        statsM = [mode+'_'+i for i in stats]
        dc = dx[locvars+statsM].merge(dc, on = locvars, how = 'inner')

        for key in stats:
            dc[key+'C'+mode] = dc[key] - dc[mode+'_'+key]

        stats = [key+'C'+mode for key in stats]
    else:
        if mode != '':
            raise ValueError('unknown mode')
        else:
            pass

    n_kbps = np.round(dc['end'].iloc[-1]/avglen).astype('int')

    for key in stats:
        dc = dc[locvars].merge(cpd_genomescan(dc, key, popid = popid, n_bkps = n_kbps), on = locvars, how = 'inner')

    keep = locvars + stats + ['s'+key for key in stats]

    dc = dc[keep]
    
    save_path = folder+'s'+str(popid)+'_'+str(chrno)+mode+'.csv'
    dc.to_csv(save_path, index = False)
    print('success!')