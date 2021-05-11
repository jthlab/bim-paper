import sys, os
import numpy as np
import pandas as pd
from scipy.stats import norm

def plot_pscan(df, cov, key, thres = 1e-5, save_loc = None):
    mu_g = df[key].mean()
    r = len(cov)
    segment = 's'+key
    
    
    def cF(k):
        if k >= r:
            return 0
        else:
            return cov[k]

    def Zscan(x, weights = 'gaussian'):
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
    
    #plt.figure(figsize = (15,5))
    
    key = 'p'+key
    
    for tail in ['upper', 'lower']:
        import matplotlib.pyplot as plt
        plt.style.use('ggplot')
        fig, ax = plt.subplots(22, 1, figsize = (15,70), sharex=True)
        for i in range(1, 23):
            da = dx[dx['Chromosome'] == i].copy()
            
            if tail == 'upper':
                da[key] = 1-da[key].to_numpy()

            write = da[da[key]<thres]
            xpos = 0.5*(write['start']+write['end'])
            ypos = write[key]
            word = write[segment]

            cps = np.sort(np.r_[da['start'].to_numpy(), da['end'].to_numpy()])

            ps = da[key].to_numpy()
            ax[i-1].plot(cps, np.repeat(ps,2), linewidth = 0.5)
            ax[i-1].set_ylabel('Chromosome '+str(i))
            for x, y, w in zip(xpos, ypos, word):
                ax[i-1].text(x, y, str(w), rotation = 90)

            #ax[i,j].set_ylim(-1,1)
            ax[i-1].set_yscale('log')
        ax[i-1].set_xlabel('Genomic Position')
        if save_loc is None:
            plt.show()
        else:
            plt.tight_layout()
            fo, fi = os.path.split(save_loc)
            plt.savefig(os.path.join(fo, tail+'_'+fi), dpi = 100)
            plt.close()
    return dx

def genes_in_segment(df, segment):
  
    dfo = df.copy()
    df = df.copy()
    
    ann = pd.read_csv('gene_result.txt', sep="\t")
    ann = ann.rename(columns = {'start_position_on_the_genomic_accession': 'start',
                                'end_position_on_the_genomic_accession': 'end'})[['chromosome', 'start', 'end', 'Symbol']]
    ann = ann[~((ann['chromosome'] == 'X') | ((ann['chromosome'] == 'Y')))]

    ann['chromosome'] = ann['chromosome'].astype('int')
    ann = ann.sort_values(['chromosome', 'start', 'end']).reset_index(drop = True)
    
    
    df['start'] = df['Chromosome']+(df['start']/1e9)
    df['end'] = df['Chromosome']+(df['end']/1e9)

    dfg = df.groupby(segment, sort=False)
    
    ss = dfg['start'].min().to_numpy()
    se = dfg['end'].max().to_numpy()
    segs = iter(list(dfg.groups.keys()))
    
    gss = (ann['chromosome']+(ann['start']/1e9)).to_numpy()
    ges = (ann['chromosome']+(ann['end']/1e9)).to_numpy()
    gnames = ann['Symbol'].to_list()

    gss = np.r_[gss, 1e12]
    ges = np.r_[ges, 1e12]
    gnames = gnames + ['']
    segnames = []
    ni = len(ss)
    nk = len(ges)
    
    k = 0
    name_list = []

    for i in range(ni):
        segnames.append(next(segs))
        rs = ss[i] + 0.5/1e9
        re = se[i] - 0.5/1e9

        gs = gss[k]
        ge = ges[k]

        genes = []

        while(rs>ge):
            k = k+1
            gs = gss[k]
            ge = ges[k]

        while(re>gs):
            genes.append(gnames[k])
            k = k + 1
            gs = gss[k]
            ge = ges[k]

        name_list.append(', '.join(genes))   
    
    dfo['genes'] = name_list
    dfo = dfo.set_index(segment)
    return dfo

if __name__ == '__main__':
    
    pop = int(sys.argv[1])
    key = sys.argv[2]
    
    try:
        mode = sys.argv[3]
    except:
        mode = ''
    
    cov = pd.read_csv('lags/'+key+'.csv').loc[:,str(pop)].to_numpy()
    dz = pd.read_csv('merged/s'+str(pop)+mode+'.csv')
    dx = plot_pscan(dz, cov, key, thres = 0.05/dz['s'+key].unique().size, save_loc = 'plots/'+str(pop)+'_'+key+'.jpg')
    genes_in_segment(dx, 's'+key).to_csv('gene_scores/'+str(pop)+'_'+key+'.csv')                 