import re, subprocess, os, sys, json
import tskit
import numpy as np
import pandas as pd

from random import randint
from sklearn.metrics import roc_curve, roc_auc_score

sys.path.append('/home/enes/bim/')
from utils import InferEta

def get_simIDs(cID, setID, nrep):
    if any([cID>9, setID>99, nrep>10000]):
        raise ValueError('Wrong IDs')
    start = cID*1_000_000 + setID*10_000  
    ids = np.arange(start, start+nrep)
    
    setID = str(setID); cID = str(cID)
    if len(setID) == 1:
        setID = '0'+setID
    regex = 'r'+cID+setID+'[0-9]{4}.trees'
    return ids, regex

def treepaths(regex, folder = 'trees'):
    r = re.compile(regex)
    files = os.listdir(folder)
    tpaths = []
    for file in files:
        m = r.match(file)
        if m is not None:
            tpaths.append(os.path.join(folder,m.group()))
    return tpaths

def get_eps(Ne, gr, t, constant_gen = 0):
    dipNe = 2*Ne
    rep = gr+1

    epshist = []
    for i in range(t): 
        
        if i>constant_gen:
            dipNe1 = dipNe*rep
        else:
            dipNe1 = dipNe                        
        
        epshist.append(dipNe)
        dipNe = dipNe1
        
    return epshist

def ROC(ax, y0, y1, score_ascending = False, label = None):
    if not score_ascending:
        y0 = -y0
        y1 = -y1
    
    leny1 = len(y1)
    leny0 = len(y0)
    
    y_true = np.r_[np.zeros(leny0), np.ones(leny1)]
    y_score = np.r_[y0, y1]

    fpr, tpr, _ = roc_curve(y_true, y_score)
    
    auc = roc_auc_score(y_true, y_score)
    ax.plot(fpr, tpr, label = label +' ('+str(round(auc, 3))+')', linewidth = 3)
    
    return auc
    
    
def hist(y, color = 'blue', label = None):
    plt.hist(y, bins = 50, density = True, alpha = 0.5, color = color, label = label)

class PSlim:

    def __init__(self, Args):
        '''
        Parameters
        ----------
        Args : dict
            This is dictionary for the variables that will used for the simulator.
            Check slim and msprime methods to know which vars they need
        '''
        
        self.Args = Args
        self.gen_slimFile()

    def gen_slimFile(self):
        '''
        This is the slim file gnerator method. Args should have a key called slimTxt
        this is the text file for slim file. The variables coded as [<var>] in this
        text will be replaced by Args[<var>] then the resulting file will be saved 
        to Args['slimTxt'].replace('.txt', '.slim'). 
        '''
        Args = self.Args
        file = open(Args['slimTxt'])
        slimFile = file.read()
        file.close()
        L = re.split(r'(\[[^\]]+\])', slimFile)
        
        name = []
        for i in range(1, len(L), 2):
            varname = L[i][1:-1]
            L[i] = str(Args[varname])
            
            if len(L[i]) < 20:
                name.append(L[i])
        
        Args['slimFile'] = 'slimfiles/'+'_'.join(name)+Args['slimTxt'][:-4]+'.slim'
        self.Args = Args
        
        file1 = open(Args['slimFile'], 'w')
        file1.writelines(L)
        file1.close()
    
    def sim(self, simID, itrees = True):
        '''
        This is the method to run the slim simulation. This will run the slimFile
        

        Parameters
        ----------
        simID : int
            unique simulaiton id. will be used as a SEED
        savelog : bool, optional
            save the log to simInfo. The default is True.

        Returns
        -------
        int
            jobid from the slurm if srun provided, if srun is None then system out.

        '''
        Args = self.Args.copy()
        
        slimFile = Args['slimFile']
        simID = str(simID)
        
        job1 = 'slim -s '+simID+' -d simID='+simID+' '+slimFile
        job2 = 'python recap.py '+' '.join([str(i) for i in [simID, Args['N'], Args['Ne'], Args['r'], Args['mu']]])
        jobs = [job1, job2]
        
        if itrees:
            rtrees = 'trees/r'+simID+'.trees'
            jobs.append('python infertrees.py '+rtrees)
        
        srun = Args['srun']
        if srun is None:
            out = subprocess.check_output(' && '.join(jobs),
                                          stderr=subprocess.STDOUT,
                                          shell=True,
                                          encoding = 'utf-8')
        else:
            out = srun.run(jobs)
        
        return out
    
class Experiment:
    def __init__(self, cID, nrep, Args):
        self.cID = cID
        self.nrep = nrep
        self.Args = Args
        self.setids = setids = list(Args.keys())
        self.neutrals = [setid for setid in setids if all([(Args[setid]['s'] == 0), (Args[setid]['h'] == 0.5)])]
        
        self.outs = {}
        self.simIDs = {}
        self.setReg = {}
        self.df = {}
        for i, setid in enumerate(setids):
            self.simIDs[setid], self.setReg[setid] = get_simIDs(cID, i, nrep)
            self.outs[setid] = []
            self.df[setid] = 'outs/'+str(cID)+str(setid)+'.csv'  
            
        self.AFS = {setid:np.zeros(Args[setid]['N']-1, dtype = 'int') for setid in setids}
        
    def sim(self):
        setids = self.setids
        Args = self.Args
        simIDs = self.simIDs
        
        simJobs = {}
        for setid in setids:
            pslim = PSlim(Args[setid])
            simJobs[setid] = []
            for simID in simIDs[setid]:
                simJobs[setid].append(pslim.sim(simID))
        
        self.simJobs = simJobs
        print('If you are using HPC (srun is not None) check the jobs!')
        
    def calc_sfs(self):
        setReg = self.setReg
        setids = self.setids
        Args = self.Args
        AFS = {setid:np.zeros(Args[setid]['N']-1, dtype = 'int') for setid in setids}
        
        for setid in setids:
            regex = setReg[setid]
            tpaths = treepaths(regex)

            for tpath in tpaths:
                ts = tskit.load(tpath)
                AFS[setid] += ts.allele_frequency_spectrum(span_normalise=False, polarised=True)[1:-1].astype('int')
        
        self.AFS = AFS
        
    def train_eta(self):
        Args = self.Args
        AFS = self.AFS
        neutrals = self.neutrals
        
        t = np.logspace(np.log10(1), np.log10(20000), 100)
        t = np.concatenate((np.array([0]), t)) # breakpoints for the time (by generations)

        a1 = 0. #sequential l1 penalty 
        a2 = 1. #sequential l2 penalty
        
        self.ebl = {}
        for setid in neutrals:
            sfs = AFS[setid]

            inferEta = InferEta(Args[setid]['N'], t, a1 = a1, a2 = a2)
            etaout = inferEta.predict(sfs, maxiter = 1000)
            eps = etaout.x
            ebl = inferEta.get_esfs(eps)
            self.ebl[setid] = ebl
            
            ti = t
            ai = np.float16(etaout.x)
            diff = np.diff(ai).nonzero()[0]
            ti = np.r_[ti[0], ti[diff + 1]]
            ai = np.r_[ai[0], ai[diff + 1]]

            outname = Args[setid]['etapath']
            with open(outname, 'w') as fp:
                json.dump({0:{'t':ti.tolist(), 'a':(1/ai).tolist()}}, fp)
        print('Done!')       
        
    def est(self, BIM, setid, now = 2, srun = None, arg = ""):
        
        Arg = self.Args[setid]
        eta = Arg['etapath']
        simIDs = self.simIDs[setid]
        simIDs = ['trees/r'+str(i)+'.trees' for i in simIDs] + ['trees/ir'+str(i)+'.trees' for i in simIDs]
        n = len(simIDs)
        simIDs = iter(simIDs)

        jc = now*[n//now]
        remainings = n-sum(jc)
        j=0
        for i in range(remainings):
            jc[j] += 1
            j += 1

        msgs = []
        self.outs[setid] = []
        for j in jc:
            out = 'outs/'+str(randint(1e8, 1e9-1))+'.csv'
            self.outs[setid].append(out)
            trees = ','.join([next(simIDs) for _ in range(j)])
            job = ' '.join(['python', BIM, trees, str(Arg['N']), '--stats=all', '--eta='+eta, '--out='+out])
            job = job+' '+arg
            
            if srun is None:
                msg = subprocess.check_output(job,
                                              stderr=subprocess.STDOUT,
                                              shell=True,
                                              encoding = 'utf-8')
            else:
                msg = srun.run(job)

            msgs.append(msg)  
        
        print('If you are using HPC (srun is not None) check the jobs!')
        
        return msgs
    
    def merge_outs(self, setid):
        outs = self.outs[setid]
        
        name = self.df[setid]
        df = pd.concat([pd.read_csv(out, comment = '#') for out in outs]).reset_index(drop = True)      
        
        dfr = df[[i[0] == 'r' for i in df['path']]].copy()
        dfi = df[[i[0] == 'i' for i in df['path']]].copy()

        dfr['path'] = [int(re.findall('[0-9]+', i)[0]) for i in dfr['path']]
        dfi['path'] = [int(re.findall('[0-9]+', i)[0]) for i in dfi['path']]
        dfi = dfi.rename(columns = {'Colless':'iColless', 'btree':'ibtree'})[['path', 'iColless', 'ibtree']]
        df = dfr.merge(dfi, on = 'path')
        df.sort_values('path').to_csv(name, index = False) 
               
        [os.remove(out) for out in outs] 