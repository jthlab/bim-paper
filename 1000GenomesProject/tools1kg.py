def Allele_Freq(ts):
    return {i:[ts.allele_frequency_spectrum(sample_sets=[ts.samples(i)], 
                                            polarised=True, 
                                            span_normalise=False)[1:-1].astype('int')] for i in range(26)}
def path_1000G(CHR):
    return '/nfs/turbo/lsa-enes/singletonsadded/1kg_chr'+str(CHR)+'_singletons.trees'

def GWAFS():
    # genome-wide AFS
    import tskit
    import numpy as np
    from tqdm.notebook import trange

    SFS = {i:[] for i in range(26)}

    for CHR in trange(1, 23):
        ts = tskit.load(path_1000G(CHR))
        sfs_chr = Allele_Freq(ts)
        for i in range(26):
            SFS[i] = SFS[i] + sfs_chr[i]   
    
    np.save("popSFS.npy", SFS)