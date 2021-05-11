import tskit, os, sys, tsinfer

if __name__ == '__main__':
    # takes a tree file and infer another tree by using genotype matrix
    
    r_file = sys.argv[1]
    fo, fi = os.path.split(r_file)
    fi = 'i'+fi
    ir_file = os.path.join(fo, fi)
    
    ts = tskit.load(r_file)
    GM = ts.genotype_matrix()
    muts = ts.mutations()
    positions = [next(muts).position for _ in range(GM.shape[0])] 
    L = GM.shape[0]
    seqlen = ts.get_sequence_length()
    with tsinfer.SampleData(sequence_length=seqlen) as sample_data:
        for i in range(L):
            sample_data.add_site(positions[i], GM[i])
    tsinfer.infer(sample_data).simplify().dump(ir_file)