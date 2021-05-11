import pandas as pd
import re, os

SPops = {0: 'EAS',
 1: 'EAS',
 2: 'EAS',
 3: 'EAS',
 4: 'EAS',
 5: 'EUR',
 6: 'EUR',
 7: 'EUR',
 8: 'EUR',
 9: 'EUR',
 10: 'AFR',
 11: 'AFR',
 12: 'AFR',
 13: 'AFR',
 14: 'AFR',
 15: 'AFR',
 16: 'AFR',
 17: 'AMR',
 18: 'AMR',
 19: 'AMR',
 20: 'AMR',
 21: 'SAS',
 22: 'SAS',
 23: 'SAS',
 24: 'SAS',
 25: 'SAS'}



Pops = {0: 'CHB--Han Chinese in Beijing, China',
 1: 'JPT--Japanese in Tokyo, Japan',
 2: 'CHS--Southern Han Chinese',
 3: 'CDX--Chinese Dai in Xishuangbanna, China',
 4: 'KHV--Kinh in Ho Chi Minh City, Vietnam',
 5: 'CEU--Utah Residents (CEPH) with Northern and Western European Ancestry',
 6: 'TSI--Toscani in Italia',
 7: 'FIN--Finnish in Finland',
 8: 'GBR--British in England and Scotland',
 9: 'IBS--Iberian Population in Spain',
 10: 'YRI--Yoruba in Ibadan, Nigeria',
 11: 'LWK--Luhya in Webuye, Kenya',
 12: 'GWD--Gambian in Western Divisions in the Gambia',
 13: 'MSL--Mende in Sierra Leone',
 14: 'ESN--Esan in Nigeria',
 15: 'ASW--Americans of African Ancestry in SW USA',
 16: 'ACB--African Caribbeans in Barbados',
 17: 'MXL--Mexican Ancestry from Los Angeles USA',
 18: 'PUR--Puerto Ricans from Puerto Rico',
 19: 'CLM--Colombians from Medellin, Colombia',
 20: 'PEL--Peruvians from Lima, Peru',
 21: 'GIH--Gujarati Indian from Houston, Texas',
 22: 'PJL--Punjabi from Lahore, Pakistan',
 23: 'BEB--Bengali from Bangladesh',
 24: 'STU--Sri Lankan Tamil from the UK',
 25: 'ITU--Indian Telugu from the UK'}
popmsg = 'Available Populations:\n'+'\n'.join(['('+str(i)+') '+ Pops[i] for i in range(len(Pops))])

Stats = ['bsfs',
 'SS',
 'Colless',
 'btree',
 'FayH',
 'ZngECmedi',
 'FerL',
 'btreeCmedi',
 'ZngE',
 'TajD',
 'FayHCmedi',
 'FulDCmedi',
 'bsfsCmedi',
 'SSCmedi',
 'FerLCmedi',
 'TajDCmedi',
 'FulD',
 'CollessCmedi']
Stats = sorted(Stats, key=str.casefold)
Stats = {i:Stats[i] for i in range(len(Stats))}
statmsg = 'Available statistics:\n'+'\n'.join(['('+str(i)+') '+ Stats[i] for i in range(len(Stats))])

lookmsg = 'Available methods:\n'+'\n'.join(['(0) Search by a genomic position', '(1) Search by a gene'])

def welcome():
    print(80*'-')
    print('A Tree imbalance measure to detect natural selection')    
    print(80*'-')
    Print('Welcome to our results browser. This tool is designed to report our findings. We have results for all 26 populations of 1000 genomes project for '+str(len(Stats)//2)+' statistics. Our process can be read in tthe paper. Briefly, we calculate genome-wide windowed statistics and then by using a change point detection method we define segments for each genome scan, then calculate how significant the region is. We also have results for median centered statistics. In order to understand how significant a statistic for a population compared to other populations, we centered them around the median. These statistics are under the name of <stat>Cmedi. ')
    
def Print(mystr):
    if all((len(mystr)>80, len(re.findall('\n', mystr))<2)):
        n = len(mystr)
        
        i = 80
        x = mystr[:i]
        while(i<n):
            x = x+'\n'+mystr[i:(i+80)]
            i = i +80
        mystr = x
        
    print(mystr)

def Warn(msg):
    print(36*'-'+'Warning!'+36*'-')
    Print(msg)
    print(36*'-'+'Warning!'+36*'-')  
    
def evalinput(i, fn, msg):
    
    try:
        ret = fn(i)
    except:
        Warn(msg)
        ret = False
        
    return ret

def compare_stats():
    print(80*'-')
    print(popmsg)    
    cont = True
    while(cont):
        inp = input('\nEnter a single population: ')
        fn = lambda i: [int(x) for x in i.replace(' ', '').split(',') if int(x) in range(26)][0]
        msg = 'Choose population by just writing its id'
        ret = evalinput(inp, fn, msg)
        if ret+1:
            cont = False
        pop = ret
            
    print(80*'-')
    print(statmsg)    
    cont = True
    while(cont):
        inp = input('\nEnter multiple statistics (eg: 1,2,3,5): ')
        fn = lambda i: [Stats[int(x)] for x in i.replace(' ', '').split(',') if int(x) in range(len(Stats))]
        msg = 'Seperate selections with comma: 1,2,3,5'
        ret = evalinput(inp, fn, msg)
        if ret:
            cont = False
        stats = ret

    print(80*'-')
    print(lookmsg)    
    cont = True
    while(cont):
        inp = input('\nIndicate your seelction: ')
        fn = lambda i: ['position', 'gene'][int(i)] 
        msg = 'Please try again'
        ret = evalinput(inp, fn, msg)
        if ret:
            cont = False
        lookup = ret       
        
    
    
    cont = True
    print(80*'-')
    while(cont):
        
        if lookup == 'position':
            chrno = int(input('Chromosome? '))
            position = int(input('Genomic location in bp? '))
        elif lookup == 'gene':
            gene = input('Name of the gene?')

        i = 0
        for i, stat in enumerate(stats):
            df = pd.read_csv('gene_scores/'+str(pop)+'_'+stat+'.csv')
            if lookup == 'gene':
                df = df[df['genes'].apply(lambda x: any([i == gene for i in str(x).replace(' ','').split(',')]))]
            if lookup == 'position':
                df = df[(df['Chromosome'] == chrno)&(df['start']<=position)&(df['end']>=position)]
            
            if df.shape[0]<1:
                Warn('Unknown setting!\n')
                break
            else:
                start = df['start'].iloc[0]
                end = df['end'].iloc[0]
                val = df[stat].iloc[0]
                pval = df['p'+stat].iloc[0]
                if pval<0.5:
                    pval = '| p-val: '+ "{:.3e}".format(pval)
                else:
                    pval = '| p-val: 1-'+"{:.3e}".format(1-pval)

                if i == 0:
                    print(80*'-')
                    print('Population: '+Pops[pop])
                    print(80*'-')
                print('Segment region:',str(start)+'--'+str(end)+' ('+str(end-start)+' bp)')
                print('Avarage '+stat+':',"{:.3f}".format(val), pval)
                Print('Genes in the segment: '+str(df['genes'].iloc[0]))
                print(80*'-')
    
        cont = input('\nDo you want to look for another '+lookup+'? (y/n): ')
        if cont.lower()[0] == 'y':
            cont = True
        else:
            cont = False

def compare_pops():
    print(80*'-')
    print(popmsg)    
    cont = True
    while(cont):
        inp = input('\nEnter multiple populations (eg: 1,2,3,5): ')
        fn = lambda i: [int(x) for x in i.replace(' ', '').split(',') if int(x) in range(26)]
        msg = 'Choose populations by just writing their ids (eg: 1,2,3,5)'
        ret = evalinput(inp, fn, msg)
        if ret:
            cont = False
        pops = ret
            
    print(80*'-')
    print(statmsg)    
    cont = True
    while(cont):
        inp = input('\nEnter a single statistic: ')
        fn = lambda i: [Stats[int(x)] for x in i.replace(' ', '').split(',') if int(x) in range(len(Stats))][0]
        msg = 'Choose a stat by just writing its id'
        ret = evalinput(inp, fn, msg)
        if ret:
            cont = False
        stat = ret

    print(80*'-')
    print(lookmsg)    
    cont = True
    while(cont):
        inp = input('\nIndicate your seelction: ')
        fn = lambda i: ['position', 'gene'][int(i)] 
        msg = 'Please try again'
        ret = evalinput(inp, fn, msg)
        if ret:
            cont = False
        lookup = ret       
    
    cont = True
    print(80*'-')
    while(cont):

        if lookup == 'position':
            chrno = int(input('Chromosome? '))
            position = int(input('Genomic location in bp? '))
        elif lookup == 'gene':
            gene = input('Name of the gene?')

        print(80*'-')
        for pop in pops:
            df = pd.read_csv('gene_scores/'+str(pop)+'_'+stat+'.csv')
            if lookup == 'gene':
                df = df[df['genes'].apply(lambda x: any([i == gene for i in str(x).replace(' ','').split(',')]))]
            if lookup == 'position':
                df = df[(df['Chromosome'] == chrno)&(df['start']<=position)&(df['end']>=position)]
            
            if df.shape[0]<1:
                Warn('Unknown setting!')
                break
            else:
                start = df['start'].iloc[0]
                end = df['end'].iloc[0]
                val = df[stat].iloc[0]
                pval = df['p'+stat].iloc[0]
                if pval<0.5:
                    pval = '| p-val: '+ "{:.3e}".format(pval)
                else:
                    pval = '| p-val: 1-'+"{:.3e}".format(1-pval)
                print('Population: '+Pops[pop])
                print('Segment region:',str(start)+'--'+str(end)+' ('+str(end-start)+' bp)')
                print('Avarage '+stat+':',"{:.3f}".format(val), pval)
                Print('Genes in the segment: '+str(df['genes'].iloc[0]))
                print(80*'-')
            
        cont = input('\nDo you want to look for another '+lookup+'? (y/n): ')
        if cont.lower()[0] == 'y':
            cont = True
        else:
            cont = False
def browse_segment():
    
    print(80*'-')
    print(statmsg)    
    cont = True
    while(cont):
        inp = input('\nEnter a single statistic: ')
        fn = lambda i: [Stats[int(x)] for x in i.replace(' ', '').split(',') if int(x) in range(len(Stats))][0]
        msg = 'Choose a stat by just writing its id'
        ret = evalinput(inp, fn, msg)
        if ret:
            cont = False
        stat = ret
    
    cont = True
    
    print(80*'-')
    while(cont):
       
        segment = input('Enter the segment: ')
        pop = int(segment[1:segment.find('c')])
        
        df = pd.read_csv('gene_scores/'+str(pop)+'_'+stat+'.csv')
        df = df[df['s'+stat] == segment]
        if df.shape[0]<1:
            Warn('Unknown setting!')
        else:     
            print(80*'-')
            start = df['start'].iloc[0]
            end = df['end'].iloc[0]
            val = df[stat].iloc[0]
            pval = df['p'+stat].iloc[0]
            if pval<0.5:
                pval = '| p-val: '+ "{:.3e}".format(pval)
            else:
                pval = '| p-val: 1-'+"{:.3e}".format(1-pval)
            print('Population: '+Pops[pop])
            print('Segment region:',str(start)+'--'+str(end)+' ('+str(end-start)+' bp)')
            print('Avarage '+stat+':',"{:.3f}".format(val), pval)
            Print('Genes in the segment: '+str(df['genes'].iloc[0]))
            print(80*'-')
        
        cont = input('\nDo you want to look for another segment? (y/n): ')
        if cont.lower()[0] == 'y':
            cont = True
        else:
            cont = False

if __name__ == '__main__':
    
    welcome()
    
    methods = ['Compare statistics for a single population.', 
               'Compare populations for the results of a single statistic.',
               'Browse genome segments for a statistic spesific to a population.(see plots/)']
    
    cont = True
    while(cont):
           
        print('\n'.join(['('+str(i)+') '+ methods[i] for i in range(len(methods))]))
        inp = input('How do you want to browse the results?\n')
        if inp == '0':
            compare_stats()
            cont = False
        elif inp == '1':
            compare_pops()
            cont = False
        elif inp == '2':
            browse_segment()
            cont = False
        else:
            print(46*'-'+'Warning!'+46*'-')
            print('Available options:'+ ' '.join([str(i) for i in range(len(methods))]))
            print(46*'-'+'Warning!'+46*'-')




