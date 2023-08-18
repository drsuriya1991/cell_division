
#
# stem cell 
#
#########################################################################################
import numpy as np, os
import matplotlib.pyplot as plt
import scipy, sys
import random as ra
from matplotlib.lines import Line2D
from datetime import datetime

#########################################################################################
start_time = datetime.now()
print (start_time)

#________________________________________________________________________________________________________
def gene_Idivision(x):
    #P_Value=0.1; Suf='A'; Noise='noise'
    Samples=5000; Time = 5000; Step = 1.0; gen=20
        
    # Final ploting
    # ----------------------------------------------------------------
    #P_Value = 0.1
    DataG, DataT = [], []
    for P_Value in [0.1,0.5,0.8]:
        for Noise in ['noise']:
            for Suf in ['M']:
                Datat = np.genfromtxt('./generation/Time_{}_{}_{}.dat'.format(Noise,Suf,P_Value))
                Datag = []
                for ij in range(gen):
                    Input = np.genfromtxt('./generation/Gen{}_{}_{}_{}.dat'.format(ij,Noise,Suf,P_Value))
                    if np.any(Input > 0): Datag.append(Input)
                    else: print (ij, Suf, Noise, P_Value); break
                Datat[0]=Samples#; Datag[0][0]=Samples
                DataG.append(np.array(Datag))
                DataT.append(Datat)


    # Final Maximum cells ploting
    # ----------------------------------------------------------------
    #P_Value = 0.1
    Data = []
    for P_Value in [0.1,0.5,0.8]:
        for Noise in ['noise']:
            for Suf in ['M']:
                Datag = []
                for ij in range(gen):
                    Input = np.genfromtxt('./generation/Gen{}_{}_{}_{}.dat'.format(ij,Noise,Suf,P_Value))
                    if np.any(Input > 0): Datag.append(max(Input/Samples))
                    else: print (ij, Suf, Noise, P_Value); break
                Data.append(Datag)


    Lagend = ['  mixed, noise', '  mixed, nonoise']
    
    POS = ['P=0.1', 'P=0.5', 'P=0.8']
    
    TITLE = ['A', 'B', 'C',
            'D', 'E', 'F',
            'G', 'H', 'I',
            'J', 'K', 'L']

    fig = plt.figure(figsize=(8.0,8.0))
    #dummie_cax = fig.scatter(gen, gen, c=gen, cmap=plt.cm.nipy_spectral(np.linspace(0, 1, gen)))
    #plt.subplots_adjust(left=0.040, right=0.985, top=0.940, bottom=0.05, hspace = 0.65, wspace = 0.45)
    FS = 15; fs = 11; m = 3; n = 3
    for i, data, datat, datam in zip(list(range(n))*m,DataG, DataT, Data):

        plt.subplot(m,n,i+1)
        #print Index
        Time_Divi = data.sum(axis=0)
        colors = plt.cm.nipy_spectral(np.linspace(0, 1, gen))
        #for j in range(len(data)):
        #    Index1 = list(data[:j].sum(axis=0)).index(0)
            #plt.plot((data[j-1:j].sum(axis=0))[:Index1+1]/float(Samples), color=colors[j],linewidth=0.65) # Unhesh for cummulative distribution
        #    plt.fill_between(range(Index1+1), (data[:j-1].sum(axis=0))[:Index1+1]/float(Samples), (data[:j].sum(axis=0))[:Index1+1]/float(Samples), facecolor=colors[j], alpha=1.0)

        Index = list(Time_Divi).index(0)
        plt.plot(Time_Divi[:Index]/float(Samples), '--g',linewidth=1.0) # Unhesh for cummulative distribution
        
        Death = []
        for ijk in range(len(Time_Divi[:Index])-1):
            Death.append(Time_Divi[ijk+1] - Time_Divi[ijk])
            
            
        for ixy in range(len(Death)):
            if Death[ixy] < 0: Death[ixy] = abs(Death[ixy])
            else: Death[ixy] = 0
        
        #print len(Death)
        #plt.plot(range(len(Death)), np.array(Death)/float(Samples), '--g',linewidth=1.0)
        #plt.plot(range(len(Death)), np.cumsum(np.array(Death))/float(Samples), '--r',linewidth=1.0)
        
        plt.legend(['Live cell', 'Dead'], loc=2, fontsize=fs-4)
        
        Find = Time_Divi[:Index]/float(Samples)
        plt.xticks([0, len(Find)/4, 2*len(Find)/4, 3*len(Find)/4, len(Find)])
        plt.yticks([0, round(max(Find)/4,0), round(2*max(Find)/4,0), round(3*max(Find)/4,0), round(max(Find),0)])

        if i==0: plt.ylabel('Cell count', fontsize=FS)
        if i in [0,1,2]: plt.xlabel('Time', fontsize=FS)

        plt.title(r'{}'.format(TITLE[i]), fontsize=15, loc='left')
        plt.title(r'{}, {}, {}'.format(POS[i], 'mixed', 'noise'), fontsize=fs-1, loc='right')
        plt.xticks(fontsize=fs-1); plt.yticks(fontsize=fs-1)

        plt.subplot(m,n,i+4)
        #print Index
        Time_Divi = data.sum(axis=0)
        colors = plt.cm.nipy_spectral(np.linspace(0, 1, gen))
        for j in range(len(data)):
            Index1 = list(data[:j].sum(axis=0)).index(0)
            #plt.plot((data[j-1:j].sum(axis=0))[:Index1+1]/float(Samples), color=colors[j],linewidth=0.65) # Unhesh for cummulative distribution
            plt.fill_between(range(Index1+1), (data[:j-1].sum(axis=0))[:Index1+1]/float(Samples), (data[:j].sum(axis=0))[:Index1+1]/float(Samples), facecolor=colors[j], alpha=1.0)

        Index = list(Time_Divi).index(0)
        #plt.plot(Time_Divi[:Index]/float(Samples), 'k',linewidth=2.0) # Unhesh for cummulative distribution
        
        Find = Time_Divi[:Index]/float(Samples)
        plt.xticks([0, len(Find)/4, 2*len(Find)/4, 3*len(Find)/4, len(Find)])
        plt.yticks([0, round(max(Find)/4,0), round(2*max(Find)/4,0), round(3*max(Find)/4,0), round(max(Find),0)])

        if i==0: plt.ylabel('Cell count', fontsize=FS)
        if i in [0,1,2]: plt.xlabel('Time', fontsize=FS)

        plt.title(r'{}'.format(TITLE[i+3]), fontsize=15, loc='left')
        plt.title(r'{}, {}, {}'.format(POS[i], 'mixed', 'noise'), fontsize=fs-1, loc='right')
        plt.xticks(fontsize=fs-1); plt.yticks(fontsize=fs-1)



        plt.subplot(m,n,i+7)
        #print Index
        Time_Divi = data.sum(axis=0)
        colors = plt.cm.nipy_spectral(np.linspace(0, 1, gen))
        
        Maxx = []
        for j in range(len(data)):
            Index1 = list(data[:j].sum(axis=0)).index(0)
            Index0 = 0
            Maxx.append(max(data[j]))
            for ii in range(len(data[j])):
                if data[j][ii] > 0: Index0 = ii-1; break
            
            plt.plot(range(len(data[j]))[Index0:Index1+1], (data[j])[Index0:Index1+1]/float(Samples), color=colors[j],linewidth=0.65)

        Index = list(Time_Divi).index(0)
        
        Maxx = max(Maxx)/5000
        
        Find = Time_Divi[:Index]/float(Samples)
        plt.xticks([0, len(Find)/4, 2*len(Find)/4, 3*len(Find)/4, len(Find)])
        plt.yticks([0, round(Maxx/4,0), round(2*Maxx/4,0), round(3*Maxx/4,0), round(Maxx,0)])
    
        if i==0: plt.ylabel('Cell count', fontsize=FS)
        if i in [0,1,2]: plt.xlabel('Time', fontsize=FS)
        
        plt.title(r'{}'.format(TITLE[i+6]), fontsize=15, loc='left')
        plt.title(r'{}, {}, {}'.format(POS[i], 'mixed', 'noise'), fontsize=fs-1, loc='right')
        plt.xticks(fontsize=fs-1); plt.yticks(fontsize=fs-1)

    plt.tight_layout()
    plt.savefig('generation.pdf')
    #plt.show()
    
print (gene_Idivision(0))

end_time = datetime.now()
print ('Ends Here||H:M:S||{}'.format(end_time - start_time), '\n')
#######################################################################################################
