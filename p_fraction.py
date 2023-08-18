
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

#--------------------------------------------------------------------------------------------------------------
def fraction(x):
    P_Range = [0.1,0.5,0.8]; Samples=1.0
    Data1 = [np.genfromtxt('./histogram/dist_{}_{}_{}.dat'.format('noise','A',P_Value)) for P_Value in P_Range]
    Data2 = [np.genfromtxt('./histogram/dist_{}_{}_{}.dat'.format('noise','S',P_Value)) for P_Value in P_Range]
    Data3 = [np.genfromtxt('./histogram/dist_{}_{}_{}.dat'.format('noise','M',P_Value)) for P_Value in P_Range]
    Data4 = [np.genfromtxt('./histogram/dist_{}_{}_{}.dat'.format('noise','R',P_Value)) for P_Value in P_Range]
    Data5 = [np.genfromtxt('./histogram/dist_{}_{}_{}.dat'.format('nonoise','A',P_Value)) for P_Value in P_Range]
    Data6 = [np.genfromtxt('./histogram/dist_{}_{}_{}.dat'.format('nonoise','S',P_Value)) for P_Value in P_Range]
    Data7 = [np.genfromtxt('./histogram/dist_{}_{}_{}.dat'.format('nonoise','M',P_Value)) for P_Value in P_Range]
    Data8 = [np.genfromtxt('./histogram/dist_{}_{}_{}.dat'.format('nonoise','R',P_Value)) for P_Value in P_Range]

    Data = [Data1,Data2,Data3,Data4,Data5,Data6,Data7,Data8]
    Lagend = ['horizontal, noise','vertical, noise','mixed, noise','random, noise','horizontal, no-noise','vertical, no-noise','mixed, no-noise','random, no-noise']

    P_R = ['P=0.1','P=0.5','P=0.8']
    TitlE = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']

    plt.figure(figsize=(8.0,7.0))
    plt.subplots_adjust(left=0.040, right=0.985, top=0.940, bottom=0.05, hspace = 0.65, wspace = 0.35)
    for i in range(9):
        #if Suf == 'S': p,q,r,s=1,2,3,4
        #else: p,q,r,s=6,7,8,9
        SAMPLES=5000; FS=13; fs=11
        plt.subplot(3,3,i+1)
        
        if i < 8:
            ALL_Cells = [(sum(Data[i][ij][:,0])+sum(Data[i][ij][:,1])+sum(Data[i][ij][:,2])) for ij in range(len(P_Range))]
            #print max(list(np.array(ALL_Cells)/SAMPLES))
            MAX = 1.0 #= max(list(np.array(ALL_Cells)/SAMPLES))
            #plt.bar(np.array(P_Range)-0.05, (np.array(ALL_Cells)/SAMPLES)/MAX, width=0.05, color='b')
            if i==0: plt.bar(np.array(P_Range)-0.025, [146.023, 124.8346, 162.5914], width=0.05, color='r')
            elif i==6: plt.bar(np.array(P_Range)-0.025, [128.3078, 100.5184, 92.2438], width=0.05, color='r')
            else: plt.bar(np.array(P_Range)-0.025, [((sum(Data[i][ij][:,0])+sum(Data[i][ij][:,1]))/(SAMPLES*MAX)) for ij in range(len(P_Range))], width=0.05, color='r')
            
            if i==6: plt.bar(np.array(P_Range)+0.025, [218.3366, 236.1694, 256.1526], width=0.05, color='g')
            elif i==7: plt.bar(np.array(P_Range)+0.025, [212.0098, 215.3684, 217.3212], width=0.05, color='g')
            else: plt.bar(np.array(P_Range)+0.025, [sum(Data[i][ij][:,2])/(SAMPLES*MAX) for ij in range(len(P_Range))], width=0.05, color='g')
            
            if i == 0 or i==3 or i==6: plt.ylabel('Cell count', fontsize=FS+2);
            plt.xlabel('Nucleus position', fontsize=FS+2);
            plt.title(Lagend[i], loc='right', fontsize=10.5);
            plt.title(TitlE[i], loc='left', fontsize=18)
            #plt.title(P_R[jj], loc='left', fontsize=6)
            plt.xticks([0.1,0.5,0.8], fontsize=fs);
            #plt.yticks([0,0.25,0.50,0.75,1.00], fontsize=9)
            plt.yticks([0,100.0,200,300], fontsize=fs)
            #if i==0: plt.legend(['Leaves: A+B+S', 'Leaves: A+B', 'Leaves: S'], fontsize=5, loc=2)

        else:
            legend_elements = [Line2D([0], [0], color='r', lw=4, label='Leaves: A+B'),
                               Line2D([0], [0], color='g', lw=4, label='Leaves: S')]
            plt.legend(handles=legend_elements, loc=2)
            #plt.legend(['Leaves: A+B+S', 'Leaves: A+B', 'Leaves: S'], fontsize=10, loc=2)
            plt.axis('off')
        
    #plt.suptitle('Histogram of Leaves A, B, S, and A+B+S total nodes', fontsize=12)
    plt.tight_layout()
    plt.savefig('fraction.pdf')
    #plt.show()

print (fraction(0))


end_time = datetime.now()
print ('Ends Here||H:M:S||{}'.format(end_time - start_time), '\n')
#######################################################################################################
