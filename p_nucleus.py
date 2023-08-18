
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


#---------------------------
def nucleus(x):
    P_Range = [0.1,0.5,0.8]; Samples=1.0
    Data1 = [np.genfromtxt('fraction/Cut_{}_{}_{}.dat'.format('noise','M',P_Value)) for P_Value in P_Range]
    Data3 = [np.genfromtxt('fraction/Cut_{}_{}_{}.dat'.format('nonoise','M',P_Value)) for P_Value in P_Range]
            
    Lagend = ['mixed, noise','mixed, no-noise','symmetry, no-noise']

    P_R = ['P=0.1','P=0.5','P=0.8']
    TitlE = ['I', 'J', 'K', 'L', 'M', 'N']

    plt.figure(figsize=(8.0,4.0))
    plt.subplots_adjust(left=0.040, right=0.985, top=0.940, bottom=0.05, hspace = 0.65, wspace = 0.35)
    
    Data = [Data1,Data3]
    for i in range(2):
        #if Suf == 'S': p,q,r,s=1,2,3,4
        #else: p,q,r,s=6,7,8,9
        SAMPLES=5000; FS=14; fs=12
        Total_Symmetric_Division = np.array([sum(Data[i][ij][:,0]) for ij in range(len(P_Range))])/SAMPLES
        Symmetric_Division = np.array([sum(Data[i][ij][:,1])+sum(Data[i][ij][:,4]) for ij in range(len(P_Range))])/SAMPLES
        plt.subplot(1,2,i+1)
        plt.bar(np.array(P_Range)-0.025, Total_Symmetric_Division, width=0.05, color='r')
        plt.bar(np.array(P_Range)+0.025, Symmetric_Division, width=0.05, color='b')
        plt.ylabel('Division count', fontsize=FS+3);
        plt.xlabel('Nucleus position', fontsize=FS+3);
        plt.title(Lagend[i], loc='right', fontsize=11);
        plt.title(TitlE[i], loc='left', fontsize=17)
        #plt.title(P_R[jj], loc='left', fontsize=6)
        plt.xticks([0.1,0.5,0.8], fontsize=fs);
        #plt.yticks([0,0.25,0.50,0.75,1.00], fontsize=9)
        plt.yticks([0,50,100,150,200], fontsize=fs)
        plt.axis([0,0.9,0,225])
        plt.legend(['Vertical division count', 'Symmetric fate'], fontsize=9, loc=2)

    #plt.suptitle('Histogram of Leaves A, B, S, and A+B+S total nodes', fontsize=12)
    plt.tight_layout()
    plt.savefig('nucleus.pdf')
    #plt.show()

print (nucleus(0))

end_time = datetime.now()
print ('Ends Here||H:M:S||{}'.format(end_time - start_time), '\n')
#######################################################################################################
