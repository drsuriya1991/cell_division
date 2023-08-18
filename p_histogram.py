
#
# stem cell
#
#########################################################################################
import numpy as np, os
import matplotlib.pyplot as plt
import scipy, sys
import random as ra
from datetime import datetime

#########################################################################################
start_time = datetime.now()
print (start_time)


#---------------------------------------------------------------------------------------------------------------
def histogram(x):
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
    Lagend = ['  horizontal, noise','  vertical, noise','  mixed, noise','  random, noise','  horizontal, no-noise','  vertical, no-noise','  mixed, no-noise','  random, no-noise']
    TitlE = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']

    P_R = ['P=0.1','P=0.5','P=0.8']

    for jj in range(3):
        plt.figure(figsize=(14.0,9.0))
        plt.subplots_adjust(left=0.040, right=0.985, top=0.940, bottom=0.05, hspace = 0.65, wspace = 0.35)
        #plt.title('Histogram of Leaves A, B, S, and A+B+S total nodes')
        num_bins = 150
        for i in range(4):
            #if Suf == 'S': p,q,r,s=1,2,3,4
            #else: p,q,r,s=6,7,8,9
            FS = 15; fs=11; FSS=16
            p,q,r,s=0,1,2,4
            plt.subplot(4,6,6*i+1)
            n, bins, patches = plt.hist(Data[i][jj][:,p]/Samples, num_bins, facecolor='r', alpha=1, weights=np.ones(len(Data[i][jj][:,p]))/len(Data[i][jj][:,p]))
            if i == 0 or 1 or 2 or 3 or 5 or 6 or 7: plt.ylabel('Probability', fontsize=FS);
            if i == 3: plt.xlabel('A Leaves', fontsize=FS);
            plt.title(r'{}$_1$'.format(TitlE[i]), loc='left', fontsize=FSS)
            plt.title(P_R[jj], loc='right', fontsize=fs)
            plt.text(0,max(n)*(0.90),Lagend[i], fontsize=10);
            #plt.title(P_R[jj], loc='left', fontsize=6)
            Xaxis = Data[i][jj][:,p]
            #plt.xticks([0, int(max(Xaxis)/2), max(Xaxis)], fontsize=fs-1)
            plt.xticks([0, 70, 140, 210, 280], fontsize=fs-1)
            plt.yticks([0, round(max(n)/4, 2), round(2*max(n)/4,2), round(3*max(n)/4,2), round(4*max(n)/4,2)], fontsize=fs-1)
        
            plt.subplot(4,6,6*i+2)
            n, bins, patches = plt.hist(Data[i][jj][:,q]/Samples, num_bins, facecolor='b', alpha=1, weights=np.ones(len(Data[i][jj][:,q]))/len(Data[i][jj][:,q]))
            #plt.ylabel('Probability', fontsize=9);
            if i == 3: plt.xlabel('B Leaves', fontsize=FS);
            plt.title(r'{}$_2$'.format(TitlE[i]), loc='left', fontsize=FSS)
            plt.title(P_R[jj], loc='right', fontsize=fs)
            plt.text(0,max(n)*(0.90),Lagend[i], fontsize=10);
            #plt.title(P_R[jj], loc='left', fontsize=6)
            Xaxis = Data[i][jj][:,q]
            #plt.xticks([0, int(max(Xaxis)/2), max(Xaxis)], fontsize=fs-1)
            plt.xticks([0, 70, 140, 210, 280], fontsize=fs-1)
            plt.yticks([0, round(max(n)/4, 2), round(2*max(n)/4,2), round(3*max(n)/4,2), round(4*max(n)/4,2)], fontsize=fs-1)

            plt.subplot(4,6,6*i+3)
            n, bins, patches = plt.hist(Data[i][jj][:,r]/Samples, num_bins, facecolor='g', alpha=1, weights=np.ones(len(Data[i][jj][:,r]))/len(Data[i][jj][:,r]))
            #plt.ylabel('Probability', fontsize=9);
            if i == 3: plt.xlabel('S Leaves', fontsize=FS);
            plt.title(r'{}$_3$'.format(TitlE[i]), loc='left', fontsize=FSS)
            plt.title(P_R[jj], loc='right', fontsize=fs)
            plt.text(0,max(n)*(0.90),Lagend[i], fontsize=10);
            #plt.title(P_R[jj], loc='left', fontsize=6)
            Xaxis = Data[i][jj][:,r]
            #plt.xticks([0, int(max(Xaxis)/2), max(Xaxis)], fontsize=fs-1)
            plt.xticks([0, 100, 200, 300, 400], fontsize=fs-1)
            plt.yticks([0, round(max(n)/4, 2), round(2*max(n)/4,2), round(3*max(n)/4,2), round(4*max(n)/4,2)], fontsize=fs-1)

            
            plt.subplot(4,6,6*i+4)
            n, bins, patches = plt.hist(Data[i+4][jj][:,p]/Samples, num_bins, facecolor='r', alpha=1, weights=np.ones(len(Data[i+4][jj][:,p]))/len(Data[i+4][jj][:,p]))
            #if i == 0 or 1 or 2 or 3 or 5 or 6 or 7: plt.ylabel('Probability', fontsize=FS);
            if i == 3: plt.xlabel('A Leaves', fontsize=FS);
            plt.title(r'{}$_4$'.format(TitlE[i]), loc='left', fontsize=FSS)
            plt.title(P_R[jj], loc='right', fontsize=fs)
            plt.text(0,max(n)*(0.90),Lagend[i+4], fontsize=10);
            #plt.title(P_R[jj], loc='left', fontsize=6)
            Xaxis = Data[i][jj][:,p]
            #plt.xticks([0, int(max(Xaxis)/2), max(Xaxis)], fontsize=fs-1)
            plt.xticks([0, 40, 80, 120, 160], fontsize=fs-1)
            plt.yticks([0, round(max(n)/4, 2), round(2*max(n)/4,2), round(3*max(n)/4,2), round(4*max(n)/4,2)], fontsize=fs-1)
            if jj==0 and i==3: plt.axis([-10,160,0.0,0.175])


            plt.subplot(4,6,6*i+5)
            n, bins, patches = plt.hist(Data[i+4][jj][:,q]/Samples, num_bins, facecolor='b', alpha=1, weights=np.ones(len(Data[i+4][jj][:,q]))/len(Data[i+4][jj][:,q]))
            #plt.ylabel('Probability', fontsize=9);
            if i == 3: plt.xlabel('B Leaves', fontsize=FS);
            plt.title(r'{}$_5$'.format(TitlE[i]), loc='left', fontsize=FSS)
            plt.title(P_R[jj], loc='right', fontsize=fs)
            plt.text(0,max(n)*(0.90),Lagend[i+4], fontsize=10);
            #plt.title(P_R[jj], loc='left', fontsize=6)
            Xaxis = Data[i][jj][:,q]
            #plt.xticks([0, int(max(Xaxis)/2), max(Xaxis)], fontsize=fs-1)
            plt.xticks([0, 40, 80, 120, 160], fontsize=fs-1)
            plt.yticks([0, round(max(n)/4, 2), round(2*max(n)/4,2), round(3*max(n)/4,2), round(4*max(n)/4,2)], fontsize=fs-1)


            plt.subplot(4,6,6*i+6)
            n, bins, patches = plt.hist(Data[i+4][jj][:,r]/Samples, num_bins, facecolor='g', alpha=1, weights=np.ones(len(Data[i+4][jj][:,r]))/len(Data[i+4][jj][:,r]))
            #plt.ylabel('Probability', fontsize=9);
            if i == 3: plt.xlabel('S Leaves', fontsize=FS);
            plt.title(r'{}$_6$'.format(TitlE[i]), loc='left', fontsize=FSS)
            plt.title(P_R[jj], loc='right', fontsize=fs)
            plt.text(0,max(n)*(0.90),Lagend[i+4], fontsize=10);
            #plt.title(P_R[jj], loc='left', fontsize=6)
            Xaxis = Data[i][jj][:,r]
            #plt.xticks([0, int(max(Xaxis)/2), max(Xaxis)], fontsize=fs-1)
            plt.xticks([0, 100, 200, 300, 400], fontsize=fs-1)
            plt.yticks([0, round(max(n)/4, 2), round(2*max(n)/4,2), round(3*max(n)/4,2), round(4*max(n)/4,2)], fontsize=fs-1)
            if jj==1 and i==3: plt.axis([-10,400,0.0,0.33])

        #plt.suptitle('Histogram of Leaves A, B, S, and A+B+S total nodes', fontsize=12)
        plt.tight_layout()
        if jj==0: plt.savefig('histogram_1.pdf')
        if jj==1: plt.savefig('histogram_2.pdf')
        if jj==2: plt.savefig('histogram_3.pdf')
        
    #plt.show()
    
print (histogram(0))
    
    
end_time = datetime.now()
print ('Ends Here||H:M:S||{}'.format(end_time - start_time), '\n')
#######################################################################################################
