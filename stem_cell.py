#
# stem cell tree
#
#########################################################################################
import numpy as np, os
import matplotlib.pyplot as plt
from scipy import *
import random as ra
import scipy.integrate
from sympy import plot_implicit, cos, sin, symbols, Eq 
import seaborn as sns
import networkx as nx
from networkx.drawing.nx_agraph import write_dot, graphviz_layout
from multiprocessing import Pool
from datetime import datetime
#########################################################################################
start_time = datetime.now()
print start_time

def Stem_Cell((P_Value,Ind)):
	G = nx.Graph()
	fate={'s':1,'a':0,'b':0}

	#parameters
	vol=2 #mm3
	Len=1 #unit
	pos=np.linspace(0,0.5,20)# position of the nucleus

	#motif deterministic implememtation
	def toggle(ab, t,a1,a2,b1,b2,S,Mux,Muy,m, n):
		"""Right hand side for toggle ODEs."""
		a, b = ab
		return np.array([a1*a**n/ (S**n + a**n)+ b1/ (S**m + b**m)-Mux*a,
			             a2*b**m/ (S**m + b**m)+ b2/ (S**n + a**n)-Muy*b])

	# Parameters
	a1=0.60; b1=.3; Mux=0.3
	a2=0.6; b2=.3; Muy=0.3
	a=a1
	b=b1

	#a1=a2 =0.5
	#b1=b2 = 0.5
	n =3.0
	m =3.0
	S=1
	#k1=k2=1
	args = (a1,a2,b1,b2,S,Mux,Muy,m, n)

	# Initial condition
	ab0 = np.array([2, 1.5])

	# Solve
	trials=50
	fixed_pts=[]
	d=[]#()
	#dict['y']=[]
	def type(f):
		a,b=np.round(f[0],3),np.round(f[1],3)
		if (a,b)==d[0]:
			return 'b'
		if (a,b)==d[1]:
			return 's'
		else:
			return 'a'
		
	for i in range(0,trials):
		t = np.linspace(0, 200, 20000)
		ab0=np.random.uniform(0,4,2)
		ab = scipy.integrate.odeint(toggle, ab0, t, args=args)
		x=[ab[i][0] for i in range(0,len(ab))]
		y=[ab[i][1] for i in range(0,len(ab))]
		#plt.figure(0)
		#plt.plot(t,x,color='b')
		#plt.plot(t,y,color='g')
		#plt.legend(('Th1','Th2'))
		#plt.xlabel('time')
		#plt.ylabel('conc')
		fixed_full=ab[-1]
		fixed=np.round(ab[-1][0],3),np.round(ab[-1][1],3)
		if fixed not in d:
			fixed_pts.append(fixed_full)
		   #d.add((np.round(ab[-1][0],3),np.round(ab[-1][1],3)))
			d.append(((np.round(ab[-1][0],3),np.round(ab[-1][1],3))))
		   
			
	d=sorted([i for i in d])
	Stem=d[1]        
	#stochiometric matrix
	#V=[[1,0],[-1,0],[0,1],[0,-1]]
	C=np.zeros(2)
	C[0]=np.round(Stem[0]*vol)
	C[1]=np.round(Stem[1]*vol)

	#assuming volume of the cell increases to V=200
	V=2
	X=np.array(C)
	K=[[1,0],[-1,0],[0,1],[0,-1]]
	propensity=[V*((a*X[0]**n)/(1+X[0]**n)+(b/(1+X[1]**m))),Mux*X[0],V*((a*X[1]**m)/(1+X[1]**m)+(b/(1+X[0]**n))),Mux*X[1]]
	Y=[X]
	S=np.zeros(4)
	S[0]=V*((a*X[0]**n)/(1+X[0]**n)+(b/(1+X[1]**m)))
	S[1]=Mux*X[0]
	S[2]=V*((a*X[1]**m)/(1+X[1]**m)+(b/(1+X[0]**n)))
	S[3]=Mux*X[1]
	div=0
	c=np.log(2)/10
	V0=2
	v_final=2*V0
	t=0
	time=[t]
	i=0
	while V<v_final:
		i+=1
		S[0]=V*((a*X[0]**n)/(1+X[0]**n)+(b/(1+X[1]**m)))
		S[1]=Mux*X[0]
		S[2]=V*((a*X[1]**m)/(1+X[1]**m)+(b/(1+X[0]**n)))
		S[3]=Mux*X[1]
		As=np.array(S[0])+np.array(S[2])
		Aq=np.array(S[1])+np.array(S[3])
		alpha=As/Aq
		r=np.random.uniform(0,1)
		beta=c*np.log(r)/Aq
		tau_next= (1/(As+Aq))*np.log(1/r)#scipy.special.lambertw(alpha*np.exp(alpha+beta))-alpha-beta
		t+=tau_next
		V=V0*np.exp(c*t)
		S[0]=V*((a*X[0]**n)/(1+X[0]**n)+(b/(1+X[1]**m)))
		S[2]=V*((a*X[1]**m)/(1+X[1]**m)+(b/(1+X[0]**n)))
		r=np.random.uniform(0,1)
		cum_sum=[sum(np.array(S[0:i+1]))/sum(np.array(S)) for i in range(0,len(S))]
		event=int(sum([1 for i in cum_sum if r>i]))
		X=X+np.array(K[event])
		time.append(t)
		Y.append(X)
	#plt.figure(1)
	#S1=[Y[i][0] for i in range(0,len(Y))]
	#S2=[Y[i][1] for i in range(0,len(Y))]
	#plt.plot(S1)
	#plt.plot(S2)
	#plt.legend(['conc of x','conc of y'])
	#plt.xlabel('time')
	#plt.ylabel('number of molecules')
	t_divide=time[-1]
	#plt.
	#plt.plot(S3)

	pos=P_Value	#nucleus position before division
	#print pos 
	#normalize  conc of a and b based on volume after stochastic partitioning
	#normalize nuclear position
	#cell=[time for next divsion,conc of A at division,conc of B at division ,volume at divison,pos normalized]
	#symmetric or assymetric div
	#0 is symmetric cut else 1 is asymmetric cut
	vol_cell_1 = pos*v_final
	vol_cell_2 = (1-pos)*v_final
	A_final=X[0]
	B_final=X[1]
	A_1,B_1,A_2,B_2=0,0,0,0
	for i in range(0,int(A_final)):
		r=np.random.uniform(0,1,1)
		if r<vol_cell_1/(vol_cell_1+vol_cell_2):
			A_1+=1
		else:
			A_2+=1
	for i in range(0,int(B_final)):
		 r=np.random.uniform(0,1,1)
		 if r<vol_cell_1/(vol_cell_1+vol_cell_2):
			B_1+=1
		 else:
			B_2+=1
		

	random_cut=np.random.uniform(0,1,1)
	#symmetry=0+int(random_cut>pos)
	symmetry=ra.choice([0,1]) ##########################################################################################################################################################
	vol_cell_1 = pos*v_final*(symmetry==1)+ (v_final/2)*(symmetry==0)  #2 is v_final
	cell_1_time_for_next_div = (1/c)*(np.log(v_final/vol_cell_1)) + t_divide
	cell_1_conc_A= A_1/vol_cell_1
	cell_1_conc_B= B_1/vol_cell_1
	pos_norm_cell_1= pos*int(symmetry==0)+int(symmetry==1)*0.5*int(pos==0.5)+ int(symmetry==1)*int(pos<0.5)*0.5 + int(symmetry==1)*int(pos>0.5)*((3*(pos)-1)/2*pos)
	vol_cell_2 = (1-pos)*v_final*(symmetry==1)+ (v_final/2)*(symmetry==0)
	cell_2_time_for_next_div =(1/c)*(np.log(v_final/vol_cell_2))+ t_divide
	cell_2_conc_A= A_2/vol_cell_2
	cell_2_conc_B= B_2/vol_cell_2
	pos_norm_cell_2= pos*int(symmetry==0)+int(symmetry==1)*0.5*int(pos==0.5)+ int(symmetry==1)*int(pos>0.5)*0.5 + int(symmetry==1)*int(pos<0.5)*(pos/(2*(1-pos)))


	noise_1 = min(pos_norm_cell_1,1-pos_norm_cell_1)
	noise_2 = min(pos_norm_cell_2,1-pos_norm_cell_2)

	pos_norm_cell_1 = pos_norm_cell_1 + np.random.uniform(-noise_1/2.0,noise_1/2.0) #################################################################################
	pos_norm_cell_2 = pos_norm_cell_2 + np.random.uniform(-noise_2/2.0,noise_2/2.0) #################################################################################
	
	#print pos_norm_cell_1
	#print pos_norm_cell_2

	test= np.linspace(0, 100, 2000)
	#cell1 fate
	ab1 = scipy.integrate.odeint(toggle, [cell_1_conc_A,cell_1_conc_B], test, args=args)
	xx_1=[ab1[i][0] for i in range(0,len(ab1))]
	yy_1=[ab1[i][1] for i in range(0,len(ab1))]
	final_fate_1=[xx_1[-1],yy_1[-1]]

	#cell1 fate
	ab2 = scipy.integrate.odeint(toggle, [cell_2_conc_A,cell_2_conc_B], test, args=args)
	xx_2=[ab2[i][0] for i in range(0,len(ab2))]
	yy_2=[ab2[i][1] for i in range(0,len(ab2))]
	final_fate_2=[xx_2[-1],yy_2[-1]]
	#print(symmetry)
	#print(final_fate_1,final_fate_2)
	#print(pos_norm_cell_1,pos_norm_cell_2)

    #FATE1 = type(final_fate_1); FATE2 = type(final_fate_2)
	with open(os.path.join('./n3_m3/noise/symmR/Tree_{}_{}.dat'.format(P_Value,Ind)), 'w') as f_series0:
		f_series0.write('#%9s %10s %10s %10s %10s %10s %10s %10s %10s\n'% ('Death', 'Birth', 'Nuclear', 'CellAc', 'CellBc', 'Volume', 'Parent', 'Daught', 'Type'))
		G.add_node(1)
		G.nodes[1]['type']='s'
		G.nodes[1]['parent']=None
		G.nodes[1]['Time']=0  # I have added
		details_cell_1=[cell_1_time_for_next_div,t_divide,pos_norm_cell_1,cell_1_conc_A,cell_1_conc_B,vol_cell_1,1,2,type(final_fate_1)]
		details_cell_2=[cell_2_time_for_next_div,t_divide,pos_norm_cell_2,cell_2_conc_A,cell_2_conc_B,vol_cell_2,1,3,type(final_fate_2)]

		for i in range(len(details_cell_1)):
			if i != len(details_cell_1)-1: f_series0.write('%10.4f '% (details_cell_1[i]))
			else: f_series0.write('%10s\n'% (details_cell_1[i]))

		for i in range(len(details_cell_2)):
			if i != len(details_cell_2)-1: f_series0.write('%10.4f '% (details_cell_2[i]))
			else: f_series0.write('%10s\n'% (details_cell_2[i]))

		G.add_nodes_from([2,3])
		G.nodes[2]['parent']=1
		G.nodes[3]['parent']=1
		G.nodes[2]['type']=type(final_fate_1)
		fate[G.nodes[2]['type']]+=1
		G.nodes[3]['type']=type(final_fate_2)
		fate[G.nodes[3]['type']]+=1
		G.add_edge(1,2)
		G.add_edge(1,3)
		#time,pos,a,b,vol,parent
		list_nodes=[]
		nodes=3
		if type(final_fate_1)=='s':
			list_nodes.append(details_cell_1)
		if type(final_fate_2)=='s':
			list_nodes.append(details_cell_2)
		list_nodes=sorted(list_nodes)


		###########################################################################################################333333#################################
		#######################333##############################33333########################################################################################
		with open(os.path.join('./n3_m3/noise/symmR/Cut_{}_{}.dat'.format(P_Value,Ind)), 'w') as f_series1:
			f_series1.write('%4s %4s %4s %4s %4s %4s %4s %4s %4s %4s %4s\n'% ('Node','SY1', 'SS1', 'SA1', 'SB1', 'AB1','SY2', 'SS2', 'SA2', 'SB2', 'AB2'))
			SY1=0; SS1=0; SA1=0; SB1=0; AB1=0; SY2=0; SS2=0; SA2=0; SB2=0; AB2=0
			######################################################################################################
			if symmetry==0:
				SY1 = 1
				if type(final_fate_1)=='s' and type(final_fate_2)=='s' or type(final_fate_2)=='s' and type(final_fate_1)=='s': SS1=1
				elif type(final_fate_1)=='s' and type(final_fate_2)=='a' or type(final_fate_2)=='s' and type(final_fate_1)=='a': SA1=1
				elif type(final_fate_1)=='s' and type(final_fate_2)=='b' or type(final_fate_2)=='s' and type(final_fate_1)=='b': SB1=1
				else: AB1=1

			if symmetry==1:
				SY2 = 1
				if type(final_fate_1)=='s' and type(final_fate_2)=='s' or type(final_fate_2)=='s' and type(final_fate_1)=='s': SS2=1
				elif type(final_fate_1)=='s' and type(final_fate_2)=='a' or type(final_fate_2)=='s' and type(final_fate_1)=='a': SA2=1
				elif type(final_fate_1)=='s' and type(final_fate_2)=='b' or type(final_fate_2)=='s' and type(final_fate_1)=='b': SB2=1
				else: AB2=1

			#G.nodes[1]['SY1']=SY1; G.nodes[1]['SS1']=SS1; G.nodes[1]['SA1']=SA1; G.nodes[1]['SB1']=SB1; G.nodes[1]['AB1']=AB1
			#G.nodes[1]['SY2']=SY2; G.nodes[1]['SS2']=SS2; G.nodes[1]['SA2']=SA2; G.nodes[1]['SB2']=SB2; G.nodes[1]['AB2']=AB2
			f_series1.write('%4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f\n'% (1, SY1, SS1, SA1, SB1, AB1, SY2, SS2, SA2, SB2, AB2))

			#print 1, '|', SY1, SS1, SA1, SB1, AB1, '|', SY2, SS2, SA2, SB2, AB2
			#print details_cell_1
			#print details_cell_2
			while list_nodes!=[] and nodes<1000:
				new_div=list_nodes.pop(0)
				t_divide=new_div[0]
				pos=new_div[2]
				a_conc=new_div[3]
				b_conc=new_div[4]
				vol=new_div[5]
				a_n=a_conc*vol
				b_n=b_conc*vol
				node_no=new_div[-2]
				New_Node = node_no
				#print node_no
				t=0
				X=a_n,b_n
				Y=[X]
				V=vol
				while V<v_final:
					S[0]=vol*((a1*X[0]**n)/(1+X[0]**n)+(b1/(1+X[1]**m)))
					S[1]=Mux*a_conc
					S[2]=vol*((a1*X[1]**m)/(1+X[1]**m)+(b1/(1+X[0]**n)))
					S[3]=Mux*b_conc
					As=np.array(S[0])+np.array(S[2])
					Aq=np.array(S[1])+np.array(S[3])
					alpha=As/Aq
					r=np.random.uniform(0,1)
					beta=c*np.log(r)/Aq
					tau_next= (1/(As+Aq))*np.log(1/r)#scipy.special.lambertw(alpha*np.exp(alpha+beta))-alpha-beta
					t+=tau_next
					V=vol*np.exp(c*t)
					S[0]=V*((a1*X[0]**n)/(1+X[0]**n)+(b1/(1+X[1]**m)))
					S[2]=V*((a1*X[1]**m)/(1+X[1]**m)+(b1/(1+X[0]**n)))
					r=np.random.uniform(0,1)
					cum_sum=[sum(np.array(S[0:i+1]))/sum(np.array(S)) for i in range(0,len(S))]
					event=int(sum([1 for i in cum_sum if r>i]))
					X=X+np.array(K[event])
					time.append(t)
					Y.append(X)
				#t_divide=time[-1]
				a_n_final=Y[-1][0]
				b_n_final=Y[-1][1]
		        
				#cellfate after partitioning
				random_cut=np.random.uniform(0,1,1)
				#symmetry=0+int(random_cut>pos)
				symmetry=ra.choice([0,1]) ###########################################################################################################################################################
				cell_1_time_for_next_div = (1/c)*(np.log(v_final/pos*(v_final)))+ t_divide
				vol_cell_1 = pos*v_final*(symmetry==1)+ (v_final/2)*(symmetry==0) #2 is v_final
				cell_1_conc_A= A_1/vol_cell_1
				cell_1_conc_B= B_1/vol_cell_1
				pos_norm_cell_1= pos*int(symmetry==0)+int(symmetry==1)*0.5*int(pos==0.5)+ int(symmetry==1)*int(pos<0.5)*0.5 + int(symmetry==1)*int(pos>0.5)*((3*(pos)-1)/2*pos)
				cell_2_time_for_next_div =(1/c)*(np.log(v_final/(1-pos)*v_final)) + t_divide
				vol_cell_2 = (1-pos)*v_final*(symmetry==1)+ (v_final/2)*(symmetry==0)
				cell_2_conc_A= A_2/vol_cell_2
				cell_2_conc_B= B_2/vol_cell_2
				pos_norm_cell_2= pos*int(symmetry==0)+int(symmetry==1)*0.5*int(pos==0.5)+ int(symmetry==1)*int(pos>0.5)*0.5 + int(symmetry==1)*int(pos<0.5)*(pos/(2*(1-pos)))

				noise_1 = min(pos_norm_cell_1,1-pos_norm_cell_1)
				noise_2 = min(pos_norm_cell_2,1-pos_norm_cell_2)

				pos_norm_cell_1 = pos_norm_cell_1 + np.random.uniform(-noise_1/2.0,noise_1/2.0) #################################################################################
				pos_norm_cell_2 = pos_norm_cell_2 + np.random.uniform(-noise_2/2.0,noise_2/2.0) #################################################################################


				#print pos_norm_cell_1
				#print pos_norm_cell_2


				test= np.linspace(0, 100, 2000)
				#cell1 fate
				ab1 = scipy.integrate.odeint(toggle, [cell_1_conc_A,cell_1_conc_B], test, args=args)
				xx_1=[ab1[i][0] for i in range(0,len(ab1))]
				yy_1=[ab1[i][1] for i in range(0,len(ab1))]
				final_fate_1=[xx_1[-1],yy_1[-1]]
				fate_1=type(final_fate_1)

				#cell1 fate
				ab2 = scipy.integrate.odeint(toggle, [cell_2_conc_A,cell_2_conc_B], test, args=args)
				xx_2=[ab2[i][0] for i in range(0,len(ab2))]
				yy_2=[ab2[i][1] for i in range(0,len(ab2))]
				final_fate_2=[xx_2[-1],yy_2[-1]]
				fate_2=type(final_fate_2)

				SY1=0; SS1=0; SA1=0; SB1=0; AB1=0; SY2=0; SS2=0; SA2=0; SB2=0; AB2=0
				######################################################################################################
				if symmetry==0:
					SY1 = 1
					if type(final_fate_1)=='s' and type(final_fate_2)=='s' or type(final_fate_2)=='s' and type(final_fate_1)=='s': SS1=1
					elif type(final_fate_1)=='s' and type(final_fate_2)=='a' or type(final_fate_2)=='s' and type(final_fate_1)=='a': SA1=1
					elif type(final_fate_1)=='s' and type(final_fate_2)=='b' or type(final_fate_2)=='s' and type(final_fate_1)=='b': SB1=1
					else: AB1=1

				if symmetry==1:
					SY2 = 1
					if type(final_fate_1)=='s' and type(final_fate_2)=='s' or type(final_fate_2)=='s' and type(final_fate_1)=='s': SS2=1
					elif type(final_fate_1)=='s' and type(final_fate_2)=='a' or type(final_fate_2)=='s' and type(final_fate_1)=='a': SA2=1
					elif type(final_fate_1)=='s' and type(final_fate_2)=='b' or type(final_fate_2)=='s' and type(final_fate_1)=='b': SB2=1
					else: AB2=1

				#G.nodes[New_Node]['SY1']=SY1; G.nodes[New_Node]['SS1']=SS1; G.nodes[New_Node]['SA1']=SA1; G.nodes[New_Node]['SB1']=SB1; G.nodes[New_Node]['AB1']=AB1
				#G.nodes[New_Node]['SY2']=SY2; G.nodes[New_Node]['SS2']=SS2; G.nodes[New_Node]['SA2']=SA2; G.nodes[New_Node]['SB2']=SB2; G.nodes[New_Node]['AB2']=AB2
				#print New_Node, '|', SY1, SS1, SA1, SB1, AB1, '|', SY2, SS2, SA2, SB2, AB2
				f_series1.write('%4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f\n'% (New_Node, SY1, SS1, SA1, SB1, AB1, SY2, SS2, SA2, SB2, AB2))


				A_final=X[0]
				B_final=X[1]
				A_1,B_1,A_2,B_2=0,0,0,0
				for i in range(0,int(A_final)):
					r=np.random.uniform(0,1,1)
					if r<vol_cell_1/(vol_cell_1+vol_cell_2):
						A_1+=1
					else:
						A_2+=1
				for i in range(0,int(B_final)):
					r=np.random.uniform(0,1,1)
					if r<vol_cell_1/(vol_cell_1+vol_cell_2):
						B_1+=1
					else:
						B_2+=1
				p=nodes+1
				q=nodes+2

				details_cell_1=[cell_1_time_for_next_div,t_divide,pos_norm_cell_1,cell_1_conc_A,cell_1_conc_B,vol_cell_1,node_no,p,fate_1]
				details_cell_2=[cell_2_time_for_next_div,t_divide,pos_norm_cell_2,cell_2_conc_A,cell_2_conc_B,vol_cell_2,node_no,q,fate_2]
		        
				#print details_cell_1
				#print details_cell_2

				for i in range(len(details_cell_1)):
					if i != len(details_cell_1)-1: f_series0.write('%10.4f '% (details_cell_1[i]))
					else: f_series0.write('%10s\n'% (details_cell_1[i]))

				for i in range(len(details_cell_2)):
					if i != len(details_cell_2)-1: f_series0.write('%10.4f '% (details_cell_2[i]))
					else: f_series0.write('%10s\n'% (details_cell_2[i]))

				G.add_nodes_from([p,q])
				G.nodes[p]['parent']=node_no
				G.nodes[q]['parent']=node_no
				G.nodes[p]['type']=type(final_fate_1)
				fate[G.nodes[p]['type']]+=1
				G.nodes[q]['type']=type(final_fate_2)
				fate[G.nodes[q]['type']]+=1
				G.add_edge(node_no,p)
				G.add_edge(node_no,q)
				nodes+=2
				#time,pos,a,b,vol,parent
				if type(final_fate_1)=='s':
					list_nodes.append(details_cell_1)
				if type(final_fate_2)=='s':
					list_nodes.append(details_cell_2)
				list_nodes=sorted(list_nodes)

			color_map=[]
			for i in G.nodes:
				if G.nodes[i]['type']=='s':
					color_map.append('yellow')
				if G.node[i]['type']=='a':
					color_map.append('red')
				if G.node[i]['type']=='b':
					color_map.append('blue')
			#plt.figure(2)
	        
			#with open(os.path.join('./data','Tree_{}_{}.dat'.format(P_Value,Ind)), 'w') as f_series:	
	        
	        
			#nx.draw(G,node_color = color_map,with_labels = True)
			#plt.show()    
			#add to list
	        
			Color=nx.get_node_attributes(G,'type')
	
	#nx.draw(G,node_color = color_map,with_labels = True)
	#plt.show()
	
	#add edges and update graph and repeat
	
	#volume dependent gillespie
	#stochastic partitioning
	#normalize values and repeat
    
	#print Ind

#print Stem_Cell((0.1, 0))


if __name__ == '__main__':
    N_Proces = 14 ; Samples=5000; P_Range = [0.1,0.3,0.5,0.8]; P_Value = []
    for Value in P_Range: 
    	for i in range(Samples):    P_Value.append(Value)
    Pool(N_Proces).map(Stem_Cell, zip(P_Value, range(Samples)*len(P_Range)))    

end_time = datetime.now()
print 'PROGRAM IS COMPLETED||TOTAL TIME TAKEN||Duration||H:M:S||{}'.format(end_time - start_time), '\n'
#######################################################################################################


#while loop it repeat till the tree ends in terminal differnentiated cells and for fifferent cell positions
#values = [val_map.get(node, 0.25) for node in G.nodes() if G.nodes[node]['type']=='s']
#val_map={'s':0,'a':0.5,'b':1}
# nx.draw(G, cmap=plt.get_cmap('viridis'), node_color=values, with_labels=True, font_color='white')
#plt.show()
