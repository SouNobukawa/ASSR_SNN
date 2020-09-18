from brian2 import *
import numpy
import matplotlib
import sys
import networkx
import random

args = sys.argv


n_size=12000
V_L=-70*mV
V_E=0*mV
V_I=-80*mV


EPSPtoW=100.

sigma_l = 1.0
mu_l=log(0.2)+sigma_l*sigma_l
small_world_p=float(args[1]);

random_seed=int(args[2]);

seed(random_seed)



eqs = '''
dv/dt  = -ge/ms*(v-V_E)-gi/ms*(v-V_I)-(v-V_L)/(tau_mem*ms)  :volt
dge/dt = -ge/(2*ms)                : 1
dgi/dt = -gi/(2*ms)               : 1
tau_mem: 1
'''
P = NeuronGroup(int(n_size), eqs, threshold='v>-50*mV', reset='v=-60*mV', refractory=1*ms)

#P_ex = PoissonGroup(n_size, 0.0*Hz)

#corresponding 40Hz  1/(25)*1000
stimulus = TimedArray(np.tile([1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 10000)*Hz, dt=1.*ms)


#corresponding 83.3Hz  1/(12)*1000
#stimulus = TimedArray(np.tile([1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 10000)*Hz, dt=1.*ms)


#corresponding 100Hz
#stimulus = TimedArray(np.tile([1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 10000)*Hz, dt=1.*ms)

#corresponding 166.6 Hz   1/(6)*1000
#stimulus = TimedArray(np.tile([1.0, 0, 0, 0, 0, 0], 10000)*Hz, dt=1.*ms)

#corresponding 142.8 Hz   1/(7)*1000
#stimulus = TimedArray(np.tile([1, 0, 0, 0, 0, 0, 0], 10000)*Hz, dt=1.*ms)



#stimulus = TimedArray(np.tile([0.3, 0.], 1)*Hz, dt=100.*ms)
P_ex = PoissonGroup(n_size, rates='Lambda*stimulus(t)')

P.v= -70*mV

    

    
Pe = P[:int(n_size*3/4)]
Pi = P[int(n_size*3/4):]
Pe.tau_mem=20
Pi.tau_mem=10
Pe_s=Pe[:int(int(n_size*3/4)*1.0)]
Cee_s = Synapses(Pe_s, Pe_s, model="""w:1
                              p:1""",
                              on_pre='ge_post+=w*(rand()<p)')
Cee = Synapses(Pe, Pe, model="""w:1
                              p:1""",
                              on_pre='ge_post+=w*(rand()<p)')


Threshold_N=9.;
temp=numpy.random.lognormal(mu_l, sigma_l,len(Pe)*int(len(Pe)/10.));
for i in range(0,len(temp)):
    while temp[i]>20:
      temp[i]=numpy.random.lognormal(mu_l, sigma_l)

#temp_s=numpy.array(filter(lambda x: x > Threshold_N, temp));
temp_s=temp[temp > Threshold_N];
#temp_w=numpy.array(filter(lambda x: x <= Threshold_N, temp));
temp_w=temp[temp <= Threshold_N];


print(len(temp_s));
print(float(len(temp_s))/float(len(temp))*int(n_size*3/4));
print(len(Pe_s));
print(2*int(len(temp_s)/len(Pe_s)));


#Strong weights
sn_s=networkx.watts_strogatz_graph(len(Pe_s),2*int(len(temp_s)/len(Pe_s)),small_world_p,random_seed);

sn_edges_s=list(sn_s.edges);
ee_sorce_s=[[0 for i in range(2)] for j in range(len(sn_edges_s))];
ee_target_s=[[0 for i in range(2)] for j in range(len(sn_edges_s))];


for i in range(0,len(sn_edges_s)):
  if random.random()>0.5:
    ee_sorce_s[i]=sn_edges_s[i][0];
    ee_target_s[i]=sn_edges_s[i][1];
  else:
    ee_sorce_s[i]=sn_edges_s[i][1];
    ee_target_s[i]=sn_edges_s[i][0];    


Cee_s.connect(i=ee_sorce_s, j=ee_target_s);
Cee_s.w= temp_s[0:len(sn_edges_s)]/EPSPtoW;
Cee_s.delay='(1+2*rand())*ms';
Cee_s.p=1.-0.1/(0.1+EPSPtoW*Cee_s.w);


#Weak weights
sn_w=networkx.watts_strogatz_graph(len(Pe),2*int(len(temp_w)/len(Pe)),1.0,1);

sn_edges_w=list(sn_w.edges);
ee_sorce_w=[[0 for i in range(2)] for j in range(len(sn_edges_w))];
ee_target_w=[[0 for i in range(2)] for j in range(len(sn_edges_w))];


for i in range(0,len(sn_edges_w)):
  if random.random()>0.5:
    ee_sorce_w[i]=sn_edges_w[i][0];
    ee_target_w[i]=sn_edges_w[i][1];
  else:
    ee_sorce_w[i]=sn_edges_w[i][1];
    ee_target_w[i]=sn_edges_w[i][0];


Cee.connect(i=ee_sorce_w, j=ee_target_w);
Cee.w= temp_w[0:len(sn_edges_w)]/EPSPtoW
Cee.delay='(1+2*rand())*ms'
Cee.p=1.-0.1/(0.1+EPSPtoW*Cee.w)

print("E to E network has been constructed.");



Cei = Synapses(Pe, Pi, model="""w:1
                              p:1""",
                              on_pre='ge_post+=w')
Cei.connect(p=0.1)
#Ce.w='rand()'
Cei.w= 0.018
Cei.delay='(2*rand())*ms'
#for i,j in zip(Ce.i,Ce.j):
#	Ce.w[i,j]= numpy.random.lognormal(mu_l, sigma_l)/120

#Ce.w[:,:]= 0.

#for i,j in zip(Ce.i,Ce.j):
#	Ce.p[i,j]=1.-0.1/(0.1+120*Ce.w[i,j])
Cei.p=1




#Ce.p[:,:]=1
#Cii = Synapses(Pi, Pi, on_pre='gi_post+=0.0025')
Cii = Synapses(Pi, Pi, model="""w:1
                              p:1""",
                              on_pre='gi_post+=w')
Cii.connect(p=0.5)
Cii.w= 0.0025
#Cii.w= 0.005
#Ci = Synapses(Pi, P, on_pre='gi-=0./ms',delay=(2*rand())*ms)
Cii.delay='(2*rand())*ms'
Cii.p=1



#Cie = Synapses(Pi, Pe, on_pre='gi_post+=0.002')
Cie = Synapses(Pi, Pe, model="""w:1
                              p:1""",
                              on_pre='gi_post+=w')
Cie.connect(p=0.5)
Cie.w= 0.002
#Cie.w= 0.004

#Ci = Synapses(Pi, P, on_pre='gi-=0./ms',delay=(2*rand())*ms)
Cie.delay='(2*rand())*ms'
Cie.p=1

#S = Synapses(P_ex, P, on_pre='v_post+=21*mV')
S = Synapses(P_ex, P, model="""w:1
                              p:1""",
                              on_pre='ge_post+=w')
S.connect(j='i')
S.w=0.2

M = SpikeMonitor(P)
M_p_e = PopulationRateMonitor(Pe)


#Me_s =StateMonitor(Pe,'v',record=True)
M_p_i = PopulationRateMonitor(Pi)
#Mi_s =StateMonitor(Pi,'v',record=True)

print("start");

Lambda=30;
run(0.5*second,'stdout')
store('before ASSR1')

Lambda=0;
run(2.5*second,'stdout')

store('before ASSR2')

Lambda=10;
run(4*second,'stdout')
store('before ASSR')

Lambda=0;
run(3*second,'stdout')
store('After ASSR')

print("Finish");


plt.rcParams["font.size"] = 10
#subplot(3,1,1)
#xlim(8000,9000)
#plot(M.t/ms, M.i, 'o', markersize=0.1)
#ylabel('Neuron Index')

#subplot(3,1,2)
#xlim(8000,9000)
smooth_rate_e=M_p_e.smooth_rate(window='gaussian', width=1*ms)


smooth_rate_i=M_p_i.smooth_rate(window='gaussian', width=1*ms)

#plot(M_p_e.t/ms, smooth_rate_e,label='ex')
#plot(M_p_i.t/ms, smooth_rate_i,label='inh')


#ylabel('Spiking Rate (Hz)')

#subplot(3,1,3)
#xlim(8000,9000)
#ylim(-75,-50)
#plot(Me_s.t/ms, Me_s.v[0]/ms, label='ex')
#plot(Mi_s.t/ms, Mi_s.v[0]/ms, label='inh')

#ylabel('v [mV]')
#xlabel('time [ms]')

f=open(args[1]+'_'+args[2]+'_raster_t.dat','w')
for x in M.t/ms:
    f.write(str(x)+"\n")
f.close
f=open(args[1]+'_'+args[2]+'_raster.dat','w')
for x in M.i:
    f.write(str(x)+"\n")
f.close

#f=open(args[1]+'_rate_t.dat','w')
#for x in M_p_e.t/ms:
#    f.write(str(x)+"\n")
#f.close
f=open(args[1]+'_'+args[2]+'_e.dat','w')
for x in smooth_rate_e/Hz:
    f.write(str(x)+"\n")
f.close



f=open(args[1]+'_'+args[2]+'_i.dat','w')
for x in smooth_rate_i/Hz:
    f.write(str(x)+"\n")
f.close



#f=open(args[1]+'_v_time.dat','w')
#for x in Me_s.t/ms:
#    f.write(str(x)+"\n")
#f.close

#f=open(args[1]+'_e_m.dat','w')
#for x in Me_s.v[0]/ms:
#    f.write(str(x)+"\n")
#f.close
    
#f=open(args[1]+'_i_m.dat','w')
#for x in Mi_s.v[0]/ms:
#    f.write(str(x)+"\n")
#f.close


#fileformat=".eps"
#filename=args[1]+fileformat

#savefig(filename)
