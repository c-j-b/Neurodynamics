from neuron import h, gui
import pickle
from matplotlib import pyplot


# Channels
# pump.mod (pump) - Sodium Potassium Pump
# newleak.mod (leak) - Leak current including GABA_A receptor current
# newkca.mod (kca) - Calcium dependent potassium channel (small conductance SK)
# newhh3.mod (hh3) - squid sodium, potassium delayed rectifier, and potassium A channels
# newcachan.mod (cachan) - calcium channels (L, N, and T types) 
# nabalan.mod (nabalan) - Sodium ion accumulation without diffusion
# cabal.mod (cabalan) - Calcium ion accumulation without diffusion and buffering
# stim.mod (MyIClamp) - Applied current

# dopaminergic.hoc - 3 soma, 38 dendrites


### --- Define Model --- ###
# Topology, 1 soma, 1 dendrites
soma = h.Section(name='soma')
dend = h.Section(name = 'dend')
dend.connect(soma(1))
h.topology()	# print topology

# Geometry 
soma.L = soma.diam = 20 # Makes a soma 20 microns in diameter
dend.L = 1 # microns
dend.diam = 1 # microns

# Add channels to needed channels to soma/dendrites
for sec in h.allsec():
    sec.Ra = 100    # Axial resistance in Ohm * cm
    sec.cm = 1      # Membrane capacitance in micro Farads / cm^2

# Insert active Hodgkin-Huxley current in the soma
# PARAMETERS
na_cond =  550.0e-6 
kdr_cond = 665.0e-6
ca_cond = 11.196e-6
kca_cond = 59.0e-6
a_cond_s = 570.0e-6
a_cond_p = 285.0e-6
a_cond_d = 266.0e-6
#stronger gA *1.28 =  729.6, 364.8, 340.48
iapl = 0 #in nA, -0.180nA=-180pA

soma.insert('hh3')
for seg in soma:
    seg.hh3.gnabar = na_cond
    seg.hh3.gkhhbar = kdr_cond  
    seg.hh3.gkabar = a_cond_s    
    seg.hh3.qs = 56.0     
    seg.hh3.qv = 8.0

# Insert passive current in the dendrite
dend.insert('pas')
for seg in dend:
    seg.pas.g = 0.001  # Passive conductance in S/cm2
    seg.pas.e = -65    # Leak reversal potential mV

# Run simulation
h.psection()
v_vec = h.Vector()             # Membrane potential vector
t_vec = h.Vector()             # Time stamp vector
v_vec.record(soma(0.5)._ref_v)
t_vec.record(h._ref_t)
h.run(40.0)
pyplot.figure(figsize=(8,4)) # Default figsize is (8,6)
pyplot.plot(t_vec, v_vec)
pyplot.xlabel('time (ms)')
pyplot.ylabel('mV')
pyplot.show()

# Pickle
# with open('t_vec.p', 'wb') as t_vec_file:
#     pickle.dump(t_vec, t_vec_file)
# with open('v_vec.p', 'wb') as v_vec_file:
#     pickle.dump(v_vec, v_vec_file)


