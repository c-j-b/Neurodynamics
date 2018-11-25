from neuron import h, gui
import pickle
from matplotlib import pyplot

soma = h.Section(name='soma')
h.psection()
soma.insert('pas')
print("type(soma) = {}".format(type(soma)))
print("type(soma(0.5)) ={}".format(type(soma(0.5))))
asyn = h.AlphaSynapse(soma(0.5))
asyn.onset = 20
asyn.gmax = 1
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
with open('t_vec.p', 'wb') as t_vec_file:
    pickle.dump(t_vec, t_vec_file)
with open('v_vec.p', 'wb') as v_vec_file:
    pickle.dump(v_vec, v_vec_file)

# # Unpickle
# with open('t_vec.p', 'rb') as t_vec_file:
#     t_vec = pickle.load(t_vec_file)
# with open('v_vec.p', 'rb') as vec_file:
#     v_vec = pickle.load(vec_file)

# # Confirm
# pyplot.figure(figsize=(8, 4))  # Default figsize is (8,6)
# pyplot.plot(t_vec, v_vec)
# pyplot.xlabel('time (ms)')
# pyplot.ylabel('mV')
# pyplot.show()