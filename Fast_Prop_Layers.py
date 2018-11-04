from brian2 import *
from scipy import signal
import numpy as np


def back_current(mean, std, enhancement, runtime, dt, tau, N, order):
    sampling_rate = 1 / dt
    nyquist = sampling_rate / 2
    cutoff = 1 / tau
    normalized_cutoff = cutoff / nyquist

    back_noise = np.random.normal(enhancement * mean, std, (int(runtime / dt), N))
    b, a = signal.butter(order, normalized_cutoff, btype='low', analog=False)
    back_noise_filt = signal.lfilter(b, a, back_noise[:, 0])

    if N > 1:
        for j in range(1, N):
            back_noise_filt = np.vstack((back_noise_filt, signal.lfilter(b, a, back_noise[:, j])))

    return back_noise_filt.T


def input_stimulus(mean, std, runtime, dt, tau, order):
    sampling_rate = 1 / dt
    nyquist = sampling_rate / 2
    cutoff = 1 / tau
    normalized_cutoff = cutoff / nyquist

    input_noise = (np.random.normal(mean, std, int(runtime / dt)))

    b, a = signal.butter(order, normalized_cutoff, btype='low', analog=False)

    input_noise_filt = signal.lfilter(b, a, input_noise)

    for i in range(0, len(input_noise_filt)):
        if input_noise_filt[i] < 0:
            input_noise_filt[i] = 0

    return input_noise_filt


def psth(spike_t, tbin):
    t = 0
    count = 0
    freq = []
    tvec = []

    for k in range(0, size(spike_t)):
        while spike_t[k] > t:
            freq.append(count / (20 * tbin))
            count = 0
            tvec.append(t)
            t += tbin
            if spike_t[k] <= t:
                count += 1
        count += 1

    return freq, tvec

# ----- SIMULATION ----- #


start_scope()

#zero_background = np.zeros(10000)
#print(shape(back_noise_filt))
back_current = back_current(55e-12, 70e-12, 1, 1, .0001, .002, 20, 1)
input_stimulus = input_stimulus(0, 2000e-12, 1, .0001, .5, 1)
I_back = TimedArray(back_current, dt = 0.1*ms)
I_input = TimedArray(input_stimulus, dt = 0.1*ms)
Vreversal = 0
v0 = -60e-3
R = 100e6
tau = 5*ms
taug = 5*ms

eqs = '''
dv/dt = (-v + I_back(t, i)*R + v0 + I_input(t)*R)/tau : 1

'''
#+ g*(v-Vreversal)
#dg/dt = -g/taug : 1
G = NeuronGroup(20, eqs, threshold='v>-50e-3', reset='v = -60e-3', method='exact', dt = 0.1*ms)
G.v = -60e-3

M = StateMonitor(G, 'v', record=True)
Mspk = SpikeMonitor(G)

run(1000*ms)

freq, tvec = psth(Mspk.t, 5*ms)

# PLOT NEURON VOLTAGES AND RASTER PLOT

figure(figsize=(8, 10))
subplot(8, 1, 1)
title('Layer 1 Neuron Voltages')
plot(M.t/ms, M.v[0], label='Neuron 0')
subplot(8, 1, 2)
plot(M.t/ms, M.v[1], label='Neuron 1')
subplot(8, 1, 3)
plot(M.t/ms, M.v[2], label='Neuron 1')
subplot(8, 1, 4)
plot(M.t/ms, M.v[3], label='Neuron 1')
subplot(8, 1, 5)
plot(M.t/ms, M.v[4], label='Neuron 1')
subplot(8, 1, 6)
plot(M.t/ms, M.v[5], label='Neuron 1')
subplot(8, 1, 7)
plot(M.t/ms, M.v[6], label='Neuron 1')
subplot(8, 1, 8)
plot(M.t/ms, M.v[7], label='Neuron 1')
xlabel('Time (ms)')
ylabel('v')

# RASTER PLOT

figure(2)
title('Layer 1 Raster Plot')
scatter(Mspk.t/ms, Mspk.i)
ylim(-1, 20)


# VOLTAGE PLOT
figure(3)
title('Layer 1 Neuron Voltages Combined plot')
plot(M.t/ms, M.v.T)

figure(4)
plot(back_current)

figure(5)
plot(input_stimulus)

# PLOT PSTH
figure(6).subplots_adjust(hspace=.5)
subplot(2, 1, 1)
plot(divide(tvec, ms), freq)
title('Layer 1 PSTH')
xlim(-100, 1100)
subplot(2, 1, 2)
plot(np.multiply(input_stimulus, 1e12))
title('Input Stimulus')
xlim(-1000, 11000)



# show all plots
show()
