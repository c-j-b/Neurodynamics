from brian2 import *
from scipy import stats

defaultclock.dt = 0.01*ms

# Define spatial neuron geometry
morpho = Soma(diameter=20*um)   # spherical soma
morpho.dendrites = Cylinder(length=1*um, diameter=1*um, n=10)   # 10 identical cylindrical dendrites

morpho.topology()  # print topology?

# Nernst Potentials
E_Ca = 100*mV
E_K = -90*mV
E_L = -50*mV
E_Na = 55*mV
E_Ks = -90*mV  # check this value
E_CaL = 60*mV  # check this value

# Conductances
g_bar_CaL = 0.2*msiemens/cm**2
g_bar_K = 0.4*msiemens/cm**2
g_bar_KCa = 0.3*msiemens/cm**2
g_bar_Na = 150*msiemens/cm**2
g_bar_Ks = 4*msiemens/cm**2
gc = 2.5*msiemens/cm**2
gL = 0.05*msiemens/cm**2
g_bar_NMDA = 0.4*msiemens/cm**2
g_bar_Ca = 0.20

# Half activation
VHK = 20*mV
VSK = 13.8*mV

# Misc. constants
KCa = 250*nmolar
beta = 0.05
z_Ca = 2
P_Ca = 0.25*nmole/second
C = 1*ufarad/cm**2
F = 0.096485
Mgcon = 1.4*nmolar
Cacon_S = 100*nmolar  # check this unit
Cacon_D = 100*nmolar  # check this unit
r_s= 10

alpha = 1/(z_Ca*F)

#
# Equations
eqs_asu = '''
# Transmembrane current for soma, vs = soma voltage
I_Ks_S = g_Ks_S * (E_Ks - Vs) : amp
I_Na_S = g_Na_S * (E_Na - Vs) : amp
I_CaL_S = g_CaL_S * (E_CaL - Vs) : amp
I_K_S = g_K_S * (E_K - Vs) : amp
I_KCa_S = g_KCa_S*(E_K - Vs) : amp #######redo...why?
I_L_S = g_L*(E_L - Vs) : amp
I_couple_S = n_d*g_c*(r_d**2)*r_s/(l_s*(l_s*r_s**2+l_d*r_d**2))*(Vd-Vs) : amp

I_Ks_D = g_Ks_D * (E_Ks - Vd) : amp
I_Na_D = g_Na_D * (E_Na - Vd) : amp
I_CaL_D = g_CaL_D * (E_CaL - Vd) : amp
I_K_D = g_K_D * (E_K - Vd) : amp
I_KCa_D = g_KCa_D*(E_K - Vd) : amp  ####### redo...?
I_L_D = g_L*(E_L - Vd) : amp
I_couple_D = g_c*((r_s**2)*r_d/(l_d*(l_d*r_d**2 + l_s*r_s**2))*(Vs-Vd)) : amp
I_NMDA_D = g_NMDA * (E_NMDA - Vd) : amp

dV_s/dt = (I_app + I_Ks_S + I_Na_S + I_CaL_S + I_K_S + I_KCa_S +I_L_S + I_couple_S)/C : amp/meter**2
dV_d/dt = (I_Ks_D + I_Na_D + I_CaL_D + I_K_D + I_KCa_D +I_L_D + I_NMDA_D + I_couple_D)/C : amp/meter**2

I : amp (point current) # applied current

dh/dt = alpha_h * (1-h) - beta_h * h : 1
dn/dt = alpha_n * (1-n) - beta_n * n : 1
# Activation/Deactivation gating variables
alpha_m = -0.32*(v +33)/(exp(-(v+33)/4.5)-1) : Hz
beta_m = 0.28*(v+4)/(exp((v+4)/10.4)-1) : Hz
alpha_h = 0.0196*exp(-(v+47)/18) : Hz
beta_h = 2.45/(1+exp(-(v+24)/4)) : Hz
alpha_n = -0.3584*((v-2)/exp(-(v+2)/3)-1) : Hz
beta_n = 0.56*exp(-(v+20)/5.8) : Hz
alpha_c  = -0.0032*(v+50)/(exp(-(v+50)/5) - 1) : Hz
beta_c  = 0.05*exp(-(v+55)/40) : Hz
# Conductances
g_KCa_S = g_bar_KCa * (Cacon_S*Cacon_S*Cacon_S*Cacon_S)/(Cacon_S*Cacon_S*Cacon_S*Cacon_S+KCa*KCa*KCa*KCa)
g_K_S = g_bar_K*1/(1+exp(-(Vs-V_HK)/V_SK))
g_CaL_S = g_bar_CaL * (alpha_c_S/(alpha_c_S+beta_c_S)) * (alpha_c_S/(alpha_c_S+beta_c_S))*(alpha_c_S/(alpha_c_S+beta_c_S))*(alpha_c_S/(alpha_c_S+beta_c_S))
g_Na_S = g_bar_Na*m_inf_S*m_inf_S*m_inf_S*h_S
g_Ks_S = g_bar_Ks * n_S*n_S*n_S*n_S
g_KCa_D = g_bar_KCa * (Cacon_D*Cacon_D*Cacon_D*Cacon_D)/(Cacon_D*Cacon_D*Cacon_D*Cacon_D+KCa*KCa*KCa*KCa)
g_K_D = g_bar_K*1/(1+exp(-(Vd-V_HK)/V_SK))
g_CaL_D = g_bar_CaL * (alpha_c_D/(alpha_c_D+beta_c_D)) * (alpha_c_D/(alpha_c_D+beta_c_D))*(alpha_c_D/(alpha_c_D+beta_c_D))*(alpha_c_D/(alpha_c_D+beta_c_D))
g_Na_D = g_bar_Na*m_inf_D*m_inf_D*m_inf_D*h_D
g_Ks_D = g_bar_Ks * n_D*n_D*n_D*n_D
g_NMDA = g_bar_NMDA*(1/(1+(Mgcon/10)*exp(Vd/12.5)))
# Misc.
m_inf = alpha_m/(alpha_m+beta_m)
dCacon_S/dt = 2*beta*(alpha*I_CaL_S - P_Ca*Cacon_S)/r_s
'''



# From Brian EX - Typical equations
eqs = '''
# The same equations for the whole neuron, but possibly different parameter values
# distributed transmembrane current
Im = gl * (El-v) + gNa * m**3 * h * (ENa-v) + gK * n**4 * (EK-v) : amp/meter**2
I : amp (point current) # applied current
dm/dt = alpham * (1-m) - betam * m : 1
dn/dt = alphan * (1-n) - betan * n : 1
dh/dt = alphah * (1-h) - betah * h : 1
alpham = (0.1/mV) * (-v+25*mV) / (exp((-v+25*mV) / (10*mV)) - 1)/ms : Hz
betam = 4 * exp(-v/(18*mV))/ms : Hz
alphah = 0.07 * exp(-v/(20*mV))/ms : Hz
betah = 1/(exp((-v+30*mV) / (10*mV)) + 1)/ms : Hz
alphan = (0.01/mV) * (-v+10*mV) / (exp((-v+10*mV) / (10*mV)) - 1)/ms : Hz
betan = 0.125*exp(-v/(80*mV))/ms : Hz
gNa : siemens/meter**2
'''

neuron = SpatialNeuron(morphology=morpho, model=eqs_asu, method="exponential_euler",
                       refractory="m > 0.4", threshold="m > 0.5",
                       Cm=1*uF/cm**2, Ri=35.4*ohm*cm)
neuron.v = 0*mV
neuron.h = 1
neuron.m = 0
neuron.n = .5
neuron.I = 0*amp
neuron.gNa = gNa
M = StateMonitor(neuron, 'v', record=True)
spikes = SpikeMonitor(neuron)

run(50*ms, report='text')
neuron.I[0] = 1*uA  # current injection at one end
run(3*ms)
neuron.I = 0*amp
run(50*ms, report='text')

# Calculation of velocity
slope, intercept, r_value, p_value, std_err = stats.linregress(spikes.t/second, neuron.distance[spikes.i]/meter)
print("Velocity = %.2f m/s" % slope)

subplot(211)
for i in range(10):
    plot(M.t/ms, M.v.T[:, i*100]/mV)
ylabel('v')
subplot(212)
plot(spikes.t/ms, spikes.i*neuron.length[0]/cm, '.k')
plot(spikes.t/ms, (intercept+slope*(spikes.t/second))/cm, 'r')
xlabel('Time (ms)')
ylabel('Position (cm)')
show()
