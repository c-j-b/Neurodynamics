C = 1*microF/cm**2
E_Ca = 100
E_L= -50
E_K= -90
E_Na=55
E_NMDA=0
E_Ks = ?
E_CaL = ?

g_bar_CaL = 0.20
g_bar_K = 0.4
g_bar_KCa= 0.3
g_bar_Na=150
g_bar_Ks=4
g_bar_NMDA = 0.4

V_HK = 20
V_SK= 13.8
KCa= 250
Mgcon = 1.4
Cacon_S = ?
Cacon_D = ?
g_L= 0.05
g_c=2.5
beta= 0.05
r_s= 10
r_d= 0.5
z_Ca=2
P_Ca=0.25
n_d=10
l_s=1
l_d=1
F=0.096485
alpha = 1/(z_Ca*F)
n = ? # gating variable for activation of the delayed rectifier potassium conductance

g_bar_Ca = 0.20
g_Ca = g_bar_Ca*
I_Ca_S = I_Ca_L

d(Cacon_S)/dt = 2*beta*(alpha*I_Ca_S - P_Ca*Cacon_S)/r_s

alpha_m_S = -0.32*(Vs +33)/(exp(-(Vs+33)/4.5)-1)
beta_m_S = 0.28*(Vs+4)/(exp((Vs+4)/10.4)-1)
alpha_h_S = 0.0196*exp(-(Vs+47)/18)
beta_h_S = 2.45/(1+exp(-(Vs+24)/4))
alpha_n_S = -0.3584*((Vs-2)/exp(-(Vs+2)/3)-1)
beta_n_S = 0.56*exp(-(Vs+20)/5.8)
alpha_c_S = ?
beta_c_S = ?

alpha_m_D = -0.32*(Vd +33)/(exp(-(Vd+33)/4.5)-1)
beta_m_D = 0.28*(Vd+4)/(exp((Vd+4)/10.4)-1)
alpha_h_D = 0.0196*exp(-(Vd+47)/18)
beta_h_D = 2.45/(1+exp(-(Vd+24)/4))
alpha_n_D = -0.3584*((Vd-2)/exp(-(Vd+2)/3)-1)
beta_n_D = 0.56*exp(-(Vd+20)/5.8)
alpha_c_D = ?
beta_c_D = ?

m_inf_S = alpha_m_S/(alpha_m_S+beta_m_S)   #eqn 2.5
m_inf_D = alpha_m_D/(alpha_m_D+beta_m_D)

d(h_S)/dt =alpha_h_S*(1-h_S)-beta_h_S*h_S
d(n_S)/dt =alpha_n_S*(1-n_S)-beta_n_S*n_S
d(h_D)/dt =alpha_h_S*(1-h_D)-beta_h_D*h_D
d(n_D)/dt =alpha_n_S*(1-n_D)-beta_n_D*n_D

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

I_app =
I_Ks_S = g_Ks_S * (E_Ks - Vs)
I_Na_S = g_Na_S * (E_Na - Vs)
I_CaL_S = g_CaL_S * (E_CaL - Vs)
I_K_S = g_K_S * (E_K - Vs)
I_KCa_S = g_KCa_S*(E_K - Vs) #######redo
I_L_S = g_L*(E_L - Vs)
I_couple_S = n_d*g_c*(r_d**2)*r_s/(l_s*(l_s*r_s**2+l_d*r_d**2))*(Vd-Vs)

I_Ks_D = g_Ks_D * (E_Ks - Vd)
I_Na_D = g_Na_D * (E_Na - Vd)
I_CaL_D = g_CaL_D * (E_CaL - Vd)
I_K_D = g_K_D * (E_K - Vd)
I_KCa_D = g_KCa_D*(E_K - Vd)  ####### redo
I_L_D = g_L*(E_L - Vd)
I_couple_D = g_c*((r_s**2)*r_d/(l_d*(l_d*r_d**2 + l_s*r_s**2))*(Vs-Vd))
I_NMDA = g_NMDA * (E_NMDA - Vd)

d(Vs)/dt = (I_app + I_Ks_S + I_Na_S + I_CaL_S + I_K_S + I_KCa_S +I_L_S + I_couple_S)/C
d(Vd)/dt = (I_Ks_D + I_Na_D + I_Ca_D + I_K_D + I_KCa_D +I_L_D + I_NMDA_D + I_couple_D)/C
