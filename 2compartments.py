C = 1*microF/cm**2
E_Ca = 100
E_L= -50
E_K= -90
E_Na=55
E_NMDA=0
E_Ks = -90 #mV
E_CaL = 60 #mV

g_bar_CaL = 0.20 #OK
g_bar_K = 0.4 #OK
g_bar_KCa= 0.3 #OK
g_bar_Na=150 #OK
g_bar_Ks=4 #OK
g_bar_NMDA = 0.4 #OK

V_HK = 20 #OK
V_SK= 13.8 #OK
KCa= 250 #OK
Mgcon = 1.4
Cacon_S = 100 #nM
Cacon_D = 100 #nM
g_L= 0.05 #OK
g_c=2.5 #OK
beta= 0.05 #OK
r_s= 10 #OK
r_d= 0.5#OK
z_Ca=2 #OK
P_Ca=0.25 #OK
n_d=10 #OK
l_s=1 #OK
l_d=1 #OK
F=0.096485 #OK
alpha = 1/(z_Ca*F) #OK
n = ? # gating variable for activation of the delayed rectifier potassium conductance

g_bar_Ca = 0.20
g_Ca = g_bar_Ca
I_Ca_S = I_Ca_L

d(Cacon_S)/dt = 2*beta*(alpha*I_Ca_S - P_Ca*Cacon_S)/r_s

alpha_m = -0.32*(v +33)/(exp(-(v+33)/4.5)-1)
beta_m = 0.28*(v+4)/(exp((v+4)/10.4)-1)
alpha_h = 0.0196*exp(-(v+47)/18)
beta_h = 2.45/(1+exp(-(v+24)/4))
alpha_n = -0.3584*((v-2)/exp(-(v+2)/3)-1)
beta_n = 0.56*exp(-(v+20)/5.8)
alpha_c  = -0.0032*(v+50)/(exp(-(v+50)/5) - 1) # from Kuznetsov 2006
beta_c  = 0.05*exp(-(v+55)/40) # from Kuznetsov 2006

m_inf_S = alpha_m/(alpha_m+beta_m)   #eqn 2.5
m_inf_D = alpha_m/(alpha_m+beta_m)

d(h_S)/dt =alpha_h*(1-h_S)-beta_h*h_S
d(n_S)/dt =alpha_n*(1-n_S)-beta_n*n_S
d(h_D)/dt =alpha_h*(1-h_D)-beta_h*h_D
d(n_D)/dt =alpha_n*(1-n_D)-beta_n*n_D

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

I_Ks_S = g_Ks_S * (E_Ks - Vs)
I_Na_S = g_Na_S * (E_Na - Vs)
I_CaL_S = g_CaL_S * (E_CaL - Vs)
I_K_S = g_K_S * (E_K - Vs)
I_KCa_S = g_KCa_S*(E_K - Vs) #######redo...why?
I_L_S = g_L*(E_L - Vs)
I_couple_S = n_d*g_c*(r_d**2)*r_s/(l_s*(l_s*r_s**2+l_d*r_d**2))*(Vd-Vs)

I_Ks_D = g_Ks_D * (E_Ks - Vd)
I_Na_D = g_Na_D * (E_Na - Vd)
I_CaL_D = g_CaL_D * (E_CaL - Vd)
I_K_D = g_K_D * (E_K - Vd)
I_KCa_D = g_KCa_D*(E_K - Vd)  ####### redo...?
I_L_D = g_L*(E_L - Vd)
I_couple_D = g_c*((r_s**2)*r_d/(l_d*(l_d*r_d**2 + l_s*r_s**2))*(Vs-Vd))
I_NMDA_D = g_NMDA * (E_NMDA - Vd)

d(Vs)/dt = (I_app + I_Ks_S + I_Na_S + I_CaL_S + I_K_S + I_KCa_S +I_L_S + I_couple_S)/C
d(Vd)/dt = (I_Ks_D + I_Na_D + I_CaL_D + I_K_D + I_KCa_D +I_L_D + I_NMDA_D + I_couple_D)/C
