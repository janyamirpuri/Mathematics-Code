# -*- coding: utf-8 -*-
"""
Created on Sun Apr  2 14:42:56 2023

@author: janya
"""


import math
import numpy as np
import matplotlib.pyplot as plt

#plt.rcParams.update(plt.rcParamsDefault)
#plt.rcParams['text.usetex'] = True

#PART II
#define variables
g = 9.81
l = 5
#define initials
T0 = math.pi/4 #initial theta
dT0 = 0 #initial dtheta/dt

#define timescale

tau = 2*math.pi*math.pow(l/g, 0.5) #timescale in periods

tf = 5*tau #how much time to run in total

dt_choices = [0.005*tau, 0.01*tau, 0.0001*tau] #time step for eulers

plt.figure()
#plt.title("Euler's Method Results by Varying Time Step")
plt.xlabel("t / $\\psi}$")
plt.ylabel("$\\theta$")
lines = []

for dt in dt_choices:
    num_steps = math.floor(tf/dt)+1 #total number of steps to make
    t = np.linspace(0, tf, num_steps) #array of each time we evaluate at
    dt = t[1]-t[0] #redefine dt to match the length of time step in t
    
    T = np.zeros(num_steps)
    T[0] = T0
    
    dT = np.zeros(num_steps)
    dT[0] = dT0
    
    #euler's method
    for i in range(num_steps-1):
        T[i+1] = T[i] + dT[i]*dt
        dT[i+1] = dT[i] - (g/l)*(math.sin(T[i]))*dt
        
    lines += plt.plot(t/tau,T, label=str(round(dt/tau, 4)))

Y2 = T0*np.cos(math.pow(g/l, 0.5)*t) #analytical solution
lines += plt.plot(t/tau, Y2, label = "Analytic", linestyle = "dashdot")

 
labels = [l.get_label() for l in lines]
plt.legend(lines, labels, ncol = 2, frameon=False, title="$\\Delta$ t / $\\psi$")
plt.show()



#%%

#PART II.2
#variables
g = 9.81
l = 5
#initial values
dT0 = 0 #dtheta/dt
#define timescale

tau = 2*math.pi*math.pow(l/g, 0.5) #time scale

tf = 5*tau #total time to run

dt = 0.0001*tau #timestep

T0_choices =  [math.pi/6, math.pi/4, math.pi/3, 2*math.pi/5] #varying initial theta

num_steps = math.floor(tf/dt)+1 #total number of steps to make
t = np.linspace(0, tf, num_steps) #array of each time we evaluate at
dt = t[1]-t[0] #redefine dt to match the length of time step in t


plt.figure()
#plt.title("Euler's Method Results by Varying Initial Angle")
plt.xlabel("t / $\\psi$")
plt.ylabel("$\\theta$")
lines = []

plt.ylim(-math.pi/2, 0.80*math.pi)

for T0 in T0_choices:
    
    T = np.zeros(num_steps)
    T[0] = T0
    
    dT = np.zeros(num_steps)
    dT[0] = dT0
    
    #euler's method
    for i in range(num_steps-1):
        T[i+1] = T[i] + dT[i]*dt
        dT[i+1] = dT[i] - (g/l)*(math.sin(T[i]))*dt
        
    lines += plt.plot(t/tau,T, label=(str(round(T0/math.pi,4)) + "$\\pi$"))
    
labels = [l.get_label() for l in lines]
plt.legend(lines, labels,ncol = 2, frameon=False, title="$\\theta_0$", loc = "upper right")
plt.show()


figure, axis = plt.subplots(4, 1, figsize=(10,20))
#comparing euler to analytical per T0
j = 0
for T0 in T0_choices:
    
    T = np.zeros(num_steps)
    T[0] = T0
    
    dT = np.zeros(num_steps)
    dT[0] = dT0
    
    #euler's method
    for i in range(num_steps-1):
        T[i+1] = T[i] + dT[i]*dt
        dT[i+1] = dT[i] - (g/l)*(math.sin(T[i]))*dt
        
    Y2 = T0*np.cos(math.pow(g/l, 0.5)*t)
    axis[j].plot(t/tau,T, label="Euler's")
    axis[j].plot(t/tau, Y2, label = "Analytic")
    axis[j].set_ylabel("$\\theta$", fontdict = {'fontsize':18})
    axis[j].set_title(str(chr(j+97)) +") $\\theta_0 = $ " + str(round(T0/math.pi,4)) + str("$\\pi$"), fontdict = {"fontsize":18})
    axis[j].legend(loc="lower right")
    j+=1
    
plt.xlabel("t / $\\psi$", fontdict = {'fontsize': 18})
plt.show()

#%%                  
#PART III

#defining variables
g = 9.81
l = 5
a = 1
m = 5
#initial values
T0 = math.pi/4
dT0 = 0

tau = 2*math.pi*math.pow(l/g, 0.5) #timescale

tf = 5*tau #how long to run

dt = 0.0001*tau #timestep

k_choices = [0, 0.1, 0.5, 1] #drag constants to vary over

num_steps = math.floor(tf/dt)+1 #total number of steps to make
t = np.linspace(0, tf, num_steps) #array of each time we evaluate at
dt = t[1]-t[0] #redefine dt to match the length of time step in t

plt.figure()
#plt.title("Euler's Method Results by Varying Drag Constant")
plt.xlabel("t / $\\psi$")
plt.ylabel("$\\theta$")
lines = []

plt.ylim(-math.pi/3,0.5*math.pi)

for k in k_choices:
    
    T = np.zeros(num_steps)
    T[0] = T0
    
    dT = np.zeros(num_steps)
    dT[0] = dT0
    
    #euler's method
    for i in range(num_steps-1):
        T[i+1] = T[i] + dT[i]*dt
        dT[i+1] = dT[i] - (k/m)*(math.pi*a*a)*(dT[i])*dt - (g/l)*(math.sin(T[i]))*dt
        
    lines += plt.plot(t/tau,T, label=(str(k)))
    
labels = [l.get_label() for l in lines]
plt.legend(lines, labels,ncol = 2, frameon=False, title="Drag Constant k")
plt.show()

#dtheta/dt vs theta for varying drag constant
figure, axis = plt.subplots(4, 1, figsize=(10, 20))
j = 0
for k in k_choices:
    
    T = np.zeros(num_steps)
    T[0] = T0
    
    dT = np.zeros(num_steps)
    dT[0] = dT0
    
    #euler's method
    for i in range(num_steps-1):
        T[i+1] = T[i] + dT[i]*dt
        dT[i+1] = dT[i] - (k/m)*(math.pi*a*a)*(dT[i])*dt - (g/l)*(math.sin(T[i]))*dt
    
    axis[j].set_ylabel("$d\\theta/dt$", fontdict = {"fontsize":18})
    axis[j].plot(T,dT)
    axis[j].set_title(str(chr(j+97)) +") Drag Constant: " + str(k),  fontdict = {"fontsize":18})

    j+=1
    
plt.xlabel("$\\theta$", fontdict = {"fontsize":18})

plt.show()






figure, axis = plt.subplots(4, 1, figsize=(10,20))

#variables/initial values

g = 9.81
l = 5
T0 = math.pi/4
a = 1
m = 5
dT0 = 0
k = 0.5

tau = 2*math.pi*math.pow(l/g, 0.5) #timescale

tf = 5*tau #length of time to run

dt = 0.0001*tau #timestep

# defining the roots of the characteristic equation with mu and theta
LAMBDA = 0.5*(-k/m)*(math.pi*a*a)

MU = math.sqrt(4*(g/l)-(k/m)*(math.pi*a*a))/2

T01= 0

#analytical solution
def func(t):
    y = T01*math.exp(LAMBDA*t)*math.cos(MU*t) -(LAMBDA*T01/MU)*math.exp(LAMBDA*t)*math.sin(MU*t)
    return y

#eulers vs analytic for varying initial theta

j = 0

for T0 in T0_choices:
    T01 = T0
    
    T = np.zeros(num_steps)
    T[0] = T0
    
    dT = np.zeros(num_steps)
    dT[0] = dT0
    
    #euler's method
    for i in range(num_steps-1):
        T[i+1] = T[i] + dT[i]*dt
        dT[i+1] = dT[i] - (k/m)*(math.pi*a*a)*(dT[i])*dt - (g/l)*(math.sin(T[i]))*dt
        
    Y2 = [func(i) for i in t]
    axis[j].plot(t/tau,T, label="Euler's")
    axis[j].plot(t/tau, Y2, label = "Analytic")
    axis[j].set_ylabel("$\\theta$", fontdict = {'fontsize':18})
    axis[j].set_title(str(chr(j+97)) +") $\\theta_0 = $ " + str(round(T0/math.pi,4)) + str("$\\pi$"), fontdict = {"fontsize":18})
    axis[j].legend(loc="lower right")
    j+=1
    
plt.xlabel("t / $\\psi$", fontdict = {'fontsize': 18})
plt.show()


#%%

#PART IV

#variables/initials
g = 9.81
l = 5
T0 = math.pi/4
a = 1
m = 5
dT0 = 0

tau = 2*math.pi*math.pow(l/g, 0.5) #define timescale

tf = 5*tau #total length of time

dt = 0.000001*tau #timestep

k_choices = [0, 0.1, 0.5, 1] #varying drag constant

num_steps = math.floor(tf/dt)+1 #total number of steps to make
t = np.linspace(0, tf, num_steps) #array of each time we evaluate at
dt = t[1]-t[0] #redefine dt to match the length of time step in t

plt.figure()
#plt.title("Euler's Method Results by Varying Drag Constant")
plt.xlabel("t / $\\psi$")
plt.ylabel("$\\theta$")
lines = []

plt.ylim(-math.pi/3,0.5*math.pi)

for k in k_choices:
    
    T = np.zeros(num_steps)
    T[0] = T0
    
    dT = np.zeros(num_steps)
    dT[0] = dT0
    
    #euler's method
    for i in range(num_steps-1):
        T[i+1] = T[i] + dT[i]*dt
        dT[i+1] = dT[i] - (k/m)*(2*a*l)*(dT[i])*dt - (3*g/2*l)*(math.sin(T[i]))*dt
        
    lines += plt.plot(t/tau,T, label=(str(k)))
    
labels = [l.get_label() for l in lines]
plt.legend(lines, labels,ncol = 2, frameon=False, title="Drag Constant k")
plt.show()
                         
#%%
#PART V
#variables and initials
#m1, l1 corresponds to top mass, m2,l2 corresponds to bottom mass
g = 9.81
l1 = 5
m1 = 5
l2 = 5
m2 = 5
T10 = math.pi/4 #initial theta of first mass

dT0 = 0

tau = 2*math.pi*math.pow(l1/g, 0.5) #timescale

tf = 5*tau #total length of time to run

dt = 0.000001*tau #timestep

num_steps = math.floor(tf/dt)+1 #total number of steps to make
t = np.linspace(0, tf, num_steps) #array of each time we evaluate at
dt = t[1]-t[0] #redefine dt to match the length of time step in t

#define timescale

T20_choices = [math.pi/16, math.pi/8, math.pi/4, math.pi/2] #initial theta of second mass

plt.figure()
#plt.title("Euler's Method Results by Varying Second Initial Theta")
plt.xlabel("t / $\\psi$")
plt.ylabel("$\\theta_1$")
lines = []

for T20 in T20_choices:
    
    T1 = np.zeros(num_steps)
    T1[0] = T10
    T2 = np.zeros(num_steps)
    T2[0] = T20
    
    dT1 = np.zeros(num_steps)
    dT1[0] = dT0
    
    dT2 = np.zeros(num_steps)
    dT2[0] = dT0
    
    #euler's method
    for i in range(num_steps-1):
        T2[i+1] = T2[i] + dT2[i]*dt
        dT2[i+1] = dT2[i] - (g/l2)*(math.sin(T2[i]))*dt
        
        T1[i+1] = T1[i] + dT1[i]*dt
        dT1[i+1] = dT1[i] - ((g/l1)*(math.sin(T1[i])))*dt - ((m2*g/(l1*m1))*(math.sin(T1[i] - T2[i])))*dt
        
        
    lines += plt.plot(t/tau,T1, label=(str(round(T20/math.pi,4))+"$\\pi$"))
    
labels = [l.get_label() for l in lines]
plt.legend(lines, labels,ncol = 2, frameon=False, title="$\\theta_{2_0}$", loc = "lower left")
plt.show()
    
