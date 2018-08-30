import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

#System of differential equations:
# dy/dt = v
# dv/dt = -g

#Analytic solutions:
# v(t) = -gt + v_0 # speed in respect to time
# y(t) = -1/2*gt^2 + H # height in respect to time

#The initial conditions:
H = 45 # Dropping object from H=45m
v = 0 # initial speed is 0 (we just drop it, we dont use any upside or downside force)
y = H # height in respect to time (at t=0 we have y=H)

#Constants:
g = 9.81 # value of acceleration
T = 10 # Maximum value of time
N = 1000 # Number of intervals
dt = T/N # Length of one interval
t = 0 # starting at time=0


def analyticSolutions(N, t, dt):
    heightArray = []
    timeArray = []
    velocityArray = []
    y = H
    while t < T and y >= 0:
        v = -g*t
        y = -(1/2)*g*t*t+H
        heightArray.append(y)
        timeArray.append(t)
        velocityArray.append(-v) # we append (-v) because -g apperas in equations, since the acceleration vector is pointing towards the earth
        t = t + dt

    # print("\nAn object fell onto the ground.\n\n")
    # print("time t = " + str(round(t,4)) + "[s]  velocity v = " + str(round(v,4)) + "[m/s]  height y = " + str(round(y,4)) + "[m]")
    return timeArray, heightArray, velocityArray


def discreteSolutions(N, t, dt, v, y):
    heightArray = []
    timeArray = []
    velocityArray = []
    y = H
    while t < T and y >= 0:
        v = -g*dt + v # from the firt equation, Euler method
        y = v*dt + y # from the second eqation
        heightArray.append(y)
        timeArray.append(t)
        velocityArray.append(-v)
        t = t + dt

    # print("\nAn object fell onto the ground.\n\n")
    # print("time t = " + str(round(t,4)) + "[s]  velocity v = " + str(round(v,4)) + "[m/s]  height y = " + str(round(y,4)) + "[m]")
    return timeArray, heightArray, velocityArray


#SciPy Solver that solves for v(t) numerically.
y0 = 0
time = np.linspace(0, 3)

def model(y, t):
    k = 0
    p = 10
    dydt = -k * y + p
    return dydt


dataArray = discreteSolutions(N,t,dt, v, y)
dataArrayAnalytic = analyticSolutions(N, t, dt)
result = odeint(model,y0,time)


fig = plt.figure(1, figsize=(12,6))

ax1 = fig.add_subplot(211)
line_1, = ax1.plot(dataArray[0], dataArray[1], 'r--', linewidth=3, label='heightNumerical')
line_2, = ax1.plot(dataArrayAnalytic[0], dataArrayAnalytic[1], 'c', label='heightAnalytic')
ax1.axis([0, 3.1, 0, H+3])
ax1.legend(handles=[line_1, line_2])
ax1.set_xlabel('time [s]')
ax1.set_ylabel('height [m]')

ax2 = fig.add_subplot(212)
line_3, = ax2.plot(dataArray[0], dataArray[2], 'r--', linewidth=3, label='velocityNumerical')
line_4, = ax2.plot(dataArrayAnalytic[0], dataArrayAnalytic[2], 'c', label='velocityAnalytic')
line_5, = ax2.plot(time, result, 'g', label='velocityODEsolver')
ax2.axis([0, 3.1, 0, 31])
ax2.legend(handles=[line_3, line_4, line_5])
ax2.set_xlabel('time [s]')
ax2.set_ylabel('velocity [m/s]')


plt.show()


# ALL IN ONE PLOT
# plt.plot(dataArray[0], dataArray[1], 'r--', dataArray[0], dataArray[2], 'y--', dataArrayAnalytic[0], dataArrayAnalytic[1], 'c--', dataArrayAnalytic[0], dataArrayAnalytic[2], 'g--')
# line_1, = plt.plot(dataArray[0], dataArray[1], 'r--', label='heightNumerical')
# line_2, = plt.plot(dataArray[0], dataArray[2], 'b--', label='velocityNumerical')
# line_3, = plt.plot(dataArrayAnalytic[0], dataArrayAnalytic[1], 'c--', label='heightAnalytic')
# line_4, = plt.plot(dataArrayAnalytic[0], dataArrayAnalytic[2], 'y--', label='velocityAnalytic')
# plt.axis([0, 3.1, 0, H+3])
# plt.xlabel('time [s]')
# plt.ylabel('height [m]\nvelocity [m/s]')
# plt.legend(handles=[line_1, line_2, line_3, line_4])
# plt.show()
