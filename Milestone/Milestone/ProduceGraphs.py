import numpy
import matplotlib.pyplot as plt
import os
import sys
import math

def onsager(J):
    kbt = 1
    zcrit = -1 + math.sqrt(2)
    z = math.exp((-2*J)/kbt)
    if (z < zcrit):
        M = ((1+(z**2))**0.25)*((1-6*(z**2)+(z**4))**0.125)/((1-(z**2))**0.5)
        return M
    else:
        z = 2* zcrit - z
        M = 1-  ((1+(z**2))**0.25)*((1-6*(z**2)+(z**4))**0.125)/((1-(z**2))**0.5)
        return M

script_dir = os.path.dirname(os.path.abspath(__file__))
file_path = os.path.join(script_dir, sys.argv[1])
f = open(file_path,"r")
Ndata = int(f.readline())
xLabel = f.readline()
yLabel = f.readline()
xs = numpy.empty(Ndata)
ys = numpy.empty(Ndata)
yErrors = numpy.empty(Ndata)
for i in range (0, Ndata):
    data = f.readline()
    xs[i], ys[i],yErrors[i] = data.split(",")

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
#ax.spines['bottom'].set_position('center')
ax.spines['left'].set_position(('data', 0))
ax.spines["left"].set_color("white")
ax.spines["bottom"].set_color("white")
#ax.set_ylim([-1,1])
#ax.set_ylabel(yLabel)
ax.set_ylim([0,1])
ax.set_ylabel("Absolute Magnetisation",color = "white")
ax.set_xlabel(xLabel,color = "white")
ax.tick_params(axis='x', colors='white') 
ax.tick_params(axis='y', colors='white') 
plt.errorbar(xs,ys,yerr = yErrors,fmt = "xr" ,ecolor = "r",capsize = 2,label = "Simulation")


analytical = numpy.empty(Ndata)
for i in range(0,Ndata):
    analytical[i] = onsager(xs[i])
plt.plot(xs,analytical,color = "g",label ="Analytical")
plt.legend(loc='center right')
fig.savefig("VaryJ.png",transparent = "True")
plt.show()




