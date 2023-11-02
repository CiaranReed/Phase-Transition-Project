import numpy
import matplotlib.pyplot as plt
import os


script_dir = os.path.dirname(os.path.abspath(__file__))
file_path = os.path.join(script_dir, "magnetisation data.txt")
f = open(file_path,"r")
Ndata = int(f.readline())
shistory = numpy.empty(Ndata)
mhistory = numpy.empty(Ndata)
for i in range (0, Ndata):
    data = f.readline()
    shistory[i], mhistory[i] = data.split(",")



fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.spines['bottom'].set_position('center')
ax.spines['left'].set_position(('data', 0))
ax.set_ylim([-1,1])
ax.set_ylabel("Magnetisation")
ax.set_xlabel("Number of Sweeps")

plt.plot(shistory,mhistory)
plt.show()
