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
plt.plot(shistory,mhistory,color = "red")


ax.xaxis.label.set_color('white')        
ax.yaxis.label.set_color('white')     
ax.spines['bottom'].set_color('white')
ax.spines['left'].set_color('white')
ax.tick_params(axis='x', colors='white')    
ax.tick_params(axis='y', colors='white') 
extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
#fig.savefig('graph1.png', bbox_inches=extent.expanded(1.3,1.4), transparent= True)
plt.tight_layout()
fig.savefig('graph2.png',transparent= True,bbox_inches='tight')
ax.xaxis.label.set_color('black')        
ax.yaxis.label.set_color('black')     
ax.spines['bottom'].set_color('black')
ax.spines['left'].set_color('black')
ax.tick_params(axis='x', colors='black')    
ax.tick_params(axis='y', colors='black') 
plt.show()



