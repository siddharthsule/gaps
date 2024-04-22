import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rc_file("../../test/mplstyleerc")

data = np.genfromtxt('acceptveto.csv', delimiter=',')

data = data[data[:,1] > 0]

x = data[:,0]
n = data[:,1]
y1 = data[:,2]
y2 = data[:,3]

fig, ax = plt.subplots(1, 2, figsize=(9, 3.75))

ax[0].plot(x, y1, label='Vetoed Emissions', color='C3')
ax[0].plot(x, y2, label='Accepted Emissions', color='C2')

ax[1].plot(x, y1/n, label='Vetoed Emissions', color='C3')
ax[1].plot(x, y2/n, label='Accepted Emissions', color='C2')

for i in range(2):
    ax[i].set_xlabel('Number of Events')
    ax[i].set_ylabel('Emissions')
    ax[i].legend()
    ax[i].grid(True)

ax[0].set_title('Accepted and Vetoed Emissions')
ax[1].set_title('Accepted and Vetoed Emissions, (Rescaled by $N_{Active Events}$)')

fig.tight_layout()
fig.savefig('acceptveto.pdf')
fig.savefig('acceptveto.png')