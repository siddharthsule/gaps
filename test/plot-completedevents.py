import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import matplotlib as mpl
mpl.rc_file("mplstyleerc")

# Plotting Number of Completed Events per cycle and Number of newly Completed
# Events per cycle

# Plot
fig, ax = plt.subplots(1, 2, figsize=(9, 3.75))

nev = [1000, 10000, 100000, 1000000]

for n in nev:

    filename = "../results-events/gaps-cycles-" + str(n)

    temp = np.genfromtxt(filename + ".dat", delimiter='\n')
    temp /= n  # Divide by number of events

    comp = np.zeros((200))
    diff = np.zeros((200))
    max = 0

    comp[:len(temp)] = temp
    comp[len(temp):] = temp[-1]

    for i in range(1, len(temp) - 1):

        diff[i] = temp[i] - temp[i-1]

    if len(temp) > max:
        max = len(temp)

    comp = comp[:max]
    comp = np.append(comp, 0.)

    diff = diff[:max]
    diff = np.append(diff, 0.)

    cycles = np.arange(1, len(comp) + 1)

    ax[0].scatter(cycles, comp, label=str(n) + '  Events ')
    ax[0].plot(cycles, comp, alpha=0.5)
    ax[0].set_xlabel('Cycle')
    ax[0].set_ylabel('Number  of  Completed  Events / Total')
    ax[0].set_title('Number  of  Completed  Events  per  Cycle')
    ax[0].legend()
    ax[0].grid(True)

    ax[1].scatter(cycles, diff, label=str(n) + '  Events ')
    ax[1].plot(cycles, diff, alpha=0.5)
    ax[1].set_xlabel('Cycle')
    ax[1].set_ylabel('Number  of  Newly  Completed  Events / Total')
    ax[1].set_title('Number  of  Newly  Completed   Events  per  Cycle')
    ax[1].legend()
    ax[1].grid(True)

fig.tight_layout()
plt.savefig("completedEvents.pdf")
