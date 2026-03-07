import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
mpl.rc_file("mplstyleerc")


nopart_10k = np.genfromtxt(f'../cycles-10k-nopart.dat', delimiter=',')
part_10k = np.genfromtxt(f'../cycles-10k-part.dat', delimiter=',')
nopart_100k = np.genfromtxt(f'../cycles-100k-nopart.dat', delimiter=',')
part_100k = np.genfromtxt(f'../cycles-100k-part.dat', delimiter=',')
nopart_1m = np.genfromtxt(f'../cycles-1m-nopart.dat', delimiter=',')
part_1m = np.genfromtxt(f'../cycles-1m-part.dat', delimiter=',')
part25_1m = np.genfromtxt(f'../cycles-1m-part25.dat', delimiter=',')

fig, ax = plt.subplots(3, 1, figsize=(6, 9))

# plot 10k
ax[0].plot(nopart_10k[:, 0], nopart_10k[:, 1],
           label='No Partitioning', color='C2')
ax[0].plot(part_10k[:, 0], part_10k[:, 1],
           label='Partitioning', color='#92D050')
ax[0].vlines(x=nopart_10k[-1, 0], ymin=0, ymax=nopart_10k[-1, 1], color='C2')
ax[0].vlines(x=part_10k[-1, 0], ymin=0, ymax=part_10k[-1, 1], color='#92D050')
ax[0].text(nopart_10k[-1, 0] + 0.1, nopart_10k[-1, 1]/2,
           f'{nopart_10k[-1, 0]:.1f}s', rotation=90, va='center', ha='left', color='C2')
ax[0].text(part_10k[-1, 0] + 0.1, part_10k[-1, 1]/2,
           f'{part_10k[-1, 0]:.1f}s', rotation=90, va='center', ha='left', color='#92D050')
ax[0].set_title('10,000 Events')
ax[0].set_xlabel('Time (seconds)')
ax[0].set_ylabel('Completed Events')
ax[0].legend()
ax[0].grid(alpha=0.3)

# Add zoomed inset to first plot
axins = inset_axes(ax[0], width="15%", height="20%",
                   loc='lower center', borderpad=2)
axins.plot(nopart_10k[:, 0], nopart_10k[:, 1], color='C2')
axins.plot(part_10k[:, 0], part_10k[:, 1], color='#92D050')
# Set the zoom region (adjust these values to highlight your desired area)
x1, x2, y1, y2 = 1.15, 1.23, 4900, 5100  # Example zoom region
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
axins.grid(alpha=0.3)
axins.set_xticks([])
axins.set_yticks([])
# Draw rectangle and connecting lines
mark_inset(ax[0], axins, loc1=2, loc2=3, fc="none", ec="0.5")

# plot 100k
ax[1].plot(nopart_100k[:, 0], nopart_100k[:, 1],
           label='No Partitioning', color='C2')
ax[1].plot(part_100k[:, 0], part_100k[:, 1],
           label='Partitioning', color='#92D050')
ax[1].vlines(x=nopart_100k[-1, 0], ymin=0, ymax=nopart_100k[-1, 1], color='C2')
ax[1].vlines(x=part_100k[-1, 0], ymin=0,
             ymax=part_100k[-1, 1], color='#92D050')
ax[1].text(nopart_100k[-1, 0] + 0.1, nopart_100k[-1, 1]/2,
           f'{nopart_100k[-1, 0]:.1f}s', rotation=90, va='center', ha='left', color='C2')
ax[1].text(part_100k[-1, 0] + 0.1, part_100k[-1, 1]/2,
           f'{part_100k[-1, 0]:.1f}s', rotation=90, va='center', ha='left', color='#92D050')
ax[1].set_title('100,000 Events')
ax[1].set_xlabel('Time (seconds)')
ax[1].set_ylabel('Completed Events')
ax[1].legend()
ax[1].grid(alpha=0.3)

# plot 1m
ax[2].plot(nopart_1m[:, 0], nopart_1m[:, 1],
           label='No Partitioning', color='C2')
ax[2].plot(part_1m[:, 0], part_1m[:, 1], label='Partitioning', color='#92D050')
ax[2].plot(part25_1m[:, 0], part25_1m[:, 1],
           label='Partitioning, Cutoff = 25k', color='#92D050', linestyle='--')
ax[2].vlines(x=nopart_1m[-1, 0], ymin=0, ymax=nopart_1m[-1, 1], color='C2')
ax[2].vlines(x=part_1m[-1, 0], ymin=0, ymax=part_1m[-1, 1], color='#92D050')
ax[2].vlines(x=part25_1m[-1, 0], ymin=0, ymax=part25_1m[-1, 1],
             color='#92D050', linestyle='--')
ax[2].text(nopart_1m[-1, 0] + 0.5, nopart_1m[-1, 1]/2,
           f'{nopart_1m[-1, 0]:.1f}s', rotation=90, va='center', ha='left', color='C2')
ax[2].text(part_1m[-1, 0] + 0.5, part_1m[-1, 1]/2,
           f'{part_1m[-1, 0]:.1f}s', rotation=90, va='center', ha='left', color='#92D050')
ax[2].text(part25_1m[-1, 0] + 0.5, part25_1m[-1, 1]/2,
           f'{part25_1m[-1, 0]:.1f}s', rotation=90, va='center', ha='left', color='#92D050')
ax[2].set_title('1,000,000 Events')
ax[2].set_xlabel('Time (seconds)')
ax[2].set_ylabel('Completed Events')
ax[2].legend()
ax[2].grid(alpha=0.3)

# Save the figure
fig.tight_layout()
fig.savefig('partitoning-results.pdf')
