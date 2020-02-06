import numpy as np
import matplotlib.pyplot as plt

"""
Script to read and plot COM position
"""
com = open('COM.txt', 'r')
com_data = []
for line in com:
    com_data.append(line)

com_time = np.array([])
com_x = np.array([])
com_y = np.array([])
com_z = np.array([])

for i in com_data:
    com_time = np.append(com_time, (i.split(('\t\t'))[0]))
    com_x = np.append(com_x, (i.split(('\t\t'))[1]))
    com_y = np.append(com_y, (i.split(('\t\t'))[2]))
    com_z = np.append(com_z, (i.split(('\t\t'))[3]))

com_time = np.delete(com_time, 0)
com_x = np.delete(com_x, 0)
com_y = np.delete(com_y, 0)
com_z = np.delete(com_z, 0)
    
com_time = com_time.astype(np.float)
com_x = com_x.astype(np.float)
com_y = com_y.astype(np.float)
com_z = com_z.astype(np.float)

#Plots the x, y, and z coordinates of the COM vs. sim time
fig, axs = plt.subplots(2, 2)
fig.suptitle('COM Position')
axs[0, 0].plot(com_time, com_x)
axs[0, 0].set_title('$x$-component')
axs[0, 1].plot(com_time, com_y, 'tab:orange')
axs[0, 1].set_title('$y$-component')
axs[1, 0].plot(com_time, com_z, 'tab:green')
axs[1, 0].set_title('$z$-component')

for ax in axs.flat:
    ax.set(xlabel='Time (yr)', ylabel='Position (AU)')

# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()

plt.show()