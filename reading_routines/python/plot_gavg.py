###############################################
#
#  Global-Averages (G_Avgs) plotting example
#  Reads in time steps 0 through 1.3 million
#  Plots average KE vs time

from diagnostic_reading import GlobalAverage, build_file_list
import matplotlib.pyplot as plt
files = build_file_list(0,1300000,path='G_Avgs')

plt.figure(1)
alldays = []
ke = []
for f in files:
    a = GlobalAverage(filename=f,path='')
    # The quantity code for kinetic energy is 7
    # We can use a.lut[4] to find the index that this 
    # quantity is stored at
    ke_index = a.lut[7]
    days = a.time/(3600.0*24.0)
    for i,d in enumerate(days):
        alldays.append(d)
        ke.append(a.vals[i,ke_index])



plt.plot(alldays,ke)




plt.xlabel('Time (days)')
plt.ylabel('Mean Kinetic Energy Density '+r'erg g$^{-1}$ cm$^{-3}$')
#plt.show()
plt.savefig('global_average.png')
