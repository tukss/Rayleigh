from diagnostic_reading import ShellAverage
import matplotlib.pyplot as plt
a = ShellAverage(filename='01330000',path='Shell_Avgs/')
rsun = 6.96e10
# The quantity code for entropy is 4
# We can use a.lut[4] to find the index that this 
# quantity is stored at
entropy_index = a.lut[4]
time_index = 0
plt.figure(1)
plt.plot(a.radius/rsun,a.vals[:,entropy_index,time_index])
plt.xlabel('Radius r/R'+r'$_\odot$')
plt.ylabel('Specific Entropy '+r'erg g$^{-1}$ K$^{-1}$')
#plt.show()
plt.savefig('shell_average.png')
