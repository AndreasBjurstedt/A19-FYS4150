#A little testfile
import numpy as np
import matplotlib.pyplot as plt

x = np.arange(-4,4.2,0.2)

y = x**2
z = x**3

fig = plt.figure
ax1 = plt.subplot(211)
ax1.plot(x,y)

ax2 = plt.subplot(212)
ax2.plot(x,z)
plt.show()

print(x)