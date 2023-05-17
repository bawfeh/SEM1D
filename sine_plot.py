# Numpy tutorial: basics
import numpy as np
from numpy import pi

import matplotlib.pyplot as plt

x = np.linspace( 0, 2*pi, 10 )
# useful to evaluate function at lots of points

f = np.sin(x)

#print(x)

plt.plot(x, f)
plt.xlabel("x")
plt.ylabel("y")
# str = input("Enter figure ID: ")
# plt.title("Figure %s:\n The sine curve [y = sin(x)] using numpy's matplotlib." % str)
plt.title("The sine curve [y = sin(x)]\n, using scipy's (or numpy's) matplotlib.")

plt.show()
