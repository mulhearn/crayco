#! /Users/mulhearn/anaconda/bin/python

import numpy as np
import matplotlib.pyplot as plt

x, y, h = np.loadtxt("movie/frame1.txt",unpack=True)

print x.size
print y.size
print h.size

print x[h > 0]


fig = plt.figure()
ax1 = fig.add_subplot(111)

ax1.set_title("Smartphone Data")    
ax1.set_xlabel('x')
ax1.set_ylabel('y')

ax1.scatter(x[h>0],y[h>0], c='r', label='hit')
ax1.scatter(x[h==0],y[h==0], c='b', label='no hit')

leg = ax1.legend()

plt.show()
 
