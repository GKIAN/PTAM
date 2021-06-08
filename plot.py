#!/bin/env python

import numpy as np
import scipy.special as ssp
import matplotlib.pyplot as plt
import math

dk = 0.01
kl = 28.5

a = 0.001
b = 3.00

nk = math.ceil(kl / dk)

k = np.arange(nk) * dk
t = np.exp( - a * k) * ssp.jv(1, b * k) * dk
t = np.cumsum(t)

plt.plot(k, t)
plt.show()

