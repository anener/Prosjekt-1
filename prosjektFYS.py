# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

file1 = open("xu_n1000.txt", "r")
file2 = open("prosjektB_n1000.txt", "r")

nline = file1.readline()
narray = nline.split()
n = narray[-1]

hline = file1.readline()
harray = hline.split()
h = harray[-1]

file2.readline()
file2.readline()
tline = file2.readline()
tarray = tline.split()
t = tarray[-1]


file1.readline()
u = []
x = []
for line in file1:
    y = line.split(" ")
    u.append(y[-1])
    x.append(y[-2])
file1.close()

file2.readline()
file2.readline()
v = [0]
for line in file2:
    y = line.split(" ")
    v.append(y[-1])
file2.close()


plt.figure
plt.subplot(211)
plt.plot(x, u, 'b-')
plt.ylabel('u(x)')
plt.xlabel('x')
plt.subplot(212)
plt.plot(x, v, 'r-')
plt.ylabel('v(x)')
plt.xlabel('x')
plt.suptitle('u(x)/v(x) plot med n='+n+' og h='+h)
plt.savefig('xB_n1000.png')
plt.show()
