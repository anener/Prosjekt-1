import numpy as np
import matplotlib.pyplot as plt

n = 4
h = 1./(n+1)

#x verdiene
x = np.ones(n)
x[0] = 0
x[n-1] = 1
for i in range(1, n-1):
    x[i] = i*h

u = 1-(1-np.exp(-10))*x - np.exp(-10*x) #funksjonen
f = 100*np.exp(-10*x) #den deriverte funksjonen
u[n-1] = 0

"""
prosjekt 1b)
"""

"""
#dette er A matrisen
a = np.ones(n)*-1
b = np.ones(n)*2
c = np.ones(n)*-1

#dette er b=f(i)*(h**2)
w = f*(h*h)

#dette er v vektoren
v = np.ones(n)

#Gaussius eliminasjon

#Forward
bt = np.ones(n)
bt[1] = b[1]

wt = np.ones(n)
wt[1] = w[1]

for i in range(2, n):
    bt[i] = b[i] - (a[i-1]*c[i-1])/bt[i-1]
    wt[i] = w[i] - (a[i-1]*wt[i-1])/b[i-1]

#Backward
v[n-1] = wt[n-1]/bt[n-1]
for i in range(n-2, 1):
    v[i] = wt[i] - (c[i]*v[i+1])/bt[i]



plt.figure(1)
plt.subplot(211)
plt.plot(x, u, 'r--')
plt.subplot(212)
plt.plot(x, v, 'b--')
plt.show()
"""

"""
prosjekt 1c)
"""


#dette er A matrisen
b = np.ones(n)*2

#dette er b=f(i)*(h**2)
w = f*(h*h)

#dette er v vektoren
v = np.ones(n)

#Gaussius eliminasjon

#Forward
bt = np.ones(n)
bt[1] = b[1]

wt = np.ones(n)
wt[1] = w[1]

for i in range(2, n):
    bt[i] = b[i] - 1/bt[i-1]
    wt[i] = w[i] - (-wt[i-1])/b[i-1]

#Backward
v[n-1] = wt[n-1]/bt[n-1]
for i in range(n-2, 1):
    v[i] = wt[i] - (-v[i+1])/bt[i]



plt.figure(1)
plt.subplot(211)
plt.plot(x, u, 'r--')
plt.subplot(212)
plt.plot(x, v, 'b--')
plt.show()


"""
prosjekt 1d)
"""
erro = np.ones(n)
for i in range(1, n):
    erro[i] = np.log(abs(v[i] - u[i])/u[i])
