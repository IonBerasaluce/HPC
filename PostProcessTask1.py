# Post-proccessing script for task 1 of the HPC coursework
# written by: Ion Berasaluce
from sys import argv

import matplotlib
import matplotlib.pyplot as plt


#reading the file outputs from the c++ source code and the details inputted by the user
# u,x = np.loadtxt('output.txt', unpack=True,
#                            delimiter=' ')

output = open("output.txt")

lines = output.readlines()

for i in range (1, len(argv)):
	if (argv[i] == '--length'):
		L = float(argv[i + 1])
	elif (argv[i] == '--no_elements'):
		N = float(argv[i + 1])
	elif (argv[i] == '--MomInertia'):
		I = float(argv[i + 1])
	elif (argv[i] == '--Emodulus'):
		E = float(argv[i + 1])

u = []
x = []

for line in lines:
	values = line.split(' ')
	u.append(float(values[0]))
	x.append(float(values[1]))

l = L/N

Fy = -1000
qy = -1000

#adding the boundary conditions 
x = [0] + x + [L]
u = [0] + u + [0]

ua = []

#Solving the analytical solution to plot against that obtained by the solver
for i in range(0, len(x)):
	if x[i] < L/2:
		ua.append(qy*x[i]*x[i]/(24 * E * I) * (L - x[i]) * (L - x[i]) + ((Fy * x[i] * x[i])/
			(48 * E * I) * (3*L - 4*x[i])))
	else: 
		ua.append(qy*x[i]*x[i]/(24 * E * I) * (L - x[i]) * (L - x[i]) + ((Fy * (L - x[i]) * (L - x[i]))/
			(48 * E * I) * (3*L - 4*(L - x[i]))))

#making the plots
plt.plot(x,ua,'--',x,u,'bs')
plt.ylabel('Deflection in m')
plt.xlabel('Length along the beam in metres')
plt.show()
