# Post-proccessing script for task 2 of the HPC coursework
# written by: Ion Berasaluce

import matplotlib
import matplotlib.pyplot as plt


#reading the file outputs from the c++ source code and the details inputted by the user
# u,x = np.loadtxt('output.txt', unpack=True,
#                            delimiter=' ')

output = open("Out.txt")

lines = output.readlines()

u = []
t = []

with open("Out.txt") as f:
	#remove the header and footer of the text-file so the right values are added into the arrays
    for line in f.readlines()[2: -1 if 2 else None]: 
		cleanedLine = line.strip()
		values = cleanedLine.split(' ')
		u.append(float(values[0]))
		t.append(float(values[1]))

#making the plots
plt.plot(t,u)
plt.show()



