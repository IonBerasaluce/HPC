# Post-proccessing script for task 2 of the HPC coursework
# written by: Ion Berasaluce
 
import matplotlib.pyplot as plt

#reading the file outputs from the c++ source code and the details inputted by the user

output = open("Out.txt")

lines = output.readlines()

maxval = []
minval = []
Amp = []
time = []

# Read the textfile which is separated by spaces and import the values into the corresponding arrays.
with open("Out.txt") as f:
    for line in f.readlines()[2: -1 if 2 else None]:
		cleanedLine = line.strip()
		values = cleanedLine.split(' ')
		maxval.append(float(values[0]))
		minval.append(float(values[1]))
		time.append(int(values[2]))

for i in range(len(maxval)):
	Amp.append(0.5*(maxval[i] - minval[i]))

#making the plots
plt.loglog(time,Amp)
plt.show()
