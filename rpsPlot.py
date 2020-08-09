import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap


my_colors = ['white', 'royalblue', 'red', 'magenta']
cmap1 = LinearSegmentedColormap.from_list('', my_colors, len(my_colors))


initial = np.loadtxt('initials')

N = int(initial[0,0])
mcstep1 = initial[:,1]
mcstep2 = initial[:, 2]
mcstep3 = initial[:, 3]
mcstep4 = initial[:,4]
d = initial[:, 5]
r = initial[:, 6]
p = initial[:, 7]


data = np.loadtxt('density-time-astable') 
data0 = np.loadtxt('lattice0', delimiter = ',')
data1 = np.loadtxt('lattice1', delimiter = ',')
data2 = np.loadtxt('lattice2', delimiter = ',')
data3 = np.loadtxt('lattice3', delimiter = ',')
data4 = np.loadtxt('lattice4', delimiter = ',')


lattice0= np.reshape(data0, (N,N))
lattice1= np.reshape(data1, (N,N))
lattice2= np.reshape(data2, (N,N))
lattice3= np.reshape(data3, (N,N))
lattice4= np.reshape(data4, (N,N))

time = data[:, 0]
rhoA = data[:, 1]
rhoB = data[:, 2]
rhoC = data[:, 3]






def densityPlot():
	plt.plot(time, rhoA, label = "A", color= "royalblue")
	plt.plot(time, rhoB, label = "B", color= "red")
	plt.plot(time, rhoC, label = "C", color= "magenta")
	plt.title("Density vs time graph")
	plt.legend()
	plt.grid()
	
	plt.title("Evolutionary RPS Model \n" + "Lattice Size:" + "(" + str(N) + "X" + str(N) + ")" + "\nTotal MC Step: " + str(
        mcstep4[0]) + "\nNatural Death Rates: " + str(
        d) + "\nReproduction rates: " + str(r) + "\nPredation Rate: " + str(p))
	plt.show()
	
	plt.show()
	
	

def latticePlot():
	plt.subplot(2,2,1)
	plt.imshow(lattice0, cmap= cmap1)
	plt.title("MCStep= " + str(0))
	
	
	plt.subplot(2,2,2)
	plt.imshow(lattice1, cmap= cmap1)
	plt.title("MCStep= " + str(mcstep1[0]))
	
	plt.subplot(2,2,3)
	plt.imshow(lattice2, cmap= cmap1)
	plt.title("MCStep= " + str(mcstep2[0]))
	
	plt.subplot(2,2,4)
	plt.imshow(lattice4, cmap= cmap1)
	plt.title("MCStep= " + str(mcstep4[0]))
	
	plt.suptitle("Evolutionary RPS Model \n" + "Lattice Size:" + "(" + str(N) + "X" + str(N) + ")" + "\nTotal MC Step: " + str(
        mcstep4[0]) + "\nNatural Death Rates: " + str(
        d) + "\nReproduction rates: " + str(r) + "\nPredation Rate: " + str(p))
	plt.show()
	

def generalPlot():
	plt.subplot(2,2,1)
	plt.imshow(lattice0, cmap= cmap1)
	plt.title("MCStep= " + str(0))
	
	
	plt.subplot(2,2,2)
	plt.imshow(lattice2, cmap= cmap1)
	plt.title("MCStep= " + str(mcstep2[0]))	
	
	
	plt.subplot(2,2,4)
	plt.imshow(lattice4, cmap= cmap1)
	plt.title("MCStep= " + str(mcstep4[0]))
	
	plt.subplot(2,2,3)
	plt.plot(time, rhoA, label = "A", color= "royalblue")
	plt.plot(time, rhoB, label = "B", color= "red")
	plt.plot(time, rhoC, label = "C", color= "magenta")
	plt.title("Density vs time graph")
	plt.legend()
	plt.grid()
	
	
	
	plt.suptitle("Evolutionary RPS Model \n" + "Lattice Size:" + "(" + str(N) + "X" + str(N) + ")" + "\nTotal MC Step: " + str(
        mcstep4[0]) + "\nNatural Death Rates: " + str(
        d) + "\nReproduction rates: " + str(r) + "\nPredation Rate: " + str(p))
	plt.show()


latticePlot()	
densityPlot()	
generalPlot()
	
	
	
	
	
	
	
	

