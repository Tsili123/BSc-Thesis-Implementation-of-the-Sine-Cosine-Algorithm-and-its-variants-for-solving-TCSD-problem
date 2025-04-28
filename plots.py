import matplotlib.pyplot as plt
import numpy as np
#pip3 install pandas
import pandas as pd
from mpl_toolkits import mplot3d


def logistic_map(a, x):
    return a * x * (1 - x)

l = []

x = 0.87
for i in range(0,100):
    l.append(x)
    x = logistic_map(4,x)
plt.plot(range(0,100), l, 'o-')
plt.title("Chaotic sequence x[k+1] = a * x[k] * (1 - x[k]) (first 100 terms , x[0] = 0.87)")
plt.savefig("./plots/plots.png")

plt.clf()

x = np.arange(-4*np.pi,4*np.pi,0.1)   # start,stop,step
y = 2*np.sin(x)
z = 2*np.cos(x)
plt.axvline(0, color='black', linewidth=.5)
plt.axhline(0, color='black', linewidth=.5)
plt.axhline(-2, color='black', linewidth=.5)
plt.axhline(2, color='black', linewidth=.5)
plt.axhline(-1, color='black', linewidth=.5)
plt.axhline(1, color='black', linewidth=.5)
plt.plot(x,y, label='y=2sin(x)',color = 'red')
plt.plot(x,z, label='y=2cosx(x)',color = 'black')
plt.legend(loc='upper left')
y1 = 1
y2 = -1
plt.axhspan(y1, y2, color='orange', alpha=0.75, lw=0)
y1 = -1
y2 = -2
plt.axhspan(y1, y2, color='blue', alpha=0.75, lw=0)
y1 = 1
y2 = 2
plt.axhspan(y1, y2, color='blue', alpha=0.75, lw=0)
plt.savefig("./plots/plots2.png")

plt.clf()
 
th = np.linspace( 0 , 8 * np.pi , 3000 )
 
r = 2
 
cosine = r * np.cos( th ) 
sine = r * np.sin( th ) 

r = 1

cosine2 = r * np.cos( th ) 
sine2 = r * np.sin( th ) 
 
plt.plot(sine,cosine, 'c-', label='circle')
plt.fill_between(sine, 0, cosine)
plt.plot(sine2,cosine2, 'c-', label='circle')
plt.fill_between(sine2, 1, cosine2)
plt.axvline(0, color='black', linewidth=.5)
plt.axhline(0, color='black', linewidth=.5)

# annotation of the third point
plt.scatter([0], [0], color="green") # plotting single point
plt.text(0,0.1,"Destination")
plt.scatter([-1], [0], color="white") # plotting single point
plt.text(-1.2,0.1,"Solution")


plt.title( 'Parametric Circle' )
plt.savefig('./plots/circle.png')

plt.clf()

import math

time = 3
samples_ratio = 0.002
x = np.arange(0, time, samples_ratio)
wave = 2*np.pi
rvar = np.flip(x, 0)

el2 = rvar * np.sin(2*wave * x)
plt.plot(x, el2,color='red',label='y=r1*sin(x)')

el = rvar * np.cos(2*wave * x)
plt.plot(x, el,color='black',label='y=r1*cos(x)')
plt.legend(loc='upper left')
plt.axhline(0, color='black', linewidth=.5)
plt.xlabel('time')
plt.ylabel('Range')
plt.title("Decreasing r1")

plt.savefig('./plots/customwave.png')


plt.clf()

samples = 8500
a = np.linspace(1.5, 4.0, samples)
loops = 2500
x = 1e-7 * np.ones(samples)

_, ax1 = plt.subplots(1, 1, figsize=(7, 7))

for i in range(loops):
    x = logistic_map(a, x)
    # display the bifurcation diagram.
    if i >= (loops-100):
        ax1.plot(a, x, ',k' ,alpha = 0.3)

ax1.set_xlim(1.5, 4)
ax1.set_title("Bifurcation diagram")

plt.savefig('./plots/chaos.png')

plt.clf()

df = pd.read_csv("sinecosine_f4_iter.csv")
#print(df.iloc[:,0])
p = df.iloc[:,0]


# Create an empty list
column_list =[]
  
# Iterate over each row
for  rows in p:
    # append the list to the final list
    column_list.append(rows)

y = [i for i in range(0,len(column_list[0:1999]))]#[0:100]
plt.step(y,column_list[0:1999], where='post')#[0:100]
plt.title("Aquila [ dimension = 100 ,  test function 4]  ")
plt.ylabel("Min fitness")
plt.xlabel("Iterations")

plt.savefig('./plots/convergence_aquila_100.png')

plt.clf()
from matplotlib import cm

def plot_ackley_3dimension():
    figure = plt.figure()
    axis = figure.add_subplot(projection='3d')

    x = np.arange(-32, 32, 0.15)
    y = np.arange(-32, 32, 0.15)
    x, y = np.meshgrid(x, y)

    m = 2 * np.pi
    k = 20
    l = 0.2

    cosine_t = -np.exp((np.cos(m*x) + np.cos(m*y)) / 2)
    sqrt_sum_t = -k * np.exp(-l * np.sqrt(x*x + y*y) / 2)
    z =  sqrt_sum_t + cosine_t +  k + np.exp(1)

    surface = axis.plot_surface(x, y, z, cmap=cm.hsv)

    figure.colorbar(surface, shrink=0.7, aspect=5)

    plt.savefig('./plots/ackley_3d.png')

plt.clf()
import math
def plot_rastrigin_3dimension():
    figure = plt.figure()
    axis = figure.add_subplot(projection='3d')

    x = np.arange(-5.12, 5.12, 0.30)
    y = np.arange(-5.12, 5.12, 0.30)
    x, y = np.meshgrid(x, y)

    a = 10
    c = 2*np.pi

    cos_term = -10 * (np.cos(c*x) + np.cos(c*y))
    x_i = x*x+y*y
    z = a +  cos_term  +x_i

    surface = axis.plot_surface(x, y, z, cmap=cm.hsv)

    figure.colorbar(surface, shrink=0.7, aspect=5)

    axis.view_init(20, 50)

    plt.savefig('./plots/rastrigin_3d.png')

plt.clf()
def plot_dejongs_3dimension():
    figure = plt.figure()
    axis = figure.add_subplot(projection='3d')

    x = np.arange(-100, 100, 0.15)
    y = np.arange(-100, 100, 0.15)
    x, y = np.meshgrid(x, y)

    z = x*x+y*y

    surface = axis.plot_surface(x, y, z, cmap=cm.hsv)

    figure.colorbar(surface, shrink=0.7, aspect=5)

    plt.savefig('./plots/dejongs_3d.png')

plt.clf()

def plot_spring_3d():
    
    figure = plt.figure(figsize = (7,7))
    axis = plt.axes(projection = '3d')
    
    z = np.linspace(0, 80, 500)
    x = np.sin(z)
    y = np.cos(z)
    axis.plot3D(x, y, z, 'blue')

    axis.view_init(-50,90)

    plt.savefig('./plots/spring.png')

plot_ackley_3dimension() 
plot_rastrigin_3dimension() 
plot_dejongs_3dimension() 
plot_spring_3d() 
  


