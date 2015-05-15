import numpy
import random
from scipy import sparse
from math import *
from scipy.sparse.linalg import spsolve
import numpy as np
import numpy, math, random
import numpy as np
from math import *
import matplotlib.pyplot as plt

L=100


size = 50 #total size of the peptone system
initWalkers = 25 #No. of initial walkers
time = 2000
peptoneConc = 50 #food

#walker parameters
reproThresh = 10
threshEnergy = 0.0
maxUptake   = 0.2
metabolism = 0.0667
jump = 0.4
InitEnergy = 0.33
reproEnergy = 0.30
Consuming_rate = 0.3
N_strikes_thresh = 6

def Initialize_Walkers(N,L):
  
  team=[]
  Nsquare = 1
  while(N > (Nsquare*Nsquare)):
    Nsquare += 1
  if(Nsquare**2 != N):
    print("CubicInit Warning: Your particle number",N, \
          "is not a perfect square; this may result " \
          "in a lousy initialization")
  rs = float(L)/Nsquare
  roffset = float(L)/2 - rs/2
  added = 0
  for x in range(0, Nsquare):
    for y in range(0, Nsquare):
        if(added < N):
          position = []
          position.append(rs*x+ roffset) 
          position.append(rs*y+ roffset) 
          position.append(reproThresh/3.0)
          added += 1
        team.append(position)

  print(team)
  return team


def Initialize_Peptone(array, size):
    for i in range(0,size):
        for j in range(0,size):
            array[i][j] = peptoneConc
    return array

D=1
dt=.1

h=L/float(size+1)



s=D*dt/2
c=1+2*s/h**2
b=-s/h**2
num=size
d_main = np.ones([num])*c
d_sub = np.ones([num])*b   
d_super = np.ones([num])*b     
data = [d_sub, d_main, d_super]   # list of all the data
diags = [-1,0,1]                  # which diagonal each vector goes into
A = sparse.spdiags(data,diags,num,num,format='csc')  # create the matrix


def Solve_PDE(C, size):
    n=size
    Y1=np.zeros([n,n])

  
    for j in range(n):
        C[0,j]=C[1,j]
        C[n-1,j]=C[n-2,j]
    
    for i in range(1,n-1):
        for j in range(1,n-1):
            Y1[i,j]=(C[i-1,j]-2*C[i,j]+C[i+1,j])*-b+C[i,j]
    for i in range(n):
        Y1[i,0]=Y1[i,1]
        Y1[i,n-1]=Y1[i,n-2]
    X1=np.zeros([n,n])
    
    
    for i in range(0,n):
        X1[i,:] = spsolve(A,Y1[i,:]) 
    
    for i in range(n):
        X1[i,0]=X1[i,1]
        X1[i,n-1]=X1[i,n-2]
  
        
        
    X2=np.zeros([n,n])
    for i in range(1,n-1):
        for j in range(1,n-1):
            X2[i,j]=(X1[i,j-1]-2*X1[i,j]+X1[i,j+1])*-b+X1[i,j]
    for j in range(n):
        X2[0,j]=X2[1,j]
        X2[n-1,j]=X2[n-2,j]
      
    Y2=np.zeros([n,n])
    
    
    
    
    for j in range(0,n):
        Y2[:,j] = spsolve(A,X2[:,j]) 

    for i in range(n):
        Y2[i,0]=Y2[i,1]
        Y2[i,n-1]=Y2[i,n-2]
    for j in range(n):
        Y2[0,j]=Y2[1,j]
        Y2[n-1,j]=Y2[n-2,j]
            
    
    return Y2
#-------------------------Simple Food bias along the direction of maximum food-----------------------#
def check_neighbors_peptone(actives,i,Peptone):
    x = int(actives[i][0])
    y = int(actives[i][1])
    d_x = 0
    d_y = 0
    if x==50: 
        d_x = 0
        d_y = 0
    elif x==0: 
        d_x = 0
        d_y = 0
    elif y==50: 
        d_x = 0
        d_y = 0
    elif y==0: 
        d_x = 0
        d_y = 0
    else: 
        if (Peptone[x][y] == Peptone[x][y-1] and Peptone[x][y] == Peptone[x][y+1] and Peptone[x][y] == Peptone[x+1][y] and Peptone[x][y] == Peptone[x-1][y]):
            d_x = -1
            d_x = -1
            
        if Peptone[x+1][y] > Peptone[x][y]:
            d_y = 0
            d_x = (Peptone[x+1][y])+(Peptone[x][y])/float(peptoneConc)
            Peptone[x][y] = Peptone[x+1][y]
    
        if Peptone[x][y+1] > Peptone[x][y]:
            d_x = 0
            d_y = (Peptone[x][y+1])+(Peptone[x][y])/float(peptoneConc)
            Peptone[x][y] = Peptone[x][y+1]
            
        if Peptone[x-1][y] > Peptone[x][y]:
            d_y = 0
            d_x = -(Peptone[x-1][y])+(Peptone[x][y])/float(peptoneConc)
            Peptone[x][y] = Peptone[x-1][y]
            
        if Peptone[x][y-1] > Peptone[x][y]:
            d_x = 0
            d_y = -(Peptone[x][y-1])+(Peptone[x][y])/float(peptoneConc)
            
        else:
            d_x = 0
            d_y = 0
    return d_x, d_y
#------------------------------Main Function----------------------------------------#

Peptone=np.zeros([size,size])
for i in range(0,size):
    for j in range(0,size):
        Peptone[i][j] =50*np.exp(-((i-25)**2+(j-25)**2)/169)



total_walkers = Initialize_Walkers(initWalkers, 5.0)

actives = []
inactives = []

for walker in total_walkers:
    if walker[2] <= threshEnergy:
        inactives.append(walker)
    
    else:
        actives.append(walker)

N_strikes = 0
active_list = []
inactive_list = []

for t in range (0, 2000):
    
    i = 0     
    while (i<len(actives)):
        #walker = nearest_grid(walker)
        x = int(actives[i][0])
        y = int(actives[i][1])
        if x>=49 or x<=10 or y>=49 or y<=1:
            d_x = 0 
            d_y =0
        else:
            d_x, d_y = check_neighbors_peptone(actives, i, Peptone)
        if (d_x == -1 and d_x == -1):
            dx = 0
            dy = 0
        else:
            #print actives[i][0]+dx
            #print actives[i][1]+dx
            theta = random.random()*(2*pi)
            #= random.random()*(2*pi)
            d = random.random()*jump
            dx = d*cos(theta) + d_x
            dy = d*sin(theta) + d_y
            #print actives[i][0]+dx
            #print actives[i][1]+dy
            #print 0.9*size
        if (actives[i][0]+dx >=size or actives[i][0]+dx <=0 or actives[i][1]+dy>=size or actives[i][1]+dy<=0):
            break
        else:
            actives[i][0] = actives[i][0] + dx
            actives[i][1] = actives[i][1] + dy
            
            x = int(actives[i][0])
            y = int(actives[i][1])
            availablefood = Peptone[x][y]
            food = min(availablefood, Consuming_rate)
            actives[i][2] = actives[i][2] + food - metabolism
            Peptone[x][y] -= food
        
            if actives[i][2] >= reproThresh:
                #print actives[i][2]
                child = []
                actives[i][2] = actives[i][2] - (reproEnergy + InitEnergy)
                child.append(actives[i][0] + 0.05)
                child.append(actives[i][1] + 0.05)
                child.append(InitEnergy)
                #print child
                actives.append(child)
                #print child
                #i = i+1
                #print "Not moving"
                
        
            if actives[i][2] <= threshEnergy:
                walker=actives[i]
                #print actives
                #print actives
                #x = len(actives)
                inactives.append(actives[i])
                actives.remove(actives[i])
                #y = len(actives)
               # z = len(inactives)
                #print inactives
                #i = i+1
                #print x,y,z
        i = i+1
    active_list.append(len(actives))
    inactive_list.append(len(inactives))
    Peptone = Solve_PDE(Peptone, size)
    #print Peptone 

print("Total number of actives is", len(actives))
print("Total number of inactives is", len(inactives))
act=np.array(actives)
plt.scatter(act[:,0], act[:,1],color='red')
inact=np.array(inactives)
plt.scatter(inact[:,0], inact[:,1],s=1,color='blue')

import numpy, math, random
import numpy as np
from math import *
import matplotlib.pyplot as plt

import timeit

time_array = []
for x in range(2000):
    time_array.append(x)
plt.plot(time_array,active_list,'g+')
plt.plot(time_array,inactive_list,'b^')
