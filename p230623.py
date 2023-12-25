# cd c:\users\anatoliy\documents\science\grafit
# python p230623.py
#
import math
import numpy as np
import matplotlib.pyplot as plt
# Dispersion calculation
def disp(n,dx,dy):
 d=0
 q=p.copy()
 q[n][0]=p[n][0]+dx
 q[n][1]=p[n][1]+dy
 for i in range(len(sdss)):
# We take galaxy number i. We will look for closest vertice.
# Second space means that we check sdss galaxies.
  r1=10000
  for j in range(len(p)):
# Third space means that we check vertices to find closest filament.
   cx=sdss[i][0]-q[j][0]
   cy=sdss[i][1]-q[j][1]
   r2=cx*cx+cy*cy
   if r2<r1:
    r1=r2
    i1=j
# i1 is the number of vertice closest to galaxy
   r1=10000
   for k in range(len(conn)):
    if conn[k][0]==i1:
     i2=conn[k][1]
     r2=(sdss[i][0]-q[i2][0])**2+(sdss[i][1]-q[i2][1])**2
     if r2<r1:
      r1=r2
      k1=k
# k is the number of closest filament to galaxy.
# We will find distance from galaxy i to filament k.
  j1=conn[k1][0]
  j2=conn[k1][1]
  x0=q[j][0]
  y0=q[j][1]
  x1=q[j1][0]
  y1=q[j1][1]
  x2=q[j2][0]
  y2=q[j2][1]
  if i%1000==0:print(i,j1,j2)
  if j1!=j2:
   d=d+abs((x2-x1)*(y1-y0)-(x1-x0)*(y2-y1))/((x2-x1)**2+(y2-y1)**2)**0.5
 return d
 # Perevirka
def perevirka(n,dx,dy):
 global d1,dx1,dy1,n1
 a=0
 d2=disp(n,dx,dy)
 print("Dispersion is equal to ",d2)
 if d2<d1:
  d1=d2
  dx1=dx
  dy1=dy
  n1=n
  a=1
 return a
 # Main program
f=open("result.txt","a+")
f.write(" --- LSS approximation by graphs --- \n")
sdss = []
with open("sdss.txt", 'r') as file:
 for line in file:
   row = [float(num) for num in line.strip().split()]
   sdss.append(row)
# Vertice generation
p = []
r=10
dx=120
dy=2
for row in range(4):
 for col in range(4):
  x=3*r*col+dx
  y=1.7*r*row+dy
  p.append([x, y])
  x=3*r*col+2*r+dx
  y=1.7*r*row+dy
  p.append([x, y])
  x=3*r*col+0.5*r+dx
  y=1.7*r*(row+0.5)+dy
  p.append([x, y])
  x=3*r*col+1.5*r+dx
  y=1.7*r*(row+0.5)+dy
  p.append([x, y])
# Array of filaments
conn=[]
for i in range(len(p)):
 for j in range(len(p)):
  if (p[i][0]-p[j][0])**2+(p[i][1]-p[j][1])**2<r*r*2:
   conn.append((i,j))
#d1=disp(1,0,0)
dx=2
dy=2
dx1=0
dy1=0
skiki=len(p)
q=p.copy()
#skiki=2
# Skiki is number of clusters which we take into account
for i in range(skiki):
    xs=p[i][0]
    ys=p[i][1]
    n1=1
    x1,x2=xs-dx,xs+dx
    y1,y2=ys-dy,ys+dy
    for j in range(len(sdss)):
        if (x1<sdss[j][0] and sdss[j][0]<x2 and y1<sdss[j][1] and sdss[j][1]<y2):
            xs,ys=xs+sdss[j][0],ys+sdss[j][1]
            n1=n1+1
    q[i][0]=xs/n1
    q[i][1]=ys/n1
    print(x1,'<',q[i][0],'<',x2,y1,'<',q[i][1],'<',y2,'n1=',n1)
# Plotting
plt.scatter(*zip(*q), color='red')
plt.scatter(*zip(*sdss), color='black')
for i in range(len(conn)):
 n1=conn[i][0]
 n2=conn[i][1]
 plt.plot((q[n1][0],q[n2][0]),(q[n1][1],q[n2][1]),color='green')
plt.xlabel('RA')
plt.ylabel('DEC')
plt.ylim(-10,70)
plt.title('Hexagonal Grid Points')
plt.show()
#print(p)
#print(conn)
f.write("\n")
#f1.close()
f.close()
