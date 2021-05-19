# -*- coding: utf-8 -*-
"""
Created on Tue May  11 10:53:16 2021

@author: axelb
"""

import numpy as np
from matplotlib.pyplot import *

#Diférents schémas :

#Explicite centrée
def constructionmatriceEC(N,c):
    A = np.eye(N)+c/2*np.eye(N,N,-1)-c/2*np.eye(N,N,1)
    A[0,N-1] = c/2
    A[N-1,0] = -c/2
    return(A)
    
#Implicite centrée
def constructionmatriceIC(N,c):
    A = np.eye(N)+(c/2)*np.eye(N,N,1)+(-c/2)*np.eye(N,N,-1)
    A[0,N-1] = -c/2
    A[N-1,0] = c/2
    return(A)

#Explicite décentrée Amont
def constructionmatriceEDA(N,c):
    A = (1-c)*np.eye(N)+c*np.eye(N,N,-1)
    A[0,N-1] = c
    return(A)

#Lax-Friedrichs
def constructionmatriceLaxF(N,c):
    A = ((1+c)/2)*np.eye(N,N,-1)+((1-c)/2)*np.eye(N,N,1)
    A[0,N-1] = ((1+c)/2)
    A[N-1,0] = ((1-c)/2)
    return(A)

#Lax-Wendroff
def constructionmatriceLaxW(N,c):
    A = (1-c**2)*np.eye(N)+((c**2-c)/2)*np.eye(N,N,1)+((c**2+c)/2)*np.eye(N,N,-1)
    A[0,N-1] = ((c**2+c)/2)
    A[N-1,0] = ((c**2-c)/2)
    return(A)

def construitX(N):
    return(np.linspace(0,1,N))

def U0(x):
    return(np.sin(np.pi*x)**10)

#Solutions des diférents schémas :
    
#Explicite centrée    
def SolutionEC(h,N,tau,c):
    ntfinal=int(Tmax/tau)+1
    A = constructionmatriceEC(N,c)
    UT = np.zeros((ntfinal,N))
    T = np.arange(ntfinal)*tau
    X = construitX(N)
    U = U0(X[0:N])
    UT[0,0:N] = U
    for n in range(ntfinal-1) :
        U = A@U
        UT[n+1,0:N] = U
    return X,T,UT

#Implicite centrée
def SolutionIC(h,N,tau,c):
    ntfinal=int(Tmax/tau)+1
    A = constructionmatriceIC(N,c)
    UT = np.zeros((ntfinal,N))
    T = np.arange(ntfinal)*tau
    X = construitX(N)
    U = U0(X[0:N])
    UT[0,0:N] = U
    for n in range(ntfinal-1) :
        U = A@U
        UT[n+1,0:N] = U
    return X,T,UT

#Explicite décentrée Amont
def SolutionEDA(h,N,tau,c):
    ntfinal=int(Tmax/tau)+1
    A = constructionmatriceEDA(N,c)
    UT = np.zeros((ntfinal,N))
    T = np.arange(ntfinal)*tau
    X = construitX(N)
    U = U0(X[0:N])
    UT[0,0:N] = U
    for n in range(ntfinal-1) :
        U = A@U
        UT[n+1,0:N] = U
    return X,T,UT

#Lax-Friedrichs
def SolutionLaxF(h,N,tau,c):
    ntfinal=int(Tmax/tau)+1
    A = constructionmatriceLaxF(N,c)
    UT = np.zeros((ntfinal,N))
    T = np.arange(ntfinal)*tau
    X = construitX(N)
    U = U0(X[0:N])
    UT[0,0:N] = U
    for n in range(ntfinal-1) :
        U = A@U
        UT[n+1,0:N] = U
    return X,T,UT

#Lax-Wendroff
def SolutionLaxW(h,N,tau,c):
    ntfinal=int(Tmax/tau)+1
    A = constructionmatriceLaxW(N,c)
    UT = np.zeros((ntfinal,N))
    T = np.arange(ntfinal)*tau
    X = construitX(N)
    U = U0(X[0:N])
    UT[0,0:N] = U
    for n in range(ntfinal-1) :
        U=A@U
        UT[n+1,0:N]=U
    return X,T,UT

#Définitions des t, h et tau
Tmax=2
Cas = [[0.02,0.01],[0.002,0.005],[0.002,0.002],[0.005,0.0002]]
Lt = [0,1,2]

#Explicite centré

fig,(axe1,axe2,axe3,axe4)=subplots(len(Cas),sharex='all',figsize=(8,6))
gcf().subplots_adjust(left = 0.1, bottom = 0.01,right = 1, top = 1, wspace = 0, hspace = 1)
axes=(axe1,axe2,axe3,axe4)

for i in range (len(Cas)):

    compteur = i

    h = Cas[i][0]
    tau = Cas[i][1]
    N = int(1/h)+1
    c = 2*tau/h

    AXE = axes[compteur]
    AXE.set(title="Explicite Centré \nh = "+str(h)+" tau= "+str(tau))

    for t in Lt:

        X,T,U = SolutionEC(h,N,tau,c)
        n=int(t/tau)
        AXE.plot(X,U[n,:],label='t='+str(t))

    axe1.legend(loc='upper right')

show()
close()

#Implicite centré

fig,(axe1,ax2,ax3,ax4)=subplots(len(Cas),sharex='all',figsize=(8,6))
gcf().subplots_adjust(left = 0.1, bottom = 0.01,right = 1, top = 1, wspace = 0, hspace = 1)
axes=(axe1,ax2,ax3,ax4)

for i in range (len(Cas)):

    compteur = i
    h = Cas[i][0]
    tau = Cas[i][1]
    N = int(1/h)+1
    c = 2*tau/h

    AXE = axes[compteur]
    AXE.set(title="Implicite Centré \nh = "+str(h)+" tau= "+str(tau))

    for t in Lt:
        X,T,U = SolutionIC(h,N,tau,c)
        n=int(t/tau)
        AXE.plot(X,U[n,:],label='t='+str(t))
    axe1.legend(loc='upper right')

show()
close()

#Explicite décentré Amont

fig,(axe1,ax2,ax3,ax4)=subplots(len(Cas),sharex='all',figsize=(8,6))
gcf().subplots_adjust(left = 0.1, bottom = 0.01,right = 1, top = 1, wspace = 0, hspace = 1)
axes=(axe1,ax2,ax3,ax4)

for i in range (len(Cas)):
    compteur = i
    h = Cas[i][0]
    tau = Cas[i][1]
    N = int(1/h)+1
    c = 2*tau/h

    AXE = axes[compteur]
    AXE.set(title="Ecplicite Décentré Amont \nh = "+str(h)+" tau= "+str(tau))
    
    for t in Lt:
        X,T,U = SolutionEDA(h,N,tau,c)
        n=int(t/tau)
        AXE.plot(X,U[n,:],label='t='+str(t))
    axe1.legend(loc='upper right')

show()
close()

##Lax-Friedrichs

fig,(axe1,ax2,ax3,ax4)=subplots(len(Cas),sharex='all',figsize=(8,6))
gcf().subplots_adjust(left = 0.1, bottom = 0.01,right = 1, top = 1, wspace = 0, hspace = 1)
axes=(axe1,ax2,ax3,ax4)

for i in range (len(Cas)):
    compteur = i
    h = Cas[i][0]
    tau = Cas[i][1]
    N = int(1/h)+1
    c = 2*tau/h

    AXE = axes[compteur]
    AXE.set(title="Lax-Friedrichs \nh = "+str(h)+" tau= "+str(tau))
    
    for t in Lt:
        X,T,U = SolutionLaxF(h,N,tau,c)
        n=int(t/tau)
        AXE.plot(X,U[n,:],label='t='+str(t))
    axe1.legend(loc='upper right')

show()
close()

##Lax-Wendroff

fig,(axe1,ax2,ax3,ax4)=subplots(len(Cas),sharex='all',figsize=(8,6))
gcf().subplots_adjust(left = 0.1, bottom = 0.01,right = 1, top = 1, wspace = 0, hspace = 1)
axes=(axe1,ax2,ax3,ax4)

for i in range (len(Cas)):
    compteur = i
    h = Cas[i][0]
    tau = Cas[i][1]
    N = int(1/h)+1
    c = 2*tau/h

    AXE = axes[compteur]
    AXE.set(title="Lax-Wendroff \nh= "+str(h)+" tau= "+str(tau))

    for t in Lt:
        X,T,U = SolutionLaxW(h,N,tau,c)
        n=int(t/tau)
        AXE.plot(X,U[n,:],label='t='+str(t))
    axe1.legend(loc='upper right')

show()
close()
