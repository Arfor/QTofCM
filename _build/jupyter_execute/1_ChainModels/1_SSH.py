#!/usr/bin/env python
# coding: utf-8

# # The Su-Schrieffer-Heeger (SSH) model
# 
# The Su-Schrieffer-Heeger (SSH) model was constructed to describe the electronic properties of polymer chains. As a one-dimensional system it of- fers the simplest demonstration of concepts of edge states and topological electronic states, which are a central topic of condensed-matter physics in the past decades.
# 
# Both the SSH model and the Rice-Mele (RM) model are defined on the same chain of lattice sites along which an electron can hop. In the SSH model there are exactly two types of hopping whose energies are $t_{1}$, $t_2$ and they alternate along the chain, producing a unit-cell with two atoms (A and B):
# 
# ```{image} ../images/SSHchain.png
# :alt: SSH Chain
# :width: 500px
# :align: center
# ```
# The RM model adds an on-site potential +(-)$\epsilon$ on all A(B) sites. Combined, for an open chain of $N$ unit-cells the tight-binding Hamiltonian used in this notebook is in the form:
# 
# $$
# H_N=\sum_{i=1}^{2N} t_{i,i+1}\left(|\phi_{i+1}\rangle\langle\phi_i|+|\phi_i\rangle\langle\phi_{i+1}|\right)-(-1)^i \epsilon|\phi_{i}\rangle\langle\phi_i|,\;\;\;\;t_{i,i+1}\equiv\begin{cases}t_1,&i\in\textrm{odd}\\t_2,&i\in\textrm{even}\\\end{cases},
# $$
# 
# where $\epsilon\equiv0$ in the SSH model, while $|\phi_i\rangle$ are an orthonormalized basis of states of an electron localized on site $i$.
# 
# For an infinite chain the bandstructure is: 
# 
# $$
# E_\pm(k)=\pm\sqrt{t_1^2+t_2^2+2t_1 t_2 cos(k)},
# $$
# 
# with the first Brillouin zone given by $k\in[-\pi,\pi)$, and we set the lattice constant (i.e., the $A-A$ distance) $a\equiv1$ throughout.
# 
# For a finite chain with periodic boundary conditions (PBC), the discrete spectrum is given by $E_\pm(k_m)$, where $k_m=2\pi\frac{m}{N}$ with $m=0,1\ldots N-1$ (which is equivalent to $k$ values in the first BZ).
# 

# ------------
# ## Import python packages, define constants

# In[1]:


import numpy as np
from numpy import linalg as LA
# import sympy as sp
# from sympy.functions import exp
from matplotlib import pyplot as plt
# from __future__ import print_function
# from ipywidgets import interact, interactive, fixed, interact_manual
# import ipywidgets as widgets
import holoviews as hv;
hv.notebook_extension('matplotlib', logo=False)
from holoviews import opts;

Pi=np.pi


# ## Define and solve SSH models with PBC
# ### Infinite chain with PBC: Bands
# 
# **ssh_PBC_bands($t_1$,$t_2$,$\epsilon$)** is a class containing the continuous bands of an infinite chain.
# 
# **plot_bands($k_{min}$,$k_{max}$,$T_1$,$T_2$)** function plots two panels: (left) The bandstructure over the range $k\in[k_{min},k_{max}]$, after imposing the new values $t_1\equiv T_1$, $t_2\equiv T_2$; (right) The trace of $f^*_k$ in the complex plane (blue dot is the origin), over the range of $k$.

# ````{margin} 
# ```{admonition} Code Blocks
# :class: tip
# We will do our best to make the code as readable as possible, so try to understand every step which encodes some physics! Bonus points if you find ways to make the code faster or more readable ;)
# ``` 
# ````

# In[2]:


dims = dict(kx = 'k_x',ky = 'k_y', mu = 'µ', delta = 'Δ', t = 't', E='ε',mu_t = 'µ/t',E_d='ε/Δ', wf='ψ',
            azi_winding='q', radi_winding= 'p', radius = 'r', R0 = 'R_0', j='j', lw='lw', dens='$|ψ|^2$',  B='B', x='x')

piTicks = [(-Pi, '-π'),(0, '0'),(Pi, 'π')]

class ssh_PBC_bands:
    def __init__(self,t1,t2,eps=0):
        self.t1=t1
        self.t2=t2
        self.eps=eps
        self.nrofbands=2 #?
    
    def bands(self,k):
        return [-np.sqrt(self.eps**2+self.t1**2+self.t2**2+2*self.t1*self.t2*np.cos(k)), 
                np.sqrt(self.eps**2+self.t1**2+self.t2**2+2*self.t1*self.t2*np.cos(k))]

    def H12(self,k):
        return self.t1+self.t2*np.exp(1j*k)

    def plot_bands(self,kmin,kmax,newt1,newt2):
        self.t1=newt1
        self.t2=newt2
        ks = np.linspace(kmin, kmax,51)
        spectrum = self.bands(ks)/np.abs(self.t1)
        
        # plot spectrum
        plotOpts = dict(xticks=piTicks, color='black',xlim=(-Pi,Pi),ylim=(-Pi,Pi)) #define options for your plot
        spectrumPlot = hv.Curve((ks,spectrum[0]), kdims=['k'],vdims=['$E/|t_1|$']).opts(**plotOpts) #plot first band and apply options
        spectrumPlot *= hv.Curve((ks, spectrum[1])).opts(color='black') # Add second band --> (*) creates an Overlay in Holoviews
        spectrumPlot.opts(title='Bands of infinite SSH Chain ($t_1={}, t_2={}$)'.format(self.t1,self.t2)) # Add title
        
        # plot winding
        plotOpts = dict(xticks=piTicks, yticks=piTicks, xlim=(-Pi,Pi), ylim=(-Pi,Pi), color='blue')
        windingPlot = hv.Curve((np.real(self.H12(ks)),np.imag(self.H12(ks))), kdims=['$Re[H_{12}(k)]$'], vdims=['$Im[H_{12}(k)]$']).opts(**plotOpts)
        windingPlot *= hv.Scatter((0,0)).opts(color='blue')
        windingPlot.opts(title='Winding ($t_1={}, t_2={}$)'.format(self.t1,self.t2))
                         
        return (spectrumPlot+windingPlot) #--> (+) creates a Layout in Holoviews 


# In[3]:


eps=0
(ssh_PBC_bands(1,0.1,eps).plot_bands(-Pi,Pi,1,0.1)+ssh_PBC_bands(1,1.1,eps).plot_bands(-Pi,Pi,1,1.1)).cols(2)


# ### Interactive bands
# A slider to change value of $t_2$ in real-time, with fixed values $t_1,\epsilon$.

# In[4]:


t1=1
modelSSH=ssh_PBC_bands(t1,0,0)

def t2PlotBand(t2):
    return modelSSH.plot_bands(-Pi,Pi,t1,t2)

t2Values = np.arange(-2,2,0.25)
spectrumMap = hv.HoloMap({t2: t2PlotBand(t2) for t2 in t2Values}, kdims = ['t2']).collate()
spectrumMap


# -------------------
