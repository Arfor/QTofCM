#!/usr/bin/env python
# coding: utf-8

# # Exercises
# ----------
# The exercises below are supposed to be done *before* class and will be discussed *in* class. To learn most from the exercises you are encouraged to try, fail and succeed as many of them as you can.
# 
# ## SSH Bulk & Topology
# 1. The Born-von Karman implementation of PBC says that translation by $Na$ has to be the identity operator, hence $e^{ik(Na)}= 1$. Find the allowed momenta in the first Brillouin Zone (BZ). How many are there?
# 
# ```{tip} 
# The first Brillouin Zone is where $k = ]-\pi,\pi]
# ```
# 2. Rewrite $H^{PBC}$ in the momentum basis using the momentum states
# 
# $$
# \ket{k\alpha} \equiv \frac{1}{\sqrt{N}} \sum^N_{i=1}e^{ikna}\ket{n,\alpha}, \ket{n,\alpha} \equiv \frac{1}{\sqrt{N}} \sum_{k\in BZ}e^{-ikna}\ket{k\alpha},
# $$
# and identity $\sum^N_{n=1}e^{i2\pi mn/N}=N\delta_{m,0}$ for $m=0,...,N-1$.
# 
# 3. Plot the path of $f_k$ in the complex plane as one varies $k \in BZ$, contrasting the cases $t'>t$ and $t'<t$. What happens to the path and to th bandstructure when $t=t'$? Identify the origin of the complex plane as the key point and contrast the two cases in terms of a *winding number* $W(t,t')$ around the origin.
# 
# ```{admonition} Spoiler!
# :class: dropdown
# TD: Insert plot of the path of $f_k$
# ```
# 
# 4. A way to characterize the bandstructure is through the *Berry phase* $\gamma$ of the band eigenstates:
# 
# $$
# \gamma_{\pm} = i \int^{2\pi /a}_{0}dk<\psi_{k\pm}|\partial_k|\psi_{k\pm}>
# $$ {berryPhase}
# 
# What is the relationship of $W(t,t')$ and $\gamma_{\pm}(t,t')$?
# 
# ## SSH Open Chain & Edge States
# 
# ## Rice Mele Chain: Symmetry & Edge Modes
# 
# ## Numerical & Graphical Exploration
# 
# ## Generalization of the SSH Chain
