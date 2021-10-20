#!/usr/bin/env python
# coding: utf-8

# # The Rice-Mele (RM) model
# 
# The Rice-Mele (RM) model adds another term to the SSH model.

# ---------------
# ## The RM model on open chain
# These are the same results as for the open SSH model chain as seen before, except we choose a non-zero value of $\epsilon$ to study the RM model.

# In[1]:


op=ssh_open(6,1,1.6,0.2)
op.ham_print()


# In[ ]:


L=28
t1=1
t2=1.6
eps=0.2
index=int(np.floor(L/2)-1)
op=ssh_open(L,t1,t2,eps)
plot=1
op.spectrum(L, t1, t2,plot);
#op.wfs(index,L,t1,t2)

