#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt


# In[57]:
tk=300.0

fig = plt.figure()
# syntax for 3-D projection
ax = plt.axes(projection ='3d')
ax.set_xlim3d([0.0, 1000.0])
ax.set_xlabel('speed-x')

ax.set_ylim3d([50.0, 300.0])
ax.set_ylabel('speed-y')

ax.set_zlabel('MB')
ax.set_title("Maxwell-Distribution at T = "+str(tk))

# In[63]:


k = 1.38e-23
u = 1.6605e-27
m = 15.966*u
n =1000

#set temp range from 20K at least for 3-d plot 
vx= np.linspace(0,1000.0,n)
vy= np.linspace(0,1000,n)


# In[64]:


def MB(vx,vy,T):
    v = (vx**2+vy**2)**0.5	
    f1 = 4*np.pi*(m/(2*np.pi*k*T))**(3.0/2.0)
    f2 = v**2
    f3 = (np.exp(-((m*v**2)/(2*k*T))))
    return f1*f2*f3


# for i in range(0,10):
#     print(MB(i,373.0))

# In[65]:


#plt.plot(v,MB(v,T[100]))


# In[66]:

X,Y = np.meshgrid(vx,vy)
Z = MB(X,Y,tk) 
ax.plot_surface(X,Y,Z,cmap='viridis')
plt.show()

# In[ ]:


