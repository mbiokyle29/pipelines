
# coding: utf-8

# In[43]:

import matplotlib.pyplot as plt
import numpy

fig = plt.figure()

depths = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
counts = [107138,48086,12866,3947,1260,542,121,87,192,171,100,32,16,16,20,19]
summed_counts = []

for i in range(0,len(counts)):
    el = counts.pop(0)
    el += sum(counts)
    summed_counts.append(el)

fig.xlim([0, max(depths)+1])
fig.plot(depths, summed_counts)
fig.fill_between(depths, summed_counts)
fig.grid(True)
fig.xlabel("Depth")
fig.ylabel("BP count")
fig.title("Depth of Coverage for sample #1")
fig.show()

