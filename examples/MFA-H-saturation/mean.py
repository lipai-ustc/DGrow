import numpy as np
a=np.loadtxt("nn")
a[:,2]=a[:,2]*2
a[:,3]=a[:,3]*2
print(a.mean(axis=0))
print(a.mean())
