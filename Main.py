from render import Render
import numpy as np


from blackhole import KerrOrbitPair

#dataN = np.array([[0,0,0,0,1.25,-0.75,0,0,0,0,0,1],[0,100,0,0,2.125,0,1.875,0,0,1,0,0],[0,30,40,120,1,0,0,0,0,0,0,0]])
#N=3
#data0 = Nbody(dataN,N)
#print(data0)

print(KerrOrbitPair(np.array([[0,0,0,0,1,0,0,0,0,0,0,0],[0,100,0,0,1,0,0.1,0,0,0,0,0]])))





Render()

#print()