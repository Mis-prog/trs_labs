import matplotlib.pyplot as ml
import numpy as np

x=np.linspace(1,2,100)
y=np.linspace(1,2,100)

ml.plot(x,y)
ml.show()

ml.savefig("test.png")