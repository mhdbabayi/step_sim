import flexring
import os
import numpy as np
os.system("clear")
from matplotlib import pyplot as plt
road = flexring.Road(step_width=0.5, step_height=0.1, step_profile_phase=2*np.pi)
plt.plot(road.x, road.y)
plt.show()