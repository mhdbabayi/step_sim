import flexring
import os
import numpy as np
os.system("clear")
from matplotlib import pyplot as plt
road = flexring.Road(step_width=0.3, step_height=0.1,step_profile_phase=np.pi)



main_fig = plt.figure()
y0 = 1.01
for x0 in np.linspace(2 , 3, 20):    
    tyre = flexring.Tyre(initial_x=x0, initial_y=y0,road=road)
    plt.plot(road.x, road.y,'.-')
    tyre.update_penetrations()
    tyre.update_contacts()
    tyre.update_deformation()
    tyre.draw()
    plt.gca().set_aspect('equal')
    plt.draw()
    while not plt.waitforbuttonpress():
        pass
    plt.cla()
