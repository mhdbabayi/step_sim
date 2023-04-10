import flexring
import os
import numpy as np
os.system("clear")
from matplotlib import pyplot as plt
road = flexring.Road(step_width=0.3, step_height=0.5,step_profile_phase=np.pi*2)
main_fig = plt.figure()
x0 = 2.5
for y0 in np.linspace(1.5 , 0, 10):    
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
