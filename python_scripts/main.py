import flexring
import time
import os
import numpy as np
from matplotlib import pyplot as plt
os.system("clear")

road = flexring.Road(step_width=0.05, step_height=0.1,step_profile_phase=np.pi)
tyre = flexring.Tyre(initial_x=2, initial_y=0.3,road=road, free_radius=0.35,
                     x_speed=1)
speed_y = 0
speed_x = 1
main_fig = plt.figure()
y0 = 0.36
for i in range(50):    
    '''
    tyre = flexring.Tyre(initial_x=x0, initial_y=y0,road=road, free_radius=0.35)
    tyre.update_penetrations()
    tyre.update_contacts()
    tyre.update_deformation()
    '''
    st = time.time()
    tyre.update_states()
    print(f'{1000*(time.time() - st):.1f} ms')
    tyre.draw()
    plt.plot(road.x, road.y,'.-')
    plt.gca().set_aspect('equal')
    plt.draw()
    if len(tyre.contacts) > 0:
        c = tyre.contacts[0]
        n = c.aft_penetration_node
        plt.plot((tyre.states['x'], c.centre_node.x),
                 (tyre.states['x'] , c.centre_node.y))
        plt.xlim((c.centre_node.x - tyre.free_radius
                  ,c.centre_node.x + tyre.free_radius))
        plt.ylim((c.centre_node.y - tyre.free_radius
                  ,c.centre_node.y + tyre.free_radius))
        '''
        while (n:=n.next) is not c.fore_penetration_node:
            plt.plot((n.penetration_point[0], n.x), (n.penetration_point[1], n.y))
        '''
    while not plt.waitforbuttonpress():
        pass
    plt.cla()
