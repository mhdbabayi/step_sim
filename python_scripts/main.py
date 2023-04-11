import flexring
import os
import numpy as np
os.system("clear")
from matplotlib import pyplot as plt
road = flexring.Road(step_width=0.1, step_height=0.05,step_profile_phase=np.pi*2)
tyre = flexring.Tyre(initial_x=2.4, initial_y=0.36,road=road, free_radius=0.35, node_res_deg=0.5)
speed_y = 0
speed_x = 1
main_fig = plt.figure()
y0 = 0.36
for i in range(10):    
    '''
    tyre = flexring.Tyre(initial_x=x0, initial_y=y0,road=road, free_radius=0.35)
    tyre.update_penetrations()
    tyre.update_contacts()
    tyre.update_deformation()
    '''
    tyre.update_state(speed_x=speed_x, speed_y=speed_y)
    tyre.draw()
    plt.plot(road.x, road.y,'.-')
    plt.gca().set_aspect('equal')
    plt.draw()
    c = tyre.contacts[0]
    n = c.aft_penetration_node
    plt.plot((tyre.centre_x, c.centre_node.x),
             (tyre.centre_y , c.centre_node.y))
    '''
    while (n:=n.next) is not c.fore_penetration_node:
        plt.plot((n.penetration_point[0], n.x), (n.penetration_point[1], n.y))
    '''
    if len(tyre.contacts) > 0:
        plt.xlim((tyre.contacts[0].centre_node.x - 0.1,
               tyre.contacts[0].centre_node.x + 0.1))
        plt.ylim((-0.05, tyre.centre_y/2))
    while not plt.waitforbuttonpress():
        pass
    plt.cla()
