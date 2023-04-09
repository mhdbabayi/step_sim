import flexring
import os
import numpy as np
os.system("clear")
from matplotlib import pyplot as plt
road = flexring.Road(step_width=0.1, step_height=0.5, )
tyre = flexring.Tyre(initial_x=1, initial_y=0.85,road=road)
main_fig = plt.figure()
plt.plot(road.x, road.y)

tyre.update_penetrations()
tyre.update_contacts()
tyre.update_deformation()
tyre.draw()
plt.gca().set_aspect('equal')
if len(tyre.contacts) >0: 
    c = tyre.contacts[0]
    n = c.centre_node
    while n is not c.fore_penetration_node:
        print((f'{-2*(tyre.beta**2)*(n.road_dr_dtheta/c.tyre.beta + n.road_dr):.5f},\t'
            f'{n.road_ddr_dtheta:.5f}'))
        if n is c.aft_separation_node:
            print('here')
        n = n.next
plt.show()