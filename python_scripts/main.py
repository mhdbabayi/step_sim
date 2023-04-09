import flexring
import os
import numpy as np
os.system("clear")
from matplotlib import pyplot as plt
road = flexring.Road(step_width=0.1, step_height=0.5, )
tyre = flexring.Tyre(initial_x=1.7, initial_y=0.95,road=road)
main_fig = plt.figure()
plt.plot(road.x, road.y)
plt.plot(tyre.centre_x , tyre.centre_y, 'r*')
tyre.update_penetrations()
tyre.update_contacts()
n = tyre.node_zero.next
while n is not tyre.node_zero:
    plt.plot(n.x , n.y, 'r.')
    if n.penetration_point:
        plt.plot(n.penetration_point[0], n.penetration_point[1] , 'm.')
    n = n.next
plt.plot(tyre.node_zero.x , tyre.node_zero.y , marker='o', color = 'black')
[plt.plot(c.centre_node.penetration_point[0], 
          c.centre_node.penetration_point[1],'o') for c in tyre.contacts]
[plt.plot(c.fore_separation_node.penetration_point[0],
          c.fore_separation_node.penetration_point[1], '*', color='black', markersize=10) for c in tyre.contacts]
[plt.plot(c.aft_separation_node.penetration_point[0],
          c.aft_separation_node.penetration_point[1], '*', color='green', markersize = 10) for c in tyre.contacts]
plt.gca().set_aspect('equal')
'''second_fig = plt.figure()
n = tyre.node_zero.next
while n is not tyre.node_zero:
    if n.road_dr_dtheta is not None:
        plt.plot(n.theta,n.road_dr_dtheta, '.')
    n = n.next
'''
c = tyre.contacts[-1]
n = c.centre_node
while n is not c.aft_penetration_node:
    print((f'{-2*tyre.beta*n.road_dr_dtheta:.5f},\t'
           f'{n.road_ddr_dtheta:.5f}'))
    if n is c.fore_separation_node:
        print('here')
    n = n.prev
plt.show()