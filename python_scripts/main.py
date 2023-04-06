import flexring
import os
import numpy as np
os.system("clear")
from matplotlib import pyplot as plt
road = flexring.Road(step_width=0.1, step_height=0.2, )
tyre = flexring.Tyre(initial_x=2, initial_y=0.8,road=road)
plt.plot(road.x, road.y)
plt.plot(tyre.centre_x , tyre.centre_y, 'r*')
tyre.update_penetrations()
n = tyre.node_zero.next
while n is not tyre.node_zero:
    plt.plot(n.x , n.y, 'r.')
    if n.penetration_point:
        plt.plot(n.penetration_point[0], n.penetration_point[1] , 'm.')
        plt.plot(np.array([n.x, n.penetration_point[0]]),
                 np.array([n.y, n.penetration_point[1]]))
    n = n.next
plt.gca().set_aspect('equal')
plt.show()