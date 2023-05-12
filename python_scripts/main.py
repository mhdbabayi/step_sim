import flexring as flx
import time
import os
import numpy as np
from matplotlib import pyplot as plt
os.system("clear")

road = flx.Road(step_width=0.05, step_height=0.1,step_profile_phase=np.pi)
tyre = flx.Tyre(initial_x=2.1, initial_y=0.31,road=road, free_radius=0.35,
                     x_speed=1)
q_car = flx.SprungMass(tyre_inst=tyre,mass = 500, speed_x=1,speed_y=1, spring_neutral_length=1)
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
    tyre.update_states(-(q_car.damper_force + q_car.spring_force))
    q_car.update_states()
    print(f'{1000*(time.time() - st):.1f} ms/t {q_car.states["y_dot"]:0.3f}')
    tyre.draw()
    q_car.draw()
    plt.plot(road.x, road.y,'.-')
    plt.gca().set_aspect('equal')
    plt.draw()
    if len(tyre.contacts) > 0:
        c = tyre.contacts[0]
        n = c.aft_penetration_node
        plt.plot((tyre.states['x'], c.centre_node.x),
                 (tyre.states['x'] , c.centre_node.y))
        plt.xlim((tyre.states['x'] - tyre.free_radius*1.2
                  ,tyre.states['x'] + tyre.free_radius*1.2))
        plt.ylim((tyre.states['y'] - tyre.free_radius*1.1
                  ,q_car.states['y'] + tyre.free_radius*1.1))
        '''
        while (n:=n.next) is not c.fore_penetration_node:
            plt.plot((n.penetration_point[0], n.x), (n.penetration_point[1], n.y))
        '''
    while not plt.waitforbuttonpress():
        pass
    plt.cla()
