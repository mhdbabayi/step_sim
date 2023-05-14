import flexring as flx
import time
import os
import numpy as np
from matplotlib import pyplot as plt
os.system("clear")

# defining sim objects, all moving objects inherit from rigid body
forward_speed = 3.

road = flx.Road(
                step_width=0.1,
                step_height=0.2,
                step_profile_phase=np.pi
                )
tyre = flx.Tyre(initial_x=2.1,
                initial_y=0.31,
                road=road,
                free_radius=0.35,
                x_speed=forward_speed
                )
q_car = flx.SprungMass(tyre_inst=tyre,
                       mass = 500,
                       speed_x=tyre.states.velocity.x,
                       speed_y=1,
                       spring_neutral_length=1,
                       damping_ratio=0.1
                       )
qmain_fig = plt.figure()


for i in range(150):    
    st = time.time() # for timing the main operations
    # main dynamics updates, should really be done in a function
    # but because of the way the external forces are implemented, it's done 
    # explicityly here
    
    tyre.update_states(-(q_car.damper_force + q_car.spring_force))
    q_car.update_states()

    print(f'{1000*(time.time() - st):.1f} ms/t {q_car.states.velocity.y:0.3f}')

    # draw results
    tyre.draw()
    q_car.draw()
    plt.plot(road.x, road.y,'.-')
    plt.gca().set_aspect('equal')
    plt.draw()
    if len(tyre.contacts) > 0:
        c = tyre.contacts[0]
        n = c.aft_penetration_node
        plt.plot((tyre.states.position.x, c.centre_node.x),
                 (tyre.states.position.y , c.centre_node.y))
    plt.xlim((tyre.states.position.x - tyre.free_radius*1.1,
              tyre.states.position.x + tyre.free_radius*1.5))
    plt.ylim((tyre.states.position.y - tyre.free_radius*1.1, 
              q_car.states.position.y + tyre.free_radius*1.1))
    while not plt.waitforbuttonpress():
        pass
    plt.cla()
