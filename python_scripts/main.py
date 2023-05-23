import flexring as flx
import time
import os
import numpy as np
from matplotlib import pyplot as plt
os.system("clear")

# defining sim objects, all moving objects inherit from rigid body
forward_speed = 2.7

road = flx.Road(
                step_width=0.01,
                step_height=0.08,
                step_profile_phase=np.pi,
                high_res=True
                )
tyre = flx.Tyre_Continous(initial_x=2.2,
                initial_y=0.30,
                road=road,
                free_radius=0.35,
                x_speed=forward_speed,
                y_speed=0,
                )
q_car = flx.SprungMass(tyre_inst=tyre,
                       mass = 500,
                       speed_x=tyre.states.velocity.x,
                       speed_y=1,
                       spring_neutral_length=1,
                       damping_ratio=0.1
                       )
qmain_fig = plt.figure()


for i in range(500):    
    st = time.time() # For timing the main operations
    # main dynamics updates, should really be done in a function
    # but because of the way the external forces are implemented, it's done 
    # explicityly here
    
    q_car.update_states()
    tyre.update_states(-(q_car.spring_force + q_car.damper_force))
    print(f'{1000*(time.time() - st):.1f} ms/t {q_car.states.velocity.y:0.3f}')
    # draw results
    tyre.draw()
    q_car.draw()
    plt.plot(road.x, road.y, color="black")
    plt.gca().set_aspect('equal')
    plt.draw()
    plt.xlim((tyre.states.position.x - tyre.free_radius*1.1,
              tyre.states.position.x + tyre.free_radius*1.5))
    plt.ylim((tyre.states.position.y - tyre.free_radius*1.1, 
              q_car.states.position.y + tyre.free_radius*1.1))
    while not plt.waitforbuttonpress():
        pass
    plt.cla()
