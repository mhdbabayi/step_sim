import flexring as flx
import time
import os
import numpy as np
from matplotlib import pyplot as plt
os.system("clear")
from scipy import io

# defining sim objects, all moving objects inherit from rigid body
forward_speed = 10

road = flx.Road(
                step_width=0.1,
                step_height=0.15,
                step_profile_phase=np.pi,
                high_res=True
                )
tyre = flx.Tyre_Continous(initial_x=1.5,
                initial_y=0.34,
                road=road,
                free_radius=0.788/2,
                x_speed=forward_speed,
                y_speed=0,
                )
q_car = flx.SprungMass(tyre_inst=tyre,
                       mass = 630,
                       speed_x=tyre.states.velocity.x,
                       speed_y=0,
                       spring_neutral_length=1,
                       damping_ratio=0.1
                       )
qmain_fig, Ax = plt.subplots(2 , 1)
tyre.find_new_contacts()

#for i in range(500): 
logged_data = []   
while tyre.states.position.x < 3:
    st = time.time() # For timing the main operations
    # main dynamics updates, should really be done in a function
    # but because of the way the external forces are implemented, it's done 
    # explicityly here
    plt.sca(Ax[0])
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
    plt.sca(Ax[1])
    for c in tyre.contacts:
        c.draw_pressure()
    while not plt.waitforbuttonpress():
        pass
    logged_data.append([tyre.forces, tyre.states.position])
    for ax in Ax:
        ax.cla()
print("finished")
# logged_forces = np.array([v[0] for v in logged_data])
# position = np.array([v[1] for v in logged_data])
# io.savemat("/Users/mahdibabayi/beam_tyre/step_out.mat",
#         {"position":position, "force": logged_forces})
# plt.plot(position[: , 0], logged_forces[: , 1])
