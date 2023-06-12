import flexring as flx
import time
import os
import numpy as np
from matplotlib import pyplot as plt
os.system("clear")
from scipy import io
# constants and initial values
forward_speed = 1.9*12./7.
tyre_radius = 0.788/2
initial_x = 2.5 - tyre_radius
sprung_mass = 700
unsprung_mass = 50
initial_y = tyre_radius - (sprung_mass + unsprung_mass)*10/(flx.Tyre_Continous.lump_stiffness)
# defining sim objects, all moving objects inherit from rigid body

road = flx.Road(
                step_width=0.1,
                step_height=0.2,
                step_profile_phase=np.pi,
                high_res=True
                )
tyre = flx.Tyre_Continous(initial_x=2.2,
                          initial_y=initial_y+0.03,
                          boundary_condition_file="./step_sim/beta_5.mat",
                          mass=unsprung_mass,
                          road=road,
                          free_radius=tyre_radius,
                          x_speed=forward_speed,
                          y_speed=0,
                        )
q_car = flx.SprungMass(tyre_inst=tyre,
                       mass = sprung_mass,
                       speed_x=forward_speed,
                       speed_y=0,
                       spring_neutral_length=1,
                       natural_frequency_hz=2,
                       damping_ratio=0.2
                       )
qmain_fig, Ax = plt.subplots(2 , 1)
tyre.find_new_contacts()

#for i in range(500): 
logged_data = []   
step = 0
while tyre.states.position.x < 4:
    step += 1
    st = time.time() # For timing the main operations
    # main dynamics updates, should really be done in a function
    # but because of the way the external forces are implemented, it's done 
    # explicityly here
    plt.sca(Ax[0])
    q_car.update_states()
    tyre.update_states(-(q_car.spring_force + q_car.damper_force))
    if np.mod(step , 50) == 0:
        print(f'{1000*(time.time() - st)/50:.1f} ms/t {q_car.states.velocity.y:0.3f}')
    # draw results
    for ax in Ax:
        ax.cla()
    # Ax.cla()
    tyre.draw()
    q_car.draw()
    plt.plot(road.x, road.y, color="brown")
    plt.gca().set_aspect('equal')
    plt.draw()
    # plt.xlim((tyre.states.position.x - tyre.free_radius*1.1,
    #           tyre.states.position.x + tyre.free_radius*1.5))
    # plt.ylim((tyre.states.position.y - tyre.free_radius*1.1, 
    #           q_car.states.position.y + tyre.free_radius*1.1))
    plt.xlim((tyre.contacts[-1].centre_point.x - tyre.free_radius, tyre.contacts[-1].centre_point.x + tyre.free_radius))
    plt.ylim((tyre.contacts[-1].centre_point.y - tyre.free_radius, tyre.contacts[-1].centre_point.y + tyre.free_radius))

    plt.sca(Ax[1])
    for c in tyre.contacts:
        c.draw_pressure()
    plt.pause(0.1)
    if step > 10:
        while not plt.waitforbuttonpress():
            pass
    #if np.mod(step , 50) == 0:
    #    plt.pause(0.01)
    logged_data.append([tyre.forces, tyre.external_forces, tyre.states.position])
for ax in Ax:
    ax.cla()
print("finished")
logged_forces = np.array([v[0] for v in logged_data])
ex_forces = np.array([np.double(v[1]) for v in logged_data])
position = np.array([v[2] for v in logged_data])
io.savemat("/Users/mahdibabayi/beam_tyre/step_out.mat",
        {"position":position, "force": logged_forces,"suspension_forces":ex_forces})

plt.plot(position[: , 0], logged_forces[: , 1])
plt.show()
