import numpy as np
import flexring as flx

class DynamicObject:
    simParameters = {
    "universalTime" : 0,
    "timeStep" : 0.01}
    def initialize(self):
        pass
    def iterate(self):
        pass

class QuarterCar(DynamicObject):
    def __init__(self, body_mass = 500,
                tyre_mass = 50,
                natural_frequency_hrz = 1.5,
                damping_ratio = 0.3,
                initial_x = 0,
                spring_netural_length = 0.3,
                initial_tyre_y = 0) -> None:
        self.body_mass = body_mass
        self.tyre_mass = tyre_mass
        self.natural_frequency_rad_s = np.deg2rad(360*natural_frequency_hrz)
        self.damping_ratio = damping_ratio
        self.spring_stiffness = self.natural_frequency_rad_s**2 * self.body_mass
        self.damping_coefficient = 2*body_mass*self.natural_frequency_rad_s*damping_ratio
        self.x = initial_x
        self.tyre_y = initial_tyre_y
        self.tyre_y_dot = 0
        self.tyre_y_dot_dot = 0
        self.spring_neutral_length = spring_netural_length
        self.body_y  = initial_tyre_y + spring_netural_length
        self.body_y_dot = 0
        self.body_y_dot_dot = 0
        self.spring_force = 0
        self.damper_force = 0

    def set_forces(self):
        self.spring_force = self.spring_stiffness*(
            self.spring_neutral_length - (self.body_y - self.tyre_y))
        self.damper_force = self.damping_coefficient*(
            self.tyre_y_dot - self.body_y_dot
        )
    def update_states(self):
        self.body_y_dot = self.body_y_dot_dot + 
    def iterate(self):
        self.set_force
        


