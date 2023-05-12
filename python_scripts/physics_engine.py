import numpy as np
from enum import Enum
class ConstraintType(Enum):
    '''
    Name of constraint referes to the axis that's not updated 
    in the update method. The axis can still be explicitly changed 
    '''
    X = 1 << 0
    Y= 1 << 1
    THETA = 1 <<3
    

class DynamicObject:
    simParameters = {
    "universal_time" : 0,
    "time_step" : 0.01}
    def initialize(self):
        pass
    def iterate(self):
        pass
class RigidBody(DynamicObject):
    def __init__(self,
                mass,
                initial_x,
                initial_y,
                initial_x_dot,
                initial_y_dot,
                constraint_type:str,
                initial_force_y = 0,
                initial_force_x = 0) -> None:
        self.states = {'x': initial_x, 'y': initial_y,
                        'x_dot':initial_x_dot,'y_dot':initial_y_dot,
                        'x_dot_dot':0, 'y_dot_dot': 0}
        self.forces = {'x': initial_force_x, 'y': initial_force_y, 'theta': 0}
        self.mass = mass
        '''
        the constraint is a string with three digits, each representing a bit
        from left to right, Theta, Y, X
        for example, '001' means, object is constrained in the X direction, and theta and Y are updated in the dynamics 
        '''
        self.constraint = constraint_type
    def update_states(self, external_forces =0):
        self.update_forces(external_forces)
        self.update_accelerations()
        self.states['x_dot']  = self.states['x_dot'] + self.states['x_dot_dot'] * self.simParameters['time_step']
        self.states['x']  = self.states['x'] + self.states['x_dot'] * self.simParameters['time_step']

        self.states['y_dot']  = self.states['y_dot'] + self.states['y_dot_dot'] * self.simParameters['time_step']
        self.states['y']  = self.states['y'] + self.states['y_dot'] * self.simParameters['time_step']
    def update_accelerations(self):
        if not(int(self.constraint) & ConstraintType.X.value):
            self.states['x_dot_dot'] = self.forces['x'] / self.mass
        if not(int(self.constraint) & ConstraintType.Y.value):
            self.states['y_dot_dot'] = self.forces['y']/self.mass
    def update_forces(self, external_forces):
        raise NotImplementedError
    def iterate(self):
        self.update_forces()
        self.update_accelerations()
        self.update_states()



    


