import numpy as np
from enum import Enum
from euclid3 import Vector2
from dataclasses import dataclass
class ConstraintType(Enum):
    '''
    Name of constraint referes to the axis that's not updated 
    in the update method. The axis can still be explicitly changed 
    '''
    X = 1 << 0
    Y= 1 << 1
    THETA = 1 <<2
    

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
        self.states = RigidBody.State(
                                      position=Vector2(initial_x, initial_y),
                                      velocity=Vector2(initial_x_dot, initial_y_dot),
                                      acceleration= Vector2(0 , 0))
        self.forces = Vector2(initial_force_x , initial_force_y)
        self.mass = mass
        '''
        the constraint is a string with three digits, each representing a bit
        from left to right, Theta, Y, X
        for example, '001' means, object is constrained in the X direction, and theta and Y are updated in the dynamics 
        '''
        self.constraint = constraint_type
    def update_states(self, external_forces :Vector2 =Vector2(0 , 0)):
        self.update_forces(external_forces)
        self.update_accelerations()
        self.states.velocity = self.states.velocity + self.states.acceleration * self.simParameters["time_step"]
        self.states.position = self.states.position + self.states.velocity * self.simParameters["time_step"]
    def update_accelerations(self):
        if not(int(self.constraint) & ConstraintType.X.value):
            self.states.acceleration.x = self.forces.x / self.mass
        if not(int(self.constraint) & ConstraintType.Y.value):
            self.states.acceleration.y = self.forces.y/self.mass
    def update_forces(self, external_forces:Vector2 = Vector2(0 , 0)):
        self.forces = Vector2(0, 0)
        self.forces = self.forces + external_forces
    def iterate(self):
        self.update_forces()
        self.update_accelerations()
        self.update_states()
    @dataclass
    class State:
        position = Vector2(0. , 0.)
        velocity = Vector2(0., 0.)
        acceleration = Vector2(0., 0.)
        def __init__(self, position, velocity, acceleration):
            self.position = position
            self.velocity = velocity
            self.acceleration = acceleration
            


    


