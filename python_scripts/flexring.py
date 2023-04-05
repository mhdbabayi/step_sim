import numpy as np

class Road:
    def __init__(self, step_width, step_height,step_profile_phase = np.pi, length = 5) -> None:
        self.length = length
        x1 = np.array([0, length/2- step_width/2])
        y1 = np.array([0,0])
        x_step = length/2 + np.linspace(-step_width/2 , step_width/2, np.int32(step_width/0.01))[1:]
        step_phase = np.linspace(0, step_profile_phase, len(x_step))
        y_step = step_height/2 - (step_height/2)*np.cos(step_phase)
        x2 = np.array([x_step[-1] + 0.01, self.length])
        y2 = np.array([y_step[-1], y_step[-1]])
        self.x = np.hstack((x1 , x_step , x2))
        self.y = np.hstack((y1, y_step , y2))


class Tyre:
    def __init__(self, initial_x, initial_y,road:Road,free_radius = 1, node_res_deg = 1) -> None:
        self.centre_x = initial_x
        self.centre_y = initial_y
        self.road = Road
        self.free_radius = free_radius
        self.node_zero = TyreNode(self,theta =0)
        theta = 0
        last_generated_node = self.node_zero
        while (theta<2*np.pi):
            theta = theta + np.deg2rad(node_res_deg)
            last_generated_node.next = TyreNode(theta=theta, previous_node=last_generated_node)
            last_generated_node = last_generated_node.next
        self.node_zero.prev =last_generated_node
        self.contact_inds = []

class TyreNode:
    def __init__(self,tyre, theta, next_node=None, previous_node =None):
        self.tyre:Tyre = tyre
        self.next:TyreNode = next_node
        self.prev: TyreNode = previous_node
        self.theta = theta
        self.deformation = 0
        self.x = self.tyre.centre_x + np.cos(self.theta)*self.tyre.free_radius
        self.y = self.tyre.centre_y + np.sign(self.theta)*self.tyre.free_radius

