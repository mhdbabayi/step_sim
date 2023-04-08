import numpy as np


def intersection(p1, p2, P1, P2):
    x0, y0 = p1
    x1, y1 = p2
    X0, Y0 = P1
    X1, Y1 = P2

    # calculate the denominator of the two equations
    den = (Y1 - Y0)*(x1 - x0) - (X1 - X0)*(y1-y0)

    # check if the lines are parallel
    if den == 0:
        return None

    # calculate the numerators of the two equations
    num1 = (X1 - X0)*(y0 - Y0) + (X0 - x0)*(Y1-Y0)
    num2 = (x1 - x0)*(y0 -Y0) + (X0 - x0)*(y1 - y0)

    # calculate the parameters t and u
    t = num1 / den
    u = num2 / den

    # check if the intersection point is within the line segments
    if t >= 0 and t <= 1 and u >= 0 and u <= 1:
        x_intersect = x0 + t*(x1 - x0)
        y_intersect = y0 + t*(y1 - y0)
        return (x_intersect, y_intersect)

    # no intersection found
    return None
def polar_derivative(X , Y , DY):
    #given point x,y and the deriavtive in cartesian coordiantes dy
    #calculates dr/d(theta) in polar coordinates
    x = np.asarray(X)
    y = np.asarray(Y)
    dy = np.asarray(DY)
    r = np.sqrt(x**2 + y**2)
    dr = r*(x + y*dy)/(x*dy - y)
    if np.isscalar(X):
        return dr.item()
    return dr
def polar_second_derivative(X , Y, DY , DDY):
    x = np.asarray(X)
    y = np.asarray(Y)
    dy = np.asarray(DY)
    ddy = np.asarray(DDY)
    r = np.sqrt(x**2 + y**2)
    dr_dtheta = polar_derivative(x , y , dy)
    dx_dtheta = dr_dtheta*x/r  - y
    dy_dtheta = dr_dtheta*y/r + x
    dydx_dtheta = ddy*dx_dtheta
    ddr_dtheta = ((x*dy - y)*(dx_dtheta + y*dydx_dtheta + dy*dy_dtheta) -\
                  (y*dy + x)*(dy*dx_dtheta + dydx_dtheta*x - dy_dtheta)) * r/(x*dy - y)**2\
                   +dr_dtheta*(y*dy  +x)/(x*dy - y)
    if np.isscalar(X):
        return ddr_dtheta.item()
    return ddr_dtheta


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
    beta = 1
    def __init__(self, initial_x, initial_y,road:Road,free_radius = 1., node_res_deg = 1.) -> None:
        self.centre_x = initial_x
        self.centre_y = initial_y
        self.road = road
        self.free_radius = free_radius
        self.node_zero = Tyre.node(self,theta =0)
        theta = 0
        last_generated_node = self.node_zero
        while (theta<2*np.pi):
            theta = theta + np.deg2rad(node_res_deg)
            last_generated_node.next = Tyre.node(tyre=self, theta=theta, previous_node=last_generated_node)
            last_generated_node = last_generated_node.next
        self.node_zero.prev =last_generated_node
        last_generated_node.next = self.node_zero
        self.contacts =[]
    def update_penetrations(self):
        current_node = self.node_zero.next
        while current_node is not self.node_zero:
            #road_point_idx = np.where((self.road.x > self.centre_x - self.free_radius) & 
            #                         (self.road.x < self.centre_x + self.free_radius))[0]
            current_node.penetration_point=None
            for i  in range(len(self.road.x[0:-1])):
                intersection_point =intersection((self.centre_x, self.centre_y),
                                                 (current_node.x, current_node.y),
                                                 (self.road.x[i], self.road.y[i]),
                                                 (self.road.x[i+1], self.road.y[i+1]))
                if intersection_point:
                    current_node.penetration_point = intersection_point
                    break
            current_node = current_node.next
        self.update_derivatives()
    def update_derivatives(self):
        current = self.node_zero.next
    # traverse the circular linked list
        while current is not self.node_zero:
            # check if the current node has two non-None neighbors on each side
            if current.prev.penetration_point is not None and current.next.penetration_point is not None:
                # calculate the first and second derivatives of y with respect to x
                h0 = current.penetration_point[0] - current.prev.penetration_point[0]
                h1 = current.next.penetration_point[0] - current.penetration_point[0]
                y0 = current.prev.penetration_point[1]
                y1 = current.penetration_point[1]
                y2 = current.next.penetration_point[1]
                current.dy = (y2 - y0)/(0.5*(h1+h1))
                current.ddy = (y2 + y0 - 2*y1)/(h1*h0)
                current.dr_dtheta = polar_derivative(X=current.x - self.centre_x,
                                                     Y = current.y - self.centre_y,
                                                     DY = current.dy)
                current.ddr_dtheta = polar_second_derivative(X = current.x - self.centre_x,
                                                             Y = current.y -self.centre_y,
                                                             DY = current.dy, DDY= current.ddy)
                # store the derivatives in the current node
            current = current.next
    def update_contacts(self):
        current_node = self.node_zero.next
        self.contacts = []
        while current_node is not self.node_zero:
            if current_node.penetration_point:
                self.contacts.append(Tyre.contact(tyre=self,
                                                  start_node=current_node))
                current_node = self.contacts[-1].penetration_end_node
            current_node = current_node.next
    def find_contact_limits(self):
        pass
    # subcalsses
    class contact:
        def __init__(self,tyre, start_node) -> None:
            self.penetration_start_node:Tyre.node = start_node
            self.penetration_end_node:Tyre.node = None
            self.centre_node:Tyre.node = None
            self.fore_limit:Tyre.node = None
            self.aft_limit:Tyre.node = None
            self.tyre:Tyre = tyre
            self.set_penetration_limits()
        def set_penetration_limits(self):
            self.centre_node = self.penetration_end_node = self.penetration_start_node
            min_distance = np.linalg.norm(np.array([self.tyre.centre_x, self.tyre.centre_y])
                                           - np.array(self.centre_node.penetration_point))
            while self.penetration_end_node.next.penetration_point is not None:
                new_distance = np.linalg.norm(np.array([self.tyre.centre_x, self.tyre.centre_y])
                                           - np.array(self.penetration_end_node.penetration_point))
                if min_distance > new_distance:
                    min_distance = new_distance
                    self.centre_node = self.penetration_end_node
                self.penetration_end_node = self.penetration_end_node.next
            




    class node:
        def __init__(self,tyre, theta, next_node=None, previous_node =None):
            self.tyre:Tyre = tyre
            self.next:Tyre.node = next_node
            self.prev: Tyre.node = previous_node
            self.theta = theta
            self.deformation = 0
            self.x = self.tyre.centre_x + np.cos(self.theta + np.pi/2)*self.tyre.free_radius
            self.y = self.tyre.centre_y + np.sin(self.theta + np.pi/2)*self.tyre.free_radius
            self.penetration_point = None
            self.dy = None
            self.ddy = None
            self.dr_dtheta = None
            self.ddr_dtheta = None
        def __sub__(self, other) -> np.int32:
            result = 0
            current_node = other
            if self is other:
                return 0
            while current_node.next is not other and current_node is not self:
                current_node = current_node.next
                result = result+1
            if current_node is other:
                return None
            return result
        def __add__(self, N:np.int32):
            result_node = self
            i = 0
            while i < N:
                result_node = result_node.next
                i = i+1
            return result_node
        


            
                

            
            
