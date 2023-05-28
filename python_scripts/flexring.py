import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import physics_engine as phsx
from euclid3 import Vector2
from dataclasses import dataclass
from scipy import interpolate

'''
sign convention: 
Tyre node zero is at the top and the nodes go counter-clockwise
Deflection of a node is positive towards the outside of the tyre(increase in radius)
'''
dt = 0.01

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

    # calculate the parameters the ratios of intersetion point distance from 
    # the line segment start point to line segement length
    distance_ratio_line_1 = num1 / den
    distance_ratio_line_2 = num2 / den

    # check if the intersection point is within the line segments
    if distance_ratio_line_1 >= 0 and distance_ratio_line_1 <= 1 and\
          distance_ratio_line_2 >= 0 and distance_ratio_line_2 <= 1:
        return distance_ratio_line_2

    # no intersection found
    return None
def circle_line_intersection(line_start:Vector2, line_end:Vector2, circle_centre:Vector2, radius:float):
    #v_ denotes a vector starting at line start
    v_line = line_end - line_start
    v_centre = line_start - circle_centre
    # quadratic equation constants:
    A = v_line.dot(v_line)
    B = 2*v_centre.dot(v_line)
    C = v_centre.dot(v_centre) - radius**2
    Result = [None, None]
    if (discriminant := B**2 - 4*A*C) <= 0:
        return None
    discriminant = np.sqrt(discriminant)
    t1 = (-B - discriminant)/(2*A)
    t2 = (-B + discriminant)/(2*A)
    if (t1<0) or (t1 > 1):
        if(t2<0) or (t2 > 1):
            return None
    if (t1 > 0 and t1 < 1):
        Result[0] = line_start + t1*v_line
    if(t2 > 0 and t2 < 1):
        Result[1] = line_start + t2*v_line
    return Result
def find_chord_centre(p1:Vector2, p2:Vector2, circle_centre:Vector2):
    #given two points p1 and p2, find the closest point on line segment p1-p2 to circle_centre
    # points p1 and p2 are assumed inside the circle, no checks are carried out. BE CAREFUL!
    v_line = p2 - p1
    v_centre = p1 - circle_centre
    t = -(v_line.dot(v_centre))/(v_line.magnitude_squared())
    t = min(max(t , 0),1)
    return p1 + t*v_line, t
def polar_derivative(point:Vector2, dy):
    #given point x,y and the deriavtive in cartesian coordiantes dy
    #calculates dr/d(theta) in polar coordinates
    x = point.x
    y = point.y
    r = np.sqrt(x**2 + y**2)
    dr = r*(x + y*dy)/(x*dy - y)
    return dr
def polar_second_derivative(point:Vector2, dy , ddy):
    x = point.x
    y = point.y
    r = np.sqrt(x**2 + y**2)
    dr_dtheta = polar_derivative(point, dy)
    dx_dtheta = dr_dtheta*x/r  - y
    dy_dtheta = dr_dtheta*y/r + x
    dydx_dtheta = ddy*dx_dtheta
    ddr_dtheta = ((x*dy - y)*(dx_dtheta + y*dydx_dtheta + dy*dy_dtheta) -\
                  (y*dy + x)*(dy*dx_dtheta + dydx_dtheta*x - dy_dtheta)) * r/(x*dy - y)**2\
                   +dr_dtheta*(y*dy  +x)/(x*dy - y)
    return ddr_dtheta
def fit_poly(P1, P2, P3):
    x1, y1, dy1 = (P1[0], P1[1], P1[2])
    x2, y2, dy2 = (P2[0], P2[1], P2[2])
    x3, y3, dy3 = (P3[0], P3[1], P3[2])

    A = np.array([
        [x1**5, x1**4, x1**3, x1**2, x1, 1],
        [5*x1**4, 4*x1**3, 3*x1**2, 2*x1, 1, 0],
        [x2**5, x2**4, x2**3, x2**2, x2, 1],
        [5*x2**4, 4*x2**3, 3*x2**2, 2*x2, 1, 0],
        [x3**5, x3**4, x3**3, x3**2, x3, 1],
        [5*x3**4, 4*x3**3, 3*x3**2, 2*x3, 1, 0]
    ])

    B = np.array([y1, dy1, y2, dy2, y3, dy3])
    coefficients = np.linalg.solve(A, B)
    
    return np.poly1d(coefficients)
def construct_piecewise_poly(start, end, peak):
    x1, y1 = (start[0],start[1])
    x2, y2 = (end[0], end[1])
    x3, y3 = (peak[0], peak[1])

    a1 = (y1-y3)/((x1-x3)**2)
    a2 = (y2 -y3)/(x2**2 - 2*x2*x3 + x3**2)
    b1 = -(2*x3*y1 - 2*x3*y3)/((x1-x3)**2)
    b2 = -(2*x3*y2 - 2*x3*y3)/(x2**2 - 2*x2*x3 + x3**2)
    c1 = (y3*x1**2 - 2*y3*x1*x3 + y1*x3**2)/((x1 - x3)**2)
    c2 = (y3*x2**2 - 2*y3*x2*x3 + y2*x3**2)/(x2**2 - 2*x2*x3 + x3**2)

    def piecewise_polynomial(x):
        if x <= x3:
            return a1*x**2 + b1*x + c1
        else:
            return a2*x**2 + b2*x + c2

    return piecewise_polynomial
def beam_solution(beta,theta,boundary_deformation, boundary_derivative):
    return np.exp(-beta*theta)*\
            (boundary_deformation*np.cos(beta*theta) +
            (boundary_derivative/beta + boundary_deformation)*np.sin(beta*theta))
def interpolate_boundary_condition(centre:Vector2,
                                   radius:float,
                                   point1:Vector2,
                                   point2:Vector2,
                                   dydx:float,
                                   ddydx:float,
                                   beta:float,
                                   direction:int):
    # given a chord or a semi chord on a circle, finds the 
    # point along the chord where separation happens
    dr_1 = (point1 - centre).magnitude() - radius
    dr_2 = (point2 - centre).magnitude() - radius
    dr_dtheta_1 = polar_derivative(point= point1 - centre,dy=dydx)
    dr_dtheta_2 = polar_derivative(point= point2- centre,dy= dydx)
    ddr_dtheta_1 = polar_second_derivative(point1 - centre,dydx,ddydx)
    ddr_dtheta_2 = polar_second_derivative(point2 - centre, dydx,ddydx)                                           
    # condition is true when lhs > rhs
    # we calculate lhs and rhs at point1 and point2 and linearly interpolate
    # to find the point where they cross
    """
    print(f'''r1: {dr_1}
    r2: {dr_2}\ndr1:{dr_dtheta_1}
    dr2: {dr_dtheta_1}
    ddr1: {ddr_dtheta_1}
    ddr2: {ddr_dtheta_2}''')
    """
    lhs_1 = 0.5*ddr_dtheta_1
    lhs_2 = 0.5*ddr_dtheta_2
    rhs_1 = -2*(beta**2)*(direction*dr_dtheta_1/beta + dr_1)
    rhs_2 = -2*(beta**2)*(direction*dr_dtheta_2/beta + dr_2)
    t = (rhs_1 - lhs_1)/((lhs_2-lhs_1)-(rhs_2 - rhs_1))# where rhs == lhs
    return point1 + t*(point2 - point1), t

class Road:
    def __init__(self, step_width, step_height,step_profile_phase = np.pi, length = 5, high_res=False) -> None:
        self.length = length
        self.x , self.y = self.create_profile(step_width, step_height , step_profile_phase, length)
        self.points = [Vector2(x , y) for x, y in zip(self.x , self.y)]
        self.high_res = high_res
        if high_res:
            self.points = self.over_sample(self.points)
            self.x = [p.x for p in self.points]
            self.y = [p.y for p in self.points]
        self.dydx = np.zeros(len(self.x))
        self.ddydx = np.zeros(len(self.x))
        for i in np.arange(start=1,stop=(len(self.x))-1):
            self.dydx[i] = (self.y[i+1] - self.y[i-1])/(self.x[i+1] - self.x[i-1])
            self.ddydx[i] = (self.y[i+1] + self.y[i-1] - 2*self.y[i])/\
                               ((self.x[i+1]-self.x[i])*(self.x[i] - self.x[i-1])) 
    @staticmethod
    def create_profile(step_width, step_height, step_profile_phase, length):
        x1 = np.array([0, length/2- step_width/2])
        y1 = np.array([0,0])
        x_step = length/2 + np.linspace(-step_width/2 , step_width/2,
                                        np.max((20,np.int32(step_width/0.01))))[1:]
        step_phase = np.linspace(0, step_profile_phase, len(x_step))
        y_step = step_height/2 - (step_height/2)*np.cos(step_phase)
        x2 = np.array([x_step[-1] + 0.01, length])
        y2 = np.array([y_step[-1], y_step[-1]])
        return np.hstack((x1 , x_step , x2)), np.hstack((y1, y_step , y2))
    @staticmethod
    def over_sample(points, distance_step = 0.001):
        cum_distance = [(points[idx+1] - points[idx]).magnitude()
                 for idx in range(len(points)-1)]
        cum_distance.insert(0 , 0)
        cum_distance = np.cumsum(cum_distance)
        x  = [p.x for p in points]
        y = [p.y for p in points]
        uniform_distance = np.arange(start=0, stop = cum_distance[-1], step=distance_step)
        x_interpolator = interpolate.interp1d(x = cum_distance, y = x,kind="linear")
        y_interpolator = interpolate.interp1d(x = cum_distance, y = y , kind="linear")
        return [Vector2(x_interpolator(d) , y_interpolator(d)) for d in uniform_distance]
class SmartRaod(Road):
    def __init__(self, step_width, step_height,step_profile_phase = np.pi, length = 5) -> None:
        # for now we initialize like before, with only a sine bump in the middle
        super().__init__(step_width, step_height,step_profile_phase, length)
        self.over_sampled_points = self.over_sample(self.points)
        self.node_list = []

    def initialize_nodes(self):
        self.node_list = [SmartRaod.Node(parent_road=self,
                                         position=self.points[idx],
                                         dydx=self.dydx[idx],
                                         ddydx=self.ddydx[idx]) for idx in range(len(self.points))]
        self.node_list[0].next = self.node_list[1]
        self.node_list[-1].prev = self.node_list[-2]
        for i in range(len(self.node_list)-2):
            self.node_list[i+1].next = self.node_list[i+2]
            self.node_list[i+1].prev = self.node_list[i]
        
    class Node():
        def __init__(self,
                     parent_road,
                     position:Vector2,
                     dydx:float,
                     ddydx:float,
                     next=None,
                     prev=None) -> None:
            self.parent_road = parent_road
            self.posistion:Vector2 = position
            self.dydx = dydx
            self.ddydx = ddydx
            self.next:SmartRaod.Node = next
            self.prev:SmartRaod.Node = prev

class Tyre(phsx.RigidBody):
    beta = 3
    def __init__(self, initial_x, initial_y,road:Road,
                 free_radius = 1., node_res_deg = 1.,
                 x_speed = 0, y_speed = 0) -> None:
        super().__init__(mass = 50, initial_x=initial_x, initial_y = initial_y,
                       initial_x_dot = x_speed, initial_y_dot = y_speed,constraint_type='101'
                       )
        self.stiffness = 100000.
        self.road = road
        self.free_radius = free_radius
        self.delta_theta = np.deg2rad(node_res_deg)
        self.node_zero = Tyre.node(self,theta =0)
        theta = 0
        # generate the nodes, with node zero at the top
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
            current_node.penetration_point=None
            for i  in range(len(self.road.x[0:-1])):
                t =intersection((self.states.position.x, self.states.position.y),
                                (current_node.x, current_node.y),
                                (self.road.x[i], self.road.y[i]),
                                (self.road.x[i+1], self.road.y[i+1]))
                if t is not None:
                    current_node.penetration_point = np.array(
                        (self.road.x[i] + t*(self.road.x[i+1] - self.road.x[i]),
                        self.road.y[i] + t*(self.road.y[i+1] - self.road.y[i])))
                    current_node.road_dr = np.linalg.norm(
                        [self.states.position.x , self.states.position.y] - 
                        current_node.penetration_point
                        ) - self.free_radius
                    current_node.road_dy = self.road.dydx[i] + t*(self.road.dydx[i+1] - self.road.dydx[i])
                    current_node.road_ddy = self.road.ddydx[i] + t*(self.road.ddydx[i+1] - self.road.ddydx[i])
                    current_node.road_dr_dtheta = polar_derivative(
                        X = current_node.penetration_point[0] - self.states.position.x,
                        Y = current_node.penetration_point[1] - self.states.position.y,
                        DY = current_node.road_dy)
                    current_node.road_ddr_dtheta = polar_second_derivative(
                        X = current_node.penetration_point[0] - self.states.position.x,
                        Y = current_node.penetration_point[1] - self.states.position.y,
                        DY = current_node.road_dy,
                        DDY= current_node.road_ddy)
                    break
            current_node = current_node.next
        #self.update_derivatives()
    def update_derivatives(self):
        current = self.node_zero.next

        while current is not self.node_zero:

            if current.prev.penetration_point is not None and current.next.penetration_point is not None:
                h0 = current.penetration_point[0] - current.prev.penetration_point[0]
                h1 = current.next.penetration_point[0] - current.penetration_point[0]
                y0 = current.prev.penetration_point[1]
                y1 = current.penetration_point[1]
                y2 = current.next.penetration_point[1]
                current.road_dy = (y2 - y0)/((h0+h1))
                current.road_ddy = (y2 + y0 - 2*y1)/(h1*h0)
                current.road_dr_dtheta = polar_derivative(X=current.x - self.states.position.x,
                                                     Y = current.y - self.states.position.y,
                                                     DY = current.road_dy)
                current.road_ddr_dtheta = polar_second_derivative(X = current.x - self.states.position.x,
                                                             Y = current.y -self.states.position.y,
                                                             DY = current.road_dy, DDY= current.road_ddy)
            current = current.next
    def update_contacts(self):
        current_node = self.node_zero.next
        self.contacts = []
        while current_node is not self.node_zero:
            if current_node.penetration_point is not None:
                self.contacts.append(Tyre.contact(tyre=self,
                                                  start_node=current_node))
                current_node = self.contacts[-1].fore_penetration_node
            current_node = current_node.next
    def update_deformation(self):
        current_node = self.node_zero.next
        # set all deformations to zero
        while current_node is not self.node_zero:
            current_node.deformation = 0
            current_node = current_node.next
        #fore:
        for c in self.contacts:
            #fore
            current_node = c.fore_separation_node
            delta_theta = 0
            # boundary conditions 
            bc_1 = current_node.road_dr
            bc_2 = current_node.road_dr_dtheta
            while np.exp(-self.beta*delta_theta) > 0.01 and current_node is not self.node_zero:
                current_node.deformation = current_node.deformation+\
                    beam_solution(
                    beta=self.beta,
                    theta=delta_theta,
                    boundary_deformation=bc_1,
                    boundary_derivative=bc_2
                    )
                current_node = current_node.next
                delta_theta = delta_theta + self.delta_theta
            #aft
            current_node = c.aft_separation_node
            delta_theta = 0
            # boundary conditions 
            bc_1 = current_node.road_dr
            bc_2 = -current_node.road_dr_dtheta
            while np.exp(-self.beta*delta_theta) > 0.01 and current_node is not self.node_zero:
                current_node.deformation = current_node.deformation +\
                    beam_solution(
                    beta=self.beta,
                    theta=delta_theta,
                    boundary_deformation=bc_1,
                    boundary_derivative=bc_2
                    )
                current_node = current_node.prev
                delta_theta = delta_theta + self.delta_theta
            try:
                c.set_deformation_fit()
            except Exception as e:
                print(e)
        #contact patches:
        for c in self.contacts:
            current_node = c.aft_separation_node
            while current_node is not c.fore_separation_node.next:
                current_node.deformation = current_node.road_dr
                #current_node.deformation = 0
                current_node = current_node.next
    def draw(self):
        plt.plot(self.states.position.x , self.states.position.y, 'r*')
        n = self.node_zero.next
        raw_x = []
        raw_y = []
        penetration_x = []
        penetration_y = []
        deformation_x = []
        deformation_y = []
        while n is not self.node_zero:
            raw_x.append(n.x)
            raw_y.append(n.y)
            if n.penetration_point is not None:
                penetration_x.append(n.penetration_point[0])
                penetration_y.append(n.penetration_point[1])
            if n.deformation is not None:
                deformation_x.append(self.states.position.x + np.cos(n.theta + np.pi/2)*(self.free_radius + n.deformation))
                deformation_y.append(self.states.position.y + np.sin(n.theta + np.pi/2)*(self.free_radius + n.deformation))
            n =n.next
        plt.plot(raw_x , raw_y)
        plt.plot(penetration_x ,penetration_y, 'm.')
        plt.plot(deformation_x , deformation_y,marker=".", color="blue")
        [c.draw() for c in self.contacts]
    def update_node_positions(self):
        current_node = self.node_zero
        current_node.x = self.states.position.x + np.cos(current_node.theta + np.pi/2)*self.free_radius
        current_node.y = self.states.position.y + np.sin(current_node.theta + np.pi/2)*self.free_radius
        while (current_node := current_node.next) is not self.node_zero:
            current_node.x = self.states.position.x + np.cos(current_node.theta + np.pi/2)*self.free_radius
            current_node.y = self.states.position.y + np.sin(current_node.theta + np.pi/2)*self.free_radius
    def update_states(self, external_forces=0):
        super().update_states(external_forces)
        self.contacts = []
        self.update_node_positions()
        self.update_penetrations()
        self.update_contacts()
        self.update_deformation()
    def update_forces(self, external_forces):
        super().update_forces(external_forces)
        for c in self.contacts:
            self.forces = self.forces + c.get_forces()
        
    # subcalsses
    class contact:
        def __init__(self,tyre, start_node) -> None:
            self.tyre:Tyre = tyre
            self.aft_penetration_node:Tyre.node = start_node
            self.fore_penetration_node:Tyre.node = None
            self.centre_node:Tyre.node = None
            self.fore_separation_node:Tyre.node = None
            self.aft_separation_node: Tyre.node = None
            self.set_penetration_limits()
            self.set_boundary_conditions()
        def set_penetration_limits(self):
            self.centre_node = self.fore_penetration_node = self.aft_penetration_node
            min_distance = np.linalg.norm(np.array([self.tyre.states.position.x, self.tyre.states.position.y])
                                           - np.array(self.centre_node.penetration_point))
            while self.fore_penetration_node.next.penetration_point is not None:
                new_distance = np.linalg.norm(np.array([self.tyre.states.position.x, self.tyre.states.position.y])
                                           - np.array(self.fore_penetration_node.penetration_point))
                if min_distance > new_distance:
                    min_distance = new_distance
                    self.centre_node = self.fore_penetration_node
                self.fore_penetration_node = self.fore_penetration_node.next
        def set_boundary_conditions(self):
            
            self.fore_separation_node = self.centre_node
            while not self.fore_separation_node.seperation_condition(direction=1) and\
                    self.fore_separation_node.next is not self.fore_penetration_node:
                self.fore_separation_node = self.fore_separation_node.next
            # print('fore done')
            self.aft_separation_node = self.centre_node
            while not self.aft_separation_node.seperation_condition(direction=-1) and\
                    self.aft_separation_node.prev is not self.aft_penetration_node:
                self.aft_separation_node = self.aft_separation_node.prev
        def draw(self):
            plt.plot(self.centre_node.penetration_point[0],
                     self.centre_node.penetration_point[1],
                     marker='o')
            plt.plot(self.fore_separation_node.penetration_point[0],
                     self.fore_separation_node.penetration_point[1],
                     marker='*', color = 'black', markersize=10)
            plt.plot(self.aft_separation_node.penetration_point[0],
                     self.aft_separation_node.penetration_point[1],
                     marker='*', color = 'green', markersize=10)
            plt.plot(self.aft_penetration_node.penetration_point[0],
                     self.aft_penetration_node.penetration_point[1],
                     marker='o', color = 'red', markersize=5)
            plt.plot(self.fore_penetration_node.penetration_point[0],
                     self.fore_penetration_node.penetration_point[1],
                     marker='o', color='green', markersize=5)
            n = self.aft_separation_node

            while(n:=n.next) is not self.fore_separation_node and n.road_dr is not None:
                plt.plot(self.tyre.states.position.x + np.cos(n.theta + np.pi/2)*(self.tyre.free_radius+n.road_dr),
                         self.tyre.states.position.y + np.sin(n.theta + np.pi/2)*(self.tyre.free_radius + n.road_dr),
                         marker="x", color="black")
        def set_deformation_fit(self):
            poly_evaluator = construct_piecewise_poly(
                start=np.array([self.aft_separation_node.next.theta,
                             self.aft_separation_node.next.road_dr]),
                peak = np.array([self.centre_node.theta,
                          self.centre_node.road_dr]),
                end = np.array([self.fore_separation_node.prev.theta,
                          self.fore_separation_node.prev.road_dr])
            )
            n = self.aft_separation_node
            while (n:=n.next) is not self.fore_separation_node.next:
                n.deformation_fit = poly_evaluator(n.theta)
        def get_forces(self):
            total_force = self.centre_node.road_dr * self.tyre.stiffness
            angle = self.centre_node.theta + np.pi/2
            return Vector2(total_force* np.cos(angle), total_force * np.sin(angle))
    class node:
        def __init__(self,tyre, theta, next_node=None, previous_node =None):
            self.tyre:Tyre = tyre
            self.next:Tyre.node = next_node
            self.prev: Tyre.node = previous_node
            self.theta = theta
            self.deformation = 0
            self.deformation_fit = 0
            self.x = self.tyre.states.position.x + np.cos(self.theta + np.pi/2)*self.tyre.free_radius
            self.y = self.tyre.states.position.y + np.sin(self.theta + np.pi/2)*self.tyre.free_radius
            self.penetration_point = None
            self.road_dr = None # amount of penetration
            self.road_dy = None
            self.road_ddy = None
            self.road_dr_dtheta = None
            self.road_ddr_dtheta = None
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
        def seperation_condition(self, direction):
            '''
            direction shoule be 1 for fore and -1 fore aft
            second derivative of the deformation profile at the start
            is equal to :
            -2*beta^2*(w0 + w_prim/beta)
            separation happens when this derivative is smaller than the road profile
            second derivative, i.e., the road is curaving away from the tyre 
            faster than the profile 
            '''
            # print(f'{self.road_ddr_dtheta:0.3f}\t{self.road_dr_dtheta:0.3f}\t'\
            #       f'{-2*(Tyre.beta**2)*(direction*self.road_dr_dtheta/Tyre.beta + self.road_dr):0.3f}\t'\
            #       f'{self.road_dy:0.3f}\t{self.road_ddy:0.3f}' )

            return 0.5*self.road_ddr_dtheta > \
                    -2*(Tyre.beta**2)*(direction*self.road_dr_dtheta/Tyre.beta + self.road_dr)

class Tyre_Continous(phsx.RigidBody):
    beta = 8
    stiffness = 150000.
    def __init__(self, initial_x, initial_y,road:Road,
                 free_radius = 1., node_res_deg = 1.,
                 x_speed = 0, y_speed = 0) -> None:
        super().__init__(mass = 50, initial_x=initial_x, initial_y = initial_y,
                       initial_x_dot = x_speed, initial_y_dot = y_speed,constraint_type='101'
                       )
        self.road = road
        self.free_radius = free_radius
        self.collisions = []
        self.contacts = []
        # nodes are only used for visualisation
        self.node_angles = np.deg2rad(np.linspace(0 , 360, 361)[0:-1])
        self.node_r = 0*self.node_angles + self.free_radius
    def find_new_collisions(self, start_idx=0):
        road_idx = start_idx
        while self.road.points[road_idx].x < self.states.position.x - self.free_radius:
                road_idx += 1
        while self.road.points[road_idx].x < self.states.position.x + self.free_radius:
            if (T := circle_line_intersection(self.road.points[road_idx],
                                              self.road.points[road_idx+1],
                                              self.states.position,
                                              self.free_radius)) is not None:
                self.collisions.append(Tyre_Continous.Collision(start=T[0],
                                                 end=T[1],
                                                 start_road_idx = road_idx,
                                                 end_road_idx = road_idx))
                # if the line crosses in only 1 point
                while self.collisions[-1].end is None:
                    road_idx +=1
                    if (T := circle_line_intersection(self.road.points[road_idx],
                                                      self.road.points[road_idx+1],
                                                      self.states.position,
                                                      self.free_radius)) is not None:
                        self.collisions[-1].end = T[1]
                        self.collisions[-1].end_road_idx = road_idx            
            road_idx +=1           
    def find_new_contacts(self, start_idx = 0):
        self.find_new_collisions(start_idx)
        for c in self.collisions:
            self.contacts.append(Tyre_Continous.Contact(self, c))
            self.collisions.remove(c)
    def draw(self):
        circle_obj = plt.Circle(self.states.position, self.free_radius, fill=False)
        plt.gca().add_patch(circle_obj)
        for c in self.collisions:
            plt.plot(c.start.x, c.start.y , "r*")
            plt.plot(c.end.x, c.end.y, "*", color="black")
        for c in self.contacts:
            c.draw()
    def update_forces(self, external_forces):
        super().update_forces(external_forces)
        for c in self.contacts:
            self.forces = self.forces + c.get_forces_centre_point()
    def update_states(self, external_forces):
        super().update_states(external_forces)
        #self.init_contacts()
        self.update_contacts()
    def update_contacts(self):
        self.find_new_contacts(start_idx=self.contacts[-1].collision.end_road_idx+5)
        self.contacts = [c for c in self.contacts if c.update() is not None]    
    @dataclass
    class Collision:
        start:Vector2
        end:Vector2
        road_idx:int
        def __init__(self, start, end, start_road_idx, end_road_idx):
            self.start = start
            self.end = end
            self.start_road_idx = start_road_idx
            self.end_road_idx = end_road_idx
        def update(self, centre_road_idx, tyre_centre, tyre_radius,road_inst):
            self.end_road_idx = centre_road_idx
            self.start_road_idx = centre_road_idx
            while ((tyre_centre - road_inst.points[self.end_road_idx]).magnitude() < \
                tyre_radius):
                self.end_road_idx += 1
            self.end = road_inst.points[self.end_road_idx]
            while ((tyre_centre - road_inst.points[self.start_road_idx]).magnitude() < \
                tyre_radius):
                self.start_road_idx -= 1
            self.start = road_inst.points[self.start_road_idx]          
    class Contact:
        def __init__(self,
                     tyre,
                     collision):
            self.tyre = tyre
            self.collision = collision
            # case of collision with a flat line
            '''
            # This is only needed for low rest roads
            if collision.start_road_idx == collision.end_road_idx:
                self.centre_point = 0.5*(
                    collision.start + collision.end)
                self.centre_point_idx = self.collision.start_road_idx
                self.centre_point_interp_value = 0.5
            else:
                idx = collision.start_road_idx
                min_distance = (self.tyre.road.points[idx]-
                self.tyre.states.position).magnitude()
                self.centre_point = collision.start
                self.centre_point_idx = collision.start_road_idx
                # find closest point in collision
                while idx != collision.end_road_idx+1:
                    idx += 1
                    if (self.tyre.road.points[idx] - 
                            self.tyre.states.position).magnitude() < min_distance:
                        min_distance = (self.tyre.road.points[idx] - 
                            self.tyre.states.position).magnitude()
                        self.centre_point_idx = idx
                    self.centre_point, self.centre_point_interp_value =\
                        find_chord_centre(self.tyre.road.points[self.centre_point_idx],
                                          self.tyre.road.points[self.centre_point_idx+1],
                                          self.tyre.states.position)
            '''
            idx = collision.start_road_idx
            self.centre_point_idx = idx
            self.centre_point = collision.start
            min_distance = (self.tyre.road.points[idx] - self.tyre.states.position).magnitude()
            while idx != collision.end_road_idx + 1:
                idx +=1
                if (self.tyre.road.points[idx] - 
                            self.tyre.states.position).magnitude() < min_distance:
                        min_distance = (self.tyre.road.points[idx] - 
                            self.tyre.states.position).magnitude()
                        self.centre_point_idx = idx
                        self.centre_point = self.tyre.road.points[idx]

            self.centre_point_angle = np.arctan2(self.centre_point.y - self.tyre.states.position.y,
                                                 self.centre_point.x - self.tyre.states.position.x)
            self.fore_theta = np.linspace(0 , np.deg2rad(90), 90)
            self.aft_theta = np.linspace(0 , np.deg2rad(90) , 90) 
            self.fore_deformation = self.fore_theta * 0
            self.aft_deformation = self.aft_theta * 0
            self.set_boundaries()
            self.set_deformation()
        def draw(self):
            plt.plot(self.centre_point.x , self.centre_point.y , "o")
            plt.plot(self.aft_separation.x, self.aft_separation.y,
             "x", color="magenta", markersize=15)
            plt.plot(self.fore_separation.x, self.fore_separation.y,
            "|", color="magenta", markersize= 15)
            #fore
            theta = self.fore_theta + np.arctan2(
                (self.fore_separation.y - self.tyre.states.position.y),
                (self.fore_separation.x - self.tyre.states.position.x)
             )
            plt.plot(np.cos(theta)*(self.fore_deformation + self.tyre.free_radius)+\
                    self.tyre.states.position.x,
                     np.sin(theta)*(self.fore_deformation+ self.tyre.free_radius)+\
                            self.tyre.states.position.y) 
            #aft
            theta = -self.aft_theta + np.arctan2(
                (self.aft_separation.y - self.tyre.states.position.y),
                (self.aft_separation.x - self.tyre.states.position.x)
             )
            plt.plot(np.cos(theta)*(self.aft_deformation + self.tyre.free_radius)+\
                    self.tyre.states.position.x,
                     np.sin(theta)*(self.aft_deformation+ self.tyre.free_radius)+\
                            self.tyre.states.position.y) 
        def is_boundary_condition(self, idx, direction):
            road_dr_dtheta = polar_derivative(
                                            point = self.tyre.road.points[idx] - self.tyre.states.position,
                                            dy = self.tyre.road.dydx[idx]
                                            )
            road_ddr_dtheta = polar_second_derivative(
                                                    point=self.tyre.road.points[idx] - self.tyre.states.position,
                                                    dy = self.tyre.road.dydx[idx],
                                                    ddy = self.tyre.road.ddydx[idx]
                                                    )
            road_dr = (self.tyre.road.points[idx] - self.tyre.states.position).magnitude() - self.tyre.free_radius
            return 0.5*road_ddr_dtheta > \
                    -2*(self.tyre.beta**2)*(direction*road_dr_dtheta/self.tyre.beta + road_dr)
        def set_boundaries(self):
            #fore
            idx = self.centre_point_idx
            while (not self.is_boundary_condition(idx, 1)):
                idx = idx+1
            self.fore_separation = self.tyre.road.points[idx]
            self.fore_separation_dr_dtheta = polar_derivative(self.fore_separation - self.tyre.states.position,
                                                              self.tyre.road.dydx[idx])
            # aft
            idx = self.centre_point_idx
            while (not self.is_boundary_condition(idx , -1)):
                idx = idx-1
            self.aft_separation = self.tyre.road.points[idx]
            self.aft_separation_dr_dtheta = polar_derivative(self.aft_separation - self.tyre.states.position,
                                                              self.tyre.road.dydx[idx])
        def get_forces_centre_point(self):
            total_force = (self.tyre.free_radius -\
                (self.centre_point - self.tyre.states.position).magnitude()) * self.tyre.stiffness
            return total_force*(self.tyre.states.position - self.centre_point).normalize()
        def set_deformation(self):
            #fore
            bc1 = (self.tyre.states.position - self.fore_separation).magnitude() - self.tyre.free_radius
            bc2 = self.fore_separation_dr_dtheta
            self.fore_deformation = beam_solution(self.tyre.beta,
                                                  self.fore_theta,
                                                  bc1,bc2)
            bc1 = (self.tyre.states.position - self.aft_separation).magnitude() - self.tyre.free_radius
            bc2 = -self.aft_separation_dr_dtheta
            self.aft_deformation = beam_solution(self.tyre.beta,
                                                  self.aft_theta,
                                                  bc1,bc2)                                     
        def update(self):
            if (self.centre_point - self.tyre.states.position).magnitude() > self.tyre.free_radius:
                return None
            road_tangent_vector = self.tyre.road.points[self.centre_point_idx+1] -\
                self.tyre.road.points[self.centre_point_idx-1]
            idx_increment = np.int32(np.sign(
                    road_tangent_vector.dot(self.tyre.states.velocity)))
            while ((self.tyre.road.points[self.centre_point_idx] - self.tyre.states.position).magnitude() >\
                (self.tyre.road.points[self.centre_point_idx+idx_increment] - self.tyre.states.position).magnitude()):
                self.centre_point_idx += idx_increment
            self.centre_point = self.tyre.road.points[self.centre_point_idx]
            self.collision.update(self.centre_point_idx,
                                  self.tyre.states.position,
                                  self.tyre.free_radius,
                                  self.tyre.road)
            self.set_boundaries()
            self.set_deformation()
            return 1

class SprungMass(phsx.RigidBody):
    def __init__(self,
                tyre_inst:Tyre_Continous,
                mass,
                speed_x=0,
                speed_y=0,
                spring_neutral_length = 0.3,
                natural_frequency_hz = 1.5,
                damping_ratio = 0.5):
        self.tyre_inst = tyre_inst
        super().__init__(mass=mass, initial_x=tyre_inst.states.position.x,
                        initial_y = tyre_inst.states.position.y + spring_neutral_length,
                            initial_x_dot = self.tyre_inst.states.velocity.x,
                              initial_y_dot = speed_y,
                            constraint_type = '001')
        self.natural_freq_rad = 360*np.deg2rad(natural_frequency_hz)
        self.damping_ratio = damping_ratio
        self.spring_neutral_length = spring_neutral_length
        self.spring_stiffness = (self.natural_freq_rad)**2 * self.mass
        self.damping_coefficient = 2*self.natural_freq_rad * self.mass * damping_ratio
        self.spring_preload = self.mass * 9.8
        self.spring_force = Vector2(0 ,self.spring_preload) 
        self.damper_force = Vector2(0 , 0)
    def update_forces(self, external_forces):
        super().update_forces()
        self.spring_force = Vector2(0 ,(self.spring_neutral_length - 
        (self.states.position.y - self.tyre_inst.states.position.y)) + self.spring_preload)
        self.damper_force = Vector2(0,(self.tyre_inst.states.velocity.y  - self.states.velocity.y) * self.damping_coefficient)
        gravity_force = Vector2(0 , -self.mass*9.8)
        self.forces =  self.spring_force + self.damper_force + gravity_force       
    def draw(self):
        width = self.tyre_inst.free_radius * 2
        rect = patches.Rectangle((self.states.position.x-width/2, self.states.position.y-width/2),
                                 width=width, height=width/2)
        plt.gca().add_patch(rect)