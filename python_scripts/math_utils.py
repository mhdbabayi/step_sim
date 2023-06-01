import numpy as np
from euclid3 import Vector2

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
def beam_solution(beta,theta,boundary_deformation, boundary_derivative, theta0 = None):
    A = boundary_deformation
    B = (boundary_derivative/beta + boundary_deformation)
    if theta0 is not None:
        A = np.exp(-beta*theta0)*(
            A * np.cos(beta*theta0) + B * np.sin(beta*theta0))
        B = np.exp(-beta*theta0)*(
            A * np.sin(beta*theta0) + B * np.cos(beta*x0))
    return np.exp(-beta*theta)*(A*np.cos(beta*theta) + B*np.sin(beta*theta))    
def beam_force_integral(beta, boundary_deformation, boundary_derivative):
    # radial force: 
    # integral exp(-beta*theta)*(A*cos(beta*theta) + B*sin(beta*theta))*cos(theta)
    # theta from 0 to infinity with a positive beta
    # tangential force: 
    #integral exp(-beta*theta)*(A*cos(beta*theta) + B*sin(beta*theta))*sin(theta)
    # theta from 0 to infinity with a positive beta
    A = boundary_deformation
    B = boundary_derivative/beta + boundary_deformation
    radial_force = beta*(2*(A + B)*beta**2 + A - B)/(4*beta**4+1)
    tangential_force = (A + 2*b*beta**2)/(4*beta**4+1)
    return radial_force , tangential_force
def fit_quadratic(left_point:Vector2,
                  right_point:Vector2,
                  y0:float) -> np.polynomial.Polynomial:
    # return quadratic numpy polynomial that passes through the point 
    # left_point, right_point and (0 , y0)
    #assert left_point.x < 0 and right_point.x > 0
    (xl, yl) = left_point
    (xr , yr)= right_point
    if xl == xr:
        return np.polynomial.Polynomial([y0 , 0 , 0])
    p = np.polynomial.polynomial.Polynomial([
            y0,
            -(xl**2*y0 - xr**2*y0 - xl**2*yr + xr**2*yl)/(xl*xr*(xl - xr)),
            (xl*y0 - xr*y0 - xl*yr + xr*yl)/(xl*xr*(xl - xr))])
    return p



