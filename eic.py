# # EIC Heat Map

# # Modules

import numpy as np
import os
import time
import sys
from tqdm import * #This is a progress bar module
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches


# ### Some quick Primes Code



def prime(n):
    if n == 0 or n == 1:
        return False
    elif n == 2 :
        return True
    else:
        for i in range(2,int(n**0.5)+1):
            if n%i == 0:
                return False   
    return True

def produce_primes( start = 0, end = 2016 ):
    list_of_primes = []
    while start <= end :
        if prime(start) == True :
            list_of_primes.append(start)
        start += 1
    return list_of_primes
    
    
def prime_factorisation_auxillary(number = 2013):
    number = abs(number)
    list_of_prime_factors = []
    for prime in produce_primes(1,int(number**2) +1):
        if number %prime == 0 :
            list_of_prime_factors.append(prime)
    return list_of_prime_factors

def unique_prime_factors(number = 100):
    initial_list = prime_factorisation_auxillary(number)
    complete_list = []
    exponent_list = []
    exponent = 1
    for prime in initial_list:
        while number%(prime**exponent) == 0:
            if prime not in complete_list:
                complete_list.append(prime)
            exponent += 1
        exponent_list.append(exponent-1)
        exponent = 1
    return dict(zip(complete_list,exponent_list))


# ### Polynomial Class



class Polynomial(object):
    def __init__(self, list_of_co_efficients):
        '''
        The list_of_co_efficients = [a, b, c, ..., k] correspond to the co-efficients of
        ax^n + bx^{n-1} + ... + k in that order.
        '''
        while list_of_co_efficients[0] == 0 :
            list_of_co_efficients = list_of_co_efficients[1:]
        self.co_efficients = list_of_co_efficients
        self.degree = len(list_of_co_efficients)-1
        self.equation = ''
        for i in range(0,self.degree):
            self.equation += str(self.co_efficients[i])+'*x**'+str(self.degree - (i)) +' + '
        self.equation += str(self.co_efficients[::-1][0]) 
        return None
    
    def __repr__(self):
        return self.equation
    
    def evaluate(self, x = 10):
        return eval(self.equation)
    
    def quadratic_descriminant(self):
        if self.degree > 2 :
            raise(ValueError('Polynomial has degree %r > 2' % (self.degree) ))
        a,b,c = self.co_efficients        
        return (b**2 - 4*a*c)
    
    def quadratic_roots(self):
        if self.degree != 2 :
            raise(ValueError('Polynomial has degree %r > 2' % (self.degree) )) 
        a,b,c = self.co_efficients
        if self.quadratic_descriminant() >= 0 :
            r_p = (-b + (b**2 - 4*a*c)**0.5)/(2.0*a)
            r_m = (-b - (b**2 - 4*a*c)**0.5)/(2.0*a)
        else :
            r_p = complex(-b/(2.0*a),  (-b**2 + 4*a*c)**0.5/(2.0*a))
            r_m = complex(-b/(2.0*a), - (-b**2 + 4*a*c)**0.5/(2.0*a))          
        return r_p,r_m
    
    def quadratic_rational_irreducibility(self,eps = 10e-60):
        '''
        Checks if quadratics are irreducible over Q. 
        '''
        if self.degree != 2 :
            raise(ValueError('Polynomial has degree %r not 2' % (self.degree) )) 
        descriminant = self.quadratic_descriminant()
        if descriminant < 0 :
            return True
        elif abs(round(descriminant**0.5) -  descriminant**0.5) > eps :
            return True
        return False
    


# #Creating the Grid


class Cell(object):
    def __init__(self, position, polynomial = 0, value = 0):
        self.position = position
        self.value = value
        self.polynomial = polynomial
        return None     


# ### Versatile Square Grid 




def find_extreme_value(polygons, x_y_or_z = 'x', min_or_max = 'max'):
    '''
    An auxillary function for finding mst extreme value for scaling.
    '''
    auxillary_dictionary = dict([('x',0),('y',1),('z',2)])
    list_of_x_y_or_z = []
    for polygon in polygons:
        for vertex in polygon:
            list_of_x_y_or_z.append(vertex[auxillary_dictionary[x_y_or_z]])
    if min_or_max == 'min':
        return min(list_of_x_y_or_z)
    else :
        return max(list_of_x_y_or_z)

def describe_polygon_path(polygon):
    '''
     Polygon should be a list of verticies of form [(a,b),(c,d)....(x,y)]
    '''
    vertices = []
    for vertex in polygon:
        vertices.append(vertex)
    vertices.append(polygon[0])
    codes = [Path.MOVETO] + [Path.LINETO]*(len(polygon)-1) + [Path.CLOSEPOLY]
    return Path(vertices, codes)


def plot_polygons(polygons,colours = ['r','b','g','y','cyan','darkblue','lightblue','aqua']):
    '''
    Polygons should be a list of list of verticies of form [[(a1,b1),(c1,d1)....(x1,y1)],...,[(an,bn),(cn,dn)....(xn,yn)]]
    '''
    kwds = dict(ec='k', alpha = 0.5) #alpha shows transluscency.
    figure = plt.figure()
    figure.set_size_inches(5,5*(find_extreme_value(polygons,'y','max')-find_extreme_value(polygons,'y','min'))/    (find_extreme_value(polygons,'x','max')-find_extreme_value(polygons,'x','min')))
    axis = figure.add_subplot(111)
    list_of_patches = [patches.PathPatch(describe_polygon_path(polygon),                                           facecolor = colours[polygons.index(polygon)%len(colours)],edgecolor = 'black', lw = 1.3,**kwds)                        for polygon in polygons]
    for patch in list_of_patches:
        axis.add_patch(patch)
    plt.axis('scaled')
    axis.set_xlim(find_extreme_value(polygons,'x','min')-1,find_extreme_value(polygons,'x','max')+1)
    axis.set_ylim(find_extreme_value(polygons,'y','min')-1,find_extreme_value(polygons,'y','max')+1)
    axis.grid(True)    
    plt.show()
    return None

def translate(polygon,x_displacement, y_displacement):
    new_polygon = []
    for vertex in polygon:
        new_polygon.append((vertex[0]+x_displacement,vertex[1]+y_displacement))
    return new_polygon

def scale(polygon,scalar):
    new_polygon = []
    for vertex in polygon:
        new_polygon.append((vertex[0]*scalar, vertex[1]*scalar))
    return new_polygon

def rotate(polygon, theta):
    new_polygon = []
    for vertex in polygon:
        new_polygon.append((vertex[0]*np.cos(theta)-vertex[1]*np.sin(theta), vertex[0]*np.sin(theta)+vertex[1]*np.cos(theta)))
    return new_polygon
    
def reflect(polygon, phi):
    new_polygon = []
    for vertex in polygon:
        new_polygon.append((vertex[0]*np.cos(phi/2)+vertex[1]*np.sin(phi/2), vertex[0]*np.sin(phi/2)-vertex[1]*np.cos(phi/2)))
    return new_polygon
    
def produce_regular_polygon(number_of_sides):
    if type(number_of_sides) != int or number_of_sides < 3 :
        raise ValueError( 'Error : need an integer greater than or equal to 3')
    list_of_points = [(0,0)]
    angle = 0
    for side in range(number_of_sides):
        list_of_points.append((list_of_points[::-1][0][0] + np.cos(np.pi*(side) - angle),                               list_of_points[::-1][0][1] + np.sin(np.pi*(side) - angle)))
        angle += (number_of_sides - 2)*np.pi/number_of_sides 
    return list_of_points
    









def eic_quadratic(quadratic,p, k = 0):
    [m,a,b] = quadratic.co_efficients
    if prime(p) == False:
        raise ValueError('%d is not a prime.'%p)
    if m%p != 0 and (2*m*k+a)%p == 0 and (m*k**2 + a*k+b)%p == 0 and (m*k**2 + a*k+b)%p**2 != 0:
        return True
    return False
    
        
def find_correct_primes(quadratic,k = 0):
    [m,a,b] = quadratic.co_efficients
    primes_b = unique_prime_factors(k**2 +a*k +b)
    primes_m = unique_prime_factors(m)
    primes_to_consider = [p for p in primes_b if (p not in primes_m) and primes_b[p] == 1]
    return primes_to_consider

def eic_quadratic_exhaustive_check(quadratic,max_k = 10):
    [m,a,b] = quadratic.co_efficients
    min_pos_k = -1
    consider_triggers = []
    if quadratic.quadratic_rational_irreducibility() == False:
        return 'Reducible'
    
    for k in range(max_k+1):
        for p in find_correct_primes(quadratic, k):
            if eic_quadratic(quadratic,p, k):
                min_pos_k = k
                consider_triggers.append(abs(min_pos_k))
                break
                 
    second_limit = max_k
    if min_pos_k != -1 :
        second_limit = min_pos_k
    for k in range(second_limit):
        for p in find_correct_primes(quadratic, k):
            if eic_quadratic(quadratic,p, k):
                min_neg_k = k
                consider_triggers.append(abs(min_neg_k))
                break
    if consider_triggers:
        return min(consider_triggers)
    return 'no k worked within [-%d, %d]'%(max_k,max_k)
                
            



def produce_monic_grid_eisenstein(m, span = 10,max_search = 8 ):
    '''
    taking absloute value right now!
    '''
    list_of_cells = []
    for a in range(-span,span+1):
        for b in range(-span,span+1):
            list_of_cells.append(Cell((a,b),Polynomial([m,a,b]),value = 0))
    for cell in tqdm(list_of_cells):
        x = eic_quadratic_exhaustive_check(cell.polynomial,max_search)
        if type(x) == int or type(x)== float:
            cell.value = [float(x)/max_search,0,1]
        elif max(x) == 'u':
            #reducible
            cell.value  = [0,0,0]
        elif max(x) == 'w':
            #not eic in range
            cell.value = [1,1,1]
    kwds = dict(ec='k', alpha = 0.9) #alpha shows transluscency.
    figure = plt.figure()
    figure.set_size_inches(8,8)
    axis = figure.add_subplot(111)
    list_of_patches = [patches.PathPatch(describe_polygon_path(translate(produce_regular_polygon(4),cell.position[0],cell.position[1])),                          facecolor = cell.value,edgecolor = 'black', lw = 1.3,**kwds)                        for cell in list_of_cells]       
    for patch in list_of_patches:
        axis.add_patch(patch)
    plt.axis('scaled')
    axis.set_xlim(-span,span+1)
    axis.set_ylim(-span,span+1)
    axis.grid(True)  
    plt.savefig('eisenstein'+str(m))
    plt.show()
    plt.close()


# ## Producing Image 


for k in range(1,2): # pick which leading coefficients range. 
    produce_monic_grid_eisenstein(m = k, span = 50,max_search=20)


# In[ ]:
