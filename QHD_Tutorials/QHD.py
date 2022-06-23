""" QHD Natural Variables Derivation """

from pysces import lin_mom
from sympy import Symbol, symbols, Function, init_printing, Integral, Derivative, sympify, expand, sqrt, sin, cos, pi, simplify, exp
from sympy.physics.quantum import Commutator, Operator, Bra, Ket
from sympy.plotting import plot, plot3d
from sympy.abc import *
from sympy.core import *
import math
from sympy import init_printing
init_printing()

hbar = Symbol("hbar")


def ham(p, q):
    p, q, v, mass = symbols("p q v mass")
    v = Function("v")
    return Operator(lin_mom(q, 2)/(2*mass)) + Operator(v(q))

def QHD_int(n, order, dt):
    xp, mass, p, alpha, x2 = symbols("xp, mass, p, alpha, x2")
    if n == xp:
        p2, alpha= symbols("p2, alpha")
        return 0.5*dt*(2.0*D*alpha*(x*(-2*x**2 + 3.0*x2) - x2) - alpha*x*(-2.0*p**2 + p2)/mass) + (0.5*dt*(2.0*D*alpha*(x*(-2*x**2 + 3.0*x2) - x2) - alpha*x*(-2.0*p**2 + p2)/mass) + xp)*exp(-2.0*alpha*dt*p/mass)
#    if n == x and order == 1:
#        xp, mass, alpha = symbols("xp, mass, alpha")
#        return -alpha*dt*xp/mass + x
    if n == x and order == 2:
        xp, mass, p, alpha, x2 = symbols("xp, mass, p, alpha, x2")
        return -2.0*alpha*dt*x*(-p*x + xp)/mass + (-2.0*alpha*dt*x*(-p*x + xp)/mass + x2)*exp(-2.0*alpha*dt*p/mass)
    else:
        return N(sympify(str(n**order + time_deriv(n, order)*dt).replace("hbar*i", "1").replace("I", "1").replace("f(q)", "1").replace("hbar**2*i**2*Derivative(1, q)", "p").replace("hbar**2*i**2", "1")))

def QHD_kinetic(var = None):
    return QHD_int(p, 2)/(2*m)   
    
def time_deriv1(var, order = 1):
    pq_s = Symbol("pq_s")
    aux = Operator(Function("f")(q))
    pq_s = (Operator(p)*Operator(q)+Operator(q)*Operator(p))/2
    
    h1 = (expand(((Commutator(Operator(var**(order)), ham(p, q)).expand(commutator=True))*aux).doit()))

    
    if var == V:
        der = f'((p)*(Derivative(V, q, {order}+1))+(Derivative(V, q, {order}+1))*(p))/2/m'
    
    if var == p:
        p1 = lin_mom(q, order)
        h1 = (str(simplify(h1)).replace("p", str(Operator(p1))))       
    
    else:
        h1 = str(expand(((Commutator(ham(p, q), Operator(var**(order))).expand(commutator=True))*aux).doit()))
        h1 = h1.replace("pq_s", "((p*q+q*p))/2")


    
    start_points = []
    end_points = []

    start = 0

    s = h1.find("Derivative(1, q)", start)
    while s != -1:
        s = h1.find("Derivative(1, q)", start)
                
        if var == pq_s:
            start = s + len(f"Derivative(1, q) ")
        else:
            start = s + len(f"Derivative(1, q)")
        
        start_points.append(s)
        if s == -1:
            break
               
        if var == pq_s:
            e = h1.find(" ", start + 1)
        else:
            e = h1.find("f(q)", start) + 4            
        
        end_points.append(e)
    
   
    if start_points[-1] == -1:
        start_points.remove(start_points[-1])
    
    new_derivative_function = []
    replace_spot = []

    for i in range(len(start_points)):
        func = h1[start_points[i]:end_points[i]]
        func = start_points[i] + len(f"Derivative(1, q)")
        func1 = h1.find("*", func, end_points[i])

        repl = h1[start_points[i]:end_points[i]].find("1") + start_points[i]
        replace_spot.append(repl)

        
        if var == p:
            new_func = h1[func1:end_points[i]]
            new_derivative_function.append(new_func)    

        else:
            new_func = h1[func1 - 3:end_points[i]]
            new_derivative_function.append(new_func)

    temp = h1.split(" ")

    position_list = []
    temp_list = []


    for i in range(len(temp)):
        if var == p:
            st = 1
        else:
            st = 0
            
        for r in range(st, len(new_derivative_function)):

            if temp[i].find(new_derivative_function[r]) != -1:
                temp[i] = temp[i].replace(new_derivative_function[r], "")

        if temp[i].find("Derivative") == -1:
            continue

        else:
            position = temp[i].find("1,")
            if position != -1:
                position_list.append(position)
                temp_list.append(i)


    TL = len(temp_list) - 1
    nested_list = []
    while TL > -1: 
        r = [temp_list[TL], position_list[TL]]
        nested_list.append(r) 
        TL -= 1 
        if TL == -1:
            break

    nested_list.reverse()

    if var == p:
        for i in range(len(nested_list)):
            temp[nested_list[i][0]] = \
                temp[nested_list[i][0]].replace(temp[nested_list[i][0]][nested_list[i][1]], new_derivative_function[i][1:])  
    
    else:
        for i in range(len(nested_list)):
            temp[nested_list[i][0]] = \
                temp[nested_list[i][0]].replace(temp[nested_list[i][0]][nested_list[i][1]], "1/" + new_derivative_function[i][1:])
        
        
    string = " ".join(temp)
    if "(*" in string:
        string = string.replace("(*", "(")
    s = sympify(string)
    s1 = expand(s.doit())
    if var == p:
        repl = 1
        if order == 2:
            s1 = s1/2
    else:
        repl = 0
    s2 = str(s1).replace("hbar*i", "1").replace("I", "1").replace("f(q)", f"{repl}").replace(f"hbar**2*i**2*Derivative({repl}, q)", "p").replace(f"Derivative({repl}, q)", "p").replace("hbar**2*i**2", "1").replace("Derivative(v(q), (q, 2))", "0")
  
        
    if var == V:
        return sympify(der)    #/(i*hbar)
    else:
        return sympify(s2)              #/(i*hbar)


def time_deriv(var, order = 1):
    return time_deriv1(var, order)#/(i*hbar)

def symmetrize(expr):
    p, x2, x, px, xp = symbols("p, x2, x, px, xp")
    expr1 = str(expr)
    expr1 = expr1.replace("p*x2", "(x2*p-2*x*x*p+2*x*px)").replace("p*x", "xp").replace("px", "xp")
    return sympify(expr1)