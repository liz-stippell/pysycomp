#!/usr/bin/env python
# coding: utf-8

# In[3]:


from sympy import Symbol, symbols, Function, init_printing, Integral, Derivative, sympify, expand, sqrt, sin, cos, pi, simplify, exp
from sympy.physics.quantum import Commutator, Operator, Bra, Ket
from sympy.plotting import plot, plot3d
from sympy.abc import *
from sympy.core import *
import math
from sympy import init_printing
init_printing()


I = sqrt(-1)
h_b = Symbol("hbar")



def find_variable(aux):
    aux = str(aux)
    new = aux.find("(")
    return sympify(aux[new:])



def comm_1(commutator_1, commutator_2, aux):
    """
    This function is not used directly, it is only used in the comm() function below. 

    Args:
        commutator_1: The first operator in the commutator
        commutator_2: The second operator in the commutator
        aux: The auxiliary function. This is defined below, as F(x).

    Returns:
        The commutator of commutator_1 and commutator_2 with respect to the auxiliary function.

    """
    
    return expand((Commutator(Operator(commutator_1), Operator(commutator_2))*aux).doit())




def comm(commutator_1, commutator_2, aux):
    L_z, L_y, L_x = symbols("L_z, L_y, L_x")
    """
    This function has a few different outputs, depending on the parameters.

    Args:
        commutator_1: the first operator
        commutator_2: the second operator
        aux: the auxiliary function

    Returns:
        This function automatically returns a solved commutator. For more complicated commutators involving the angular momentum operators (L_z, L_y, L_x), please make sure that the correct notation is being used (ex: if you want the angular momentum operator in the "x" direction, please use "L_x"
    
    
    Note:    
        x, p_y, y, p_x, z, p_z, L_z, L_y, L_x = symbols("x, p_y, y, p_x, z, p_z, L_z, L_y, L_x")
        L_z = x*p_y - y*p_x
        L_y = z*p_x - x*p_z
        L_x = y*p_z - z*p_y

    """
    
    if commutator_1 == L_z and commutator_2 == L_y:
        x, p_y, y, p_x, z, p_z, L_z, L_y, L_x = symbols("x, p_y, y, p_x, z, p_z, L_z, L_y, L_x")
        return (Commutator(Operator(x*p_y), Operator(z*p_x))*aux + Commutator(Operator(y*p_x), Operator(x*p_z))*aux - Commutator(Operator(x*p_y), Operator(x*p_z))*aux - Commutator(Operator(y*p_x), Operator(z*p_x))*aux)
    if commutator_1 == L_z and commutator_2 == L_x:
        x, p_y, y, p_x, z, p_z, L_z, L_y, L_x = symbols("x, p_y, y, p_x, z, p_z, L_z, L_y, L_x")
        return (Commutator(Operator(x*p_y), Operator(y*p_z))*aux + Commutator(Operator(y*p_x), Operator(z*p_y))*aux - Commutator(Operator(x*p_y), Operator(z*p_y))*aux - Commutator(Operator(y*p_x), Operator(y*p_z))*aux)
    if commutator_1 == L_y and commutator_2 == L_z:
        x, p_y, y, p_x, z, p_z, L_z, L_y, L_x = symbols("x, p_y, y, p_x, z, p_z, L_z, L_y, L_x")
        return (Commutator(Operator(z*p_x), Operator(x*p_y))*aux + Commutator(Operator(x*p_z), Operator(y*p_x))*aux - Commutator(Operator(z*p_x), Operator(y*p_x))*aux - Commutator(Operator(x*p_z), Operator(x*p_y))*aux)
    if commutator_1 == L_y and commutator_2 == L_x:
        x, p_y, y, p_x, z, p_z, L_z, L_y, L_x = symbols("x, p_y, y, p_x, z, p_z, L_z, L_y, L_x")
        return (Commutator(Operator(z*p_x), Operator(y*p_z))*aux + Commutator(Operator(x*p_z), Operator(z*p_y))*aux - Commutator(Operator(z*p_x), Operator(z*p_y))*aux - Commutator(Operator(x*p_z), Operator(y*p_z))*aux)
    if commutator_1 == L_x and commutator_2 == L_z:
        x, p_y, y, p_x, z, p_z, L_z, L_y, L_x = symbols("x, p_y, y, p_x, z, p_z, L_z, L_y, L_x")
        return (Commutator(Operator(y*p_z), Operator(x*p_y))*aux + Commutator(Operator(z*p_y), Operator(y*p_x))*aux - Commutator(Operator(y*p_z), Operator(y*p_x))*aux - Commutator(Operator(z*p_y), Operator(x*p_y))*aux)
    if commutator_1 == L_x and commutator_2 == L_y:
        x, p_y, y, p_x, z, p_z, L_z, L_y, L_x = symbols("x, p_y, y, p_x, z, p_z, L_z, L_y, L_x")
        return (Commutator(Operator(y*p_z), Operator(z*p_x))*aux + Commutator(Operator(z*p_y), Operator(x*p_z))*aux - Commutator(Operator(y*p_z), Operator(x*p_z))*aux - Commutator(Operator(z*p_y), Operator(z*p_x))*aux)
    
    if commutator_1 == L_z:
        x, p_y, y, p_x, z, p_z, L_z, L_y, L_x = symbols("x, p_y, y, p_x, z, p_z, L_z, L_y, L_x")
        return (Commutator(Operator(x*p_y), Operator(commutator_2))*aux - Commutator(Operator(y*p_x), Operator(commutator_2))*aux)
    if commutator_1 == L_y:
        x, p_y, y, p_x, z, p_z, L_z, L_y, L_x = symbols("x, p_y, y, p_x, z, p_z, L_z, L_y, L_x")
        return (Commutator(Operator(z*p_x), Operator(commutator_2))*aux - Commutator(Operator(x*p_z), Operator(commutator_2))*aux)
    if commutator_1 == L_x:
        x, p_y, y, p_x, z, p_z, L_z, L_y, L_x = symbols("x, p_y, y, p_x, z, p_z, L_z, L_y, L_x")
        return (Commutator(Operator(y*p_z), Operator(commutator_2))*aux - Commutator(Operator(z*p_y), Operator(commutator_2))*aux)
    if commutator_2 == L_z:
        x, p_y, y, p_x, z, p_z, L_z, L_y, L_x = symbols("x, p_y, y, p_x, z, p_z, L_z, L_y, L_x")
        return (Commutator(Operator(commutator_1), Operator(x*p_y))*aux - Commutator(Operator(commutator_1), Operator(y*p_x))*aux)
    if commutator_2 == L_y:
        x, p_y, y, p_x, z, p_z, L_z, L_y, L_x = symbols("x, p_y, y, p_x, z, p_z, L_z, L_y, L_x")
        return (Commutator(Operator(commutator_1), Operator(z*p_x))*aux - Commutator(Operator(commutator_1), Operator(x*p_z))*aux)
    if commutator_2 == L_x:
        x, p_y, y, p_x, z, p_z, L_z, L_y, L_x = symbols("x, p_y, y, p_x, z, p_z, L_z, L_y, L_x")
        return (Commutator(Operator(commutator_1), Operator(y*p_z))*aux - Commutator(Operator(commutator_1), Operator(z*p_y))*aux)
    else:
        return (expression_replace(comm_1(Operator(commutator_1), Operator(commutator_2), aux), find_variable(aux)))
    
    
    
    
def factorization(expr, var):
    """
    Args:    
        expr: the expression of interest
        var: the variable of interest. This is most likely going to be the same variable used in the auxiliary function.
        
    Returns:    
        The simplified commutator. This is only used when both of the operators in the commutator are angular momentum operators.
    
    """
    
    L_z, L_y, L_x = symbols("L_z, L_y, L_x")
    if var == z:
        return comm(Operator(z), lin_mom(z), f(z))*L_z
    if var == y:
        return comm(Operator(y), lin_mom(y), f(y))*L_y
    if var == x:
        return comm(Operator(x), lin_mom(x), f(x))*L_x

        
    
def comm_steps(commutator_1, commutator_2, aux):
    """
    Args:    
        commutator_1: the first operator
        commutator_2: the second operator
        aux: the auxiliary function
    
    Returns:    
        Three main steps used to solve a commutator. The first is the printed commutator in brackets, the second is the expansion, and the 
        third output will be the answer. Please note that as of now, the code does not understand the 
        addition/subtraction/multiplication/division of commutators, and so if there are multiple commutators in a line (ex: 
        comm_steps(comm(A, B) - comm(C, D))) the code will not process this and therefore each individual commutator needs to have its own 
        line.
        
    """
    
    return display(Commutator(Operator(commutator_1), Operator(commutator_2))*aux),\
                    display(comm_1(Operator(commutator_1), Operator(commutator_2), aux)),\
                            display(expression_replace(comm_1(Operator(commutator_1), Operator(commutator_2), aux), find_variable(aux)))



#def X(expr, n=1):
#    """
#    
#    Args:
#        expr: The expression of interest
#        n: Variable to call a different output when higher powers are used
#        
#    Returns:
#        The operator x
#    """
#    x = Symbol("x")
#    if n==1:
#        return simplify(x*expr)
#    else:
#        return simplify(x*X(expr, n-1))

    
    
def lin_mom(x, exp = 1):
    """
    
    Args:
        x: The variable in which the linear momentum operator is with respect to
        exp: The exponent of the linear momentum operator
    
    Returns:
        The linear momentum (p) operator for a selected variable (x, y, etc.)
    """
    
    if exp == 1:
        return (-i*h_b*Derivative("1", x))
    else:
        return (-i*h_b*Derivative( lin_mom(x, exp - 1), x))



#def lin_mom(var = None):
#    """
#
#    Args:
#        var: The variable that the linear momentum operator is with respect to.
#
#    Returns:
#        The linear momentum operator, with respect to the parameter.
#    
#    Note:
#        The "1" in the derivative is a placeholder, which will be replaced.
#        If var == None, the general linear operator is printed, with all three positional arguments "x", "y", and "z"
#
#    """
#    
#    h_b, m = symbols("h_b m")
#    if var == None:
#        return Operator(-I*h_b*(Derivative("1", x))) + Operator(-I*h_b*(Derivative("1", y))) + Operator(-I*h_b*(Derivative("1", z)))
#    else:
#        return Operator(-I*h_b*(Derivative("1", var)))



def kinetic_energy(var = None):
    """

    Args:
        var: The variable in which the derivative is with respect to.

    Returns:
        The kinetic energy operator, with respect to the chosen parameter. 
        
    Note:
        The "1" in the derivative is a placeholder, which will be replaced.
        If var == None, the general linear operator is printed, with all three positional arguments "x", "y", and "z"

    """
    
    m = Symbol("m")
    if var == None:
        return (Operator((-I*h_b*(Derivative("1", x)))**2)*(1/(2*m))) + (Operator((-I*h_b*(Derivative("1", y)))**2)*(1/(2*m))) + (Operator((-I*h_b*(Derivative("1", z)))**2)*(1/(2*m)))
    else:
        return (Operator((-I*h_b*(Derivative("1", var)))**2)*(1/(2*m)))



def potential_energy(var):
    """
    
    Args:
        var: The variable for the given function. This is usually "x", "y" or "z".
    
    Returns:
        The general potential energy operator, v(x)
    
    """
    
    return Operator(Function("v")(var))



def hamiltonian(var):
    """

    Args:
        var: The variable for the given function. This is usually "x", "y" or "z".
    
    Returns:    
        The Hamiltonian operator, made up of the kinetic energy operator and the general potential energy operator.
    
    """
    
    return kinetic_energy(var) + v(var)



def expression_replace(expr, var):
    """
    
    This is only used within the "comm()" function. 
    
    expression_replace(expr, var)

    Args:
        expr: The expanded commutator to be replaced
        var: The variable/parameter with respect to the chosen commutator.

    Returns:
        This replaces the Derivative(1, x) present in the expanded commutator with either Derivative(F(x), x) or Derivative(x*F(x), x).
        
    Note: 
        The above states "x", but can be done for x, y, or z variables.
        The "expr" parameter is a str()


    """
    
    L_x, L_z, L_y, p_x, p_y, p_z, x, y, z = symbols("L_x, L_z, L_y, p_x, p_y, p_z, x, y, z")
    
    if str('[p_x,x]') in str(expr):
        return sympify(str(sympify(str(expr).replace(str('[p_x,x]'), str(expression_replace(comm(lin_mom(var), Operator(var), f(var)), K))))).replace(str(p_y*z - p_z*y), str(-L_x)))
    
    if str('[p_y, y]') in str(expr):
        return sympify(str(sympify(str(R).replace(str('[p_y,y]'), str(expression_replace(comm(lin_mom(var), Operator(var), f(var)), K))))).replace(str(p_x*z - p_z*x), str(-L_y))) 
    
    if str('[p_z,z]') in str(expr):
        return sympify(str(sympify(str(expr).replace(str('[p_z,z]'), str(expression_replace(comm(lin_mom(var), Operator(var), f(var)), K))))).replace(str(p_x*y - p_y*x), str(-L_z)))
    
    else:
        
        expr = str(expr)
        
        starting_pos = expr.find(f"Derivative(1, {var})")
        second = expr[starting_pos + 17:].find(f"Derivative(1, {var})")


        w = expr[starting_pos + 17:].find(f"f({var})")
        wn = starting_pos + 17 + w + len(f"f({var})")
        f = expr[starting_pos + second + len(f"Derivative(1, {var})")*2:].find(f"f({var})")
        fn = starting_pos + second + len(f"Derivative(1, {var})")*2 + f + len(f"f({var})")

        k = f"Derivative(1, {var}), {var})"

        if k in expr:
            term1 = expr[starting_pos + 21:wn]
            term2 = expr[starting_pos + second + 17 + 21:fn]
            one = f"Derivative({term1}, {var}), {var})"
            two = f"Derivative({term2}, {var}), {var})"

        else:
            term1 = expr[starting_pos + len(f"Derivative(1, {var})") + 1:wn]
            if len(f"{var}") == 1:
                term2 = expr[starting_pos + second + len(f"Derivative(1, {var})")*2 + 2:fn]
            if len(f"{var}") == 2:
                term2 = expr[starting_pos + second + len(f"Derivative(1, {var})")*2 + 1:fn]
            if len(f"{var}") == 3:
                term2 = expr[starting_pos + second + len(f"Derivative(1, {var})")*2 :fn]

            one = f"Derivative({term1}, {var})"
            two = f"Derivative({term2}, {var})"

        x1 = expr.replace(expr[starting_pos + second + 17: fn], two).replace(expr[starting_pos: wn], one)
        x2 = sympify(x1).doit()
        x3 = simplify(x2.replace(i**2, -1))
    
        return x3



def f(var):
    """

    Args:
        var: This is what the auxiliary function is with respect to. It should also match the parameters of the other arguments in the 
        comm() function. (example: if one operator is lin_mom(y), the auxiliary function should be f(y).

    Returns:
        This is the auxiliary function, commonly used in the comm() function.
        This simply returns "f(x)", x being with chosen parameter.

    """
    
    return Operator(Function('f')(var))




def planewave(x):
    """

    Args:
        x: What the function is with respect to.

    Returns:
        The PlaneWave function with respect to the parameter.

    """

    k, omega, t = symbols("k omega t")
    return exp(I*(k*x-omega*t))



def PIB(x, L, n):
    """

    Args:
        x: a variable.
        L: Length of the box.
        n: an integer.

    Returns:
        The WaveFunction for Particle in a Box with respect to the chosen parameters.

    """
    
    return sin((n*pi*x)/L)




def PIB_normalized(x, L, n):
    """

    Args:
        x: a variable.
        L: Length of the box.
        n: an integer.

    Returns:
        The normalized WaveFunction for Particle in a Box, with respect to the chosen variables. This answer can also be calculated using 
        the normalization_constant() function.

    """
    
    return sqrt(2/L)*sin((n*pi*x)/L)




def PIB_2D(x, y, L_x, L_y, n):
    """
    
    Args:
        x: a variable.
        y: a variable.
        L_x: Length of the box in the x direction
        L_y: Length of the box in the y direction
        n: an integer.
    
    Returns:
        The WaveFunction for Particle in a Box in two dimensions with respect to the chosen parameters.
    
    """
    return  sin((n*pi*x)/L_x)*sin((n*pi*y)/L_y)

def PIB_2D_normalized(x, y, L_x, L_y, n):
    """
    
    Args:
        x: a variable.
        y: a variable.
        L_x: Length of the box in the x direction
        L_y: Length of the box in the y direction
        n: an integer.
    
    Returns:
        The normalized WaveFunction for Particle in a Box, with respect to the chosen variables. This answer at this time cannot be 
        calculated using the normalization_constant() function.
    
    """
    return sqrt(2/L_x)*sin((n*pi*x)/L_x)*sqrt(2/L_y)*sin((n*pi*y)/L_y)




def moving_gaussian(alpha, x_0, p):
    """

    Args:
        alpha:alpha parameter of the gaussian
        x_0: x_0 parameter of the gaussian
        p: p parameter of the gaussian
    
    Returns:
        The moving gaussian wave function.
    
    """
    
    alpha, x, p, h_b = symbols("alpha x p h_b")
    return exp(-((alpha/2)*(x-x_0)**2 + (I*p)/h_b*(x-x_0)))




def fixed_gaussian(alpha, x_0, p):
    """

    Args:
        alpha:alpha parameter of the gaussian
        x_0: x_0 parameter of the gaussian
        p: p parameter of the gaussian
    
    Returns:
        The moving gaussian wave function.
    
    """
    
    alpha, x, p, h_b = symbols("alpha x p h_b")
    return exp(-((alpha/2)*(x-x_0)**2))




def moving_gaussian_normalized(alpha, gamma, x_0, p_0):
    """

    Args:
        alpha: alpha parameter of the normalized gaussian
        gamma: gamma parameter of the normalized gaussian
        x_0: x_0 parameter of the normalized gaussian
        p_0: p_0 parameter of the normalized gaussian
    
    Returns:
        The normalized gaussian WaveFunction, with respect to the chosen variables.
    
    """
    
    alpha, gamma, x, p_0, h_b = symbols("alpha gamma x p_0 h_b")
    return ((2*alpha)/pi)**(1/4) * exp(-alpha*(x - x_0)**2 + ((I*p_0)/h_b)*(x - x_0) + (I*gamma)/h_b)




def fixed_gaussian_normalized(alpha, x_0, p):
    """

    Args:
        alpha:alpha parameter of the gaussian
        x_0: x_0 parameter of the gaussian
        p: p parameter of the gaussian
    
    Returns:
        The moving gaussian wave function.
    
    """
    
    alpha, x, p, h_b = symbols("alpha x p h_b")
    return (alpha/pi)**(1/4)*exp(-((alpha/2)*(x-x_0)**2))




def conjugate(expr):
    """

    Args:
        x: The term of interest. This is commonly a WaveFunction.

    Returns:
        The complex conjugate of a WaveFunction. If there are imaginary terms, they get negated, and if there are no imaginary terms, the 
        WaveFunction is not affected.

    """
    
    return expr.replace(I, -I)




def normalization_constant(wavefunc, lower, upper, var1, var2 = None):
    """

    Args:
        wavefunc: The WaveFunction/expression of interest
        lower: The lower bound of the integral. If bounds are not listed, this is -oo
        upper: The upper bound of the integral. If bounds are not listed, this is oo
        var1: What the integral is taken with respect to
        var2: What the second integral is taken with respect to 

    Returns:
        The normalization constant, with respect to the given parameters. To find the normalized WaveFunction, the output for 
        normalization_constant() must be multiplied by the original WaveFunction. A normalized WaveFunction indicates that the probability 
        of finding a particle within certain bounds must be equal to one.
    
    Note:
        The double integration only works if they are integrated over the same bounds. 

    """
    nreps = 2
    initial = [ sin(n*pi), cos(n*pi)]
    final = [0, 1]

    res = 1/sqrt(Integral(wavefunc*conjugate(wavefunc), (var1, lower, upper)).doit())

    if var2 == None:
        return simplify(res)
    else:
        res2 = 1/sqrt(Integral(wavefunc*conjugate(wavefunc), (var1, lower, upper), (var2, lower, upper)).doit())

    for i in range(nreps):
        res = res.replace(initial[i], final[i])
        res2 = res2.replace(initial[i], final[i])
        return simplify( res2 )
    # n = Symbol("n")
    # return simplify( 1/sqrt(Integral(wavefunc*conjugate(wavefunc), (var, lower, upper)).doit().replace(sin(n*pi), 0).replace(cos(n*pi), 1)) )





def normalization_constant_steps(wavefunc, lower, upper, var1, var2 = None):
    """
    
    Args:
        wavefunc: The WaveFunction/expression of interest
        lower: The lower bound of the integral. If bounds are not listed, this is -oo
        upper: The upper bound of the integral. If bounds are not listed, this is oo
        var1: What the (first) integral is taken with respect to
        var2: What the second integral is taken with respect to (if there are two integrals)
    
    Returns:
        The steps on how to solve a normalization problem.
    """
    if var2 == None:
        return display(1/sqrt(Integral(wavefunc*conjugate(wavefunc), (var1, lower, upper)))), \
                display(1/sqrt(Integral(wavefunc*conjugate(wavefunc), (var1, lower, upper)).doit())), \
                    display(simplify( 1/sqrt(Integral(wavefunc*conjugate(wavefunc), (var1, lower, upper)).doit().replace(sin(n*pi), 0).replace(cos(n*pi), 1)) ))
    else:
        return display(1/sqrt(Integral(wavefunc*conjugate(wavefunc), (var1, lower, upper), (var2, lower, upper)))), \
                display(1/sqrt(Integral(wavefunc*conjugate(wavefunc), (var1, lower, upper), (var2, lower, upper)).doit())), \
                    display(simplify( 1/sqrt(Integral(wavefunc*conjugate(wavefunc), (var1, lower, upper), (var2, lower, upper)).doit().replace(sin(n*pi), 0).replace(cos(n*pi), 1)) ))
    



def expectation_value(wavefunc_1, operator, wavefunc_2, lower, upper, var):
    """
    Computes a symbolic expression for an expectation value of an operator `operator`
    with the two wavefunctions: `<WaveFunc_1|operator|WaveFunc_2>` 

    Args:    
        wavefunc_1: The "bra" normalized WaveFunction
        operator: The operator of interest
        wavefunc_2: The "ket" normalized WaveFunction
        lower: The lower bound of the integral. If bounds are not listed, this is -oo
        upper: The upper bound of the integral. If bounds are not listed, this is oo
        var: What the integral is taken with respect to
    
    Returns:
        The expectation value for the given operator and normalized WaveFunction. 
        An expectation value is the average value of an operator for a given WaveFunction. 

    """
    n = Symbol("n")
    if operator == kinetic_energy(var):
        return sympify(str(sympify(str(Integral(conjugate(wavefunc_1)*operator, (var, lower, upper))).replace(str(Derivative("1", var)**2), str(Derivative(wavefunc_1, var, var)))).doit()).replace(str('sin(pi*n)'), str(0)).replace(str('cos(pi*n)'), str(0)))
    if operator == lin_mom(var):
        return simplify(Integral(conjugate(wavefunc_1)*operator, (var, lower, upper)).replace(Derivative("1", var), Derivative(wavefunc_1, var).doit()))
    else:
        return simplify(Integral(conjugate(wavefunc_1)*operator*wavefunc_2, (var, lower, upper)).doit().replace(sin(n*pi), 0).replace(cos(n*pi), 1))



    
def expectation_value_steps(wavefunc_1, operator, wavefunc_2, lower, upper, var):
    """
    Computes a symbolic expression for an expectation value of an operator `operator`
    with the two wavefunctions: `<WaveFunc_1|operator|WaveFunc_2>` 

    Args:    
        wavefunc_1: The "bra" normalized WaveFunction
        operator: The operator of interest
        wavefunc_2: The "ket" normalized WaveFunction
        lower: The lower bound of the integral. If bounds are not listed, this is -oo
        upper: The upper bound of the integral. If bounds are not listed, this is oo
        var: What the integral is taken with respect to
    
    Returns:
        The steps on how to solve an expectation value problem.
        An expectation value is the average value of an operator for a given WaveFunction. 

    """
    
    n = Symbol("n")
    if operator == kinetic_energy(var):
        return display(Integral(conjugate(wavefunc_1)*operator, (var, lower, upper))), \
                       display((Integral(conjugate(wavefunc_1)*operator, (var, lower, upper))).replace(str(Derivative("1", var)**2), str(Derivative(wavefunc_1, var, var)))).doit(), \
                            display(sympify(str(sympify(str(Integral(conjugate(wavefunc_1)*operator, (var, lower, upper))).replace(str(Derivative("1", var)**2), str(Derivative(wavefunc_1, var, var)))).doit()).replace(str('sin(pi*n)'), str(0)).replace(str('cos(pi*n)'), str(0))))
    if operator == lin_mom(var):
        return display(Integral(conjugate(wavefunc_1)*operator, (var, lower, upper))), \
                       display(Integral(conjugate(wavefunc_1)*operator, (var, lower, upper)).replace(Derivative("1", var), Derivative(wavefunc_1, var).doit())), \
                           display(simplify(Integral(conjugate(wavefunc_1)*operator, (var, lower, upper)).replace(Derivative("1", var), Derivative(wavefunc_1, var).doit())))
    else:
        return display(Integral(conjugate(wavefunc_1)*operator*wavefunc_2, (var, lower, upper))), \
                       display(simplify(Integral(conjugate(wavefunc_1)*operator*wavefunc_2, (var, lower, upper)).doit().replace(sin(n*pi), 0).replace(cos(n*pi), 1)))
    
    
    

def overlap(WaveFunc_1, WaveFunc_2, lower, upper, var):
    """
    Computes a symbolic expression for an overlap of two functions, `<WaveFunc_1|WaveFunc_2>`

    Args:    
        WaveFunc_1: The "bra" normalized WaveFunction
        WaveFunc_2: The "ket" normalized WaveFunction
        var: What the integral is taken with respect to
        lower: The lower bound of the integral. If bounds are not listed, this is -oo
        upper: The upper bound of the integral. If bounds are not listed, this is oo 
        
    Returns:    
        The overlap of the two WaveFunctions of interest over given bounds. 
    
    """
    
    return simplify(Integral(WaveFunc_1*WaveFunc_2, (var, lower, upper)).doit()).replace(sin(n*pi), 0).replace(cos(n*pi), 1)




def plot_function(func, lower1, upper1, B, lower2 = None, upper2 = None, C = None):
    """

    Args:
        func: The function/Normalized WaveFunction of interest
        lower: The lower bound of the B domain (for Particle in a Box, 0)
        upper: The upper bound of the B domain (for Particle in a Box, 1)
        B: This is "x" usually (x-axis) - the domain
        lower2: The lower bound for the C range (for Particle in a Box, 0)
        upper2: The upper bound for the C range (for Particle in a Box, 1)
        C: This is "y" usually (y-axis) - the range

    Returns:
        A plotted function of the function/Normalized WaveFunction of interest.

    Note:
        The cell often has to be run TWICE in order to print the plot/graph.

    """
    
    if C == None and lower2 == None and upper2 == None:
        plt = plot(sympify(str(func)), (B, lower1, upper1))
    else:
        plt = plot3d(sympify(str(func)), (B, lower1, upper1), (C, lower2, upper2))
    return plt




def laguerre(r, n):
    """

    Args:
        r: What the equation is with respect to. This is commonly "r" or "x"
        n: The principle Quantum Number

    Returns:
        The Laguerre polynomial. This is commonly used to solve the laguerre_assoc() function.

    """
    
    return simplify((exp(r)*Derivative(exp(-r)*(r**n), (r, n))).doit())




def laguerre_2(r, n):
    """

    Args:
        r: What the equation is with respect to. This is commonly "r" or "x"
        n: The principle Quantum Number

    Returns:
        The Laguerre polynomial, without simplification and without the derivative computed.

    """
    
    return (exp(r)*Derivative(exp(-r)*(r**n), (r, n)))




def laguerre_assoc(n, l):
    """

    Args:
        n: The principle Quantum Number
        l: The angular momentum Quantum Number

    Returns:
        The Laguerre Associated polynomial, commonly used to solve the radial function.

    """
    
    return simplify(((-1)**l)*Derivative(laguerre_2(r, n+l), (r, l)).doit())



def kronecker(i, j):
    """

    Args:
        i: this is usually "i" and the first variable of the kronecker delta.
        j: this is usually given a numerical value, and is the second variable of the kronecker delta.
    
    Returns:
        This function does not print anything, but instead is used in conjunction with the code found in the kronecker_delta file on:
        https://github.com/liz-stippell/Pysces
    
    """
    if i == j:
        return 1
    else:
        return 0

    
    
    
""" The following is for harmonic oscillators, ladder operators """

def a_raising(A = None):
    """
    
    Args:
        A: this is usually either empty, (no parameter), "normalized" for the normalized ladder operator, or a "symbol" for the symbolic 
        representation of the parameter.
    
    Returns:
        Either the mathematic representation of the "a" raising operator (commonly used in harmonic oscillator problems), the normalized 
        raising operator (in dirac notation), or the symbol notation of the operator.
    
    """
    
    m, omega, p, x, a_R, normalized, n, symbol = symbols("m omega p x a_+ normalized n symbol")
    if A == None:
        return (1/sqrt(2*h_b*m*omega)*(-I*p + m*omega*x))
    if A == normalized:
        return Bra(str(n))*sqrt(n+1)*Ket(str(n+1))
    if A == symbol:
        return a_R


    
        
def a_lowering(A = None):
    """    
    Args:
        A: this is usually either empty, (no parameter), "normalized" for the normalized ladder operator,
        or a "symbol" for the symbolic representation of the parameter.
    
    Returns:    
        Either the mathematic representation of the "a" lowering operator (commonly used in harmonic oscillator problems),
        the normalized lowering operator (in dirac notation), or the symbol notation of the operator.
    
    """
    
    m, omega, p, x, a_L, normalized, n, symbol = symbols("m omega p x a_- normalized n symbol")
    if A == None:
        return (1/sqrt(2*h_b*m*omega)*(I*p + m*omega*x)) 
    if A == normalized:
        return Bra(n)*sqrt(n)*Ket(n-1)
    if A == symbol:
        return a_L



def x_ladder():
    """
    
    Args:
        At this time, there are no parameters necessary.
    
    Returns:
        "x" with respect to the "a" raising/lowering operators
    
    """
    
    m, omega, symbol = symbols("m omega symbol")
    return (sqrt(h_b/(2*m*omega))*(a_raising(symbol)+a_lowering(symbol)))




def simplify_ladder(expr):
    """
    
    Args:
        expr: the expression of interest
    
    Returns:
        The simplified expression, replacing any a_rasing(symbol)**2 with 0, a_lowering(symbol)**2 with 0, and 
        2*a_raising(symbol)*a_lowering(symbol) with (2*n+1). 
    
    Note:
        The code does not "understand" the difference between a_raising(symbol)*a_lowering(symbol) and a_lowering(symbol)*a_raising(symbol) 
        which are mathematically different and produce different outcomes.
    
    """
    symbol, n = symbols("symbol n")
    return Bra(n)*(simplify(expand(simplify(expand(expr)))).replace((2*a_raising(symbol)*a_lowering(symbol)), (2*n+1)).replace(a_raising(symbol)**2, 0).replace(a_lowering(symbol)**2, 0))*Ket(n)




def harmonic(condition = None):
    """    
    Args:
        condition: the condition of the harmonic oscillator. Examples include: "ground_state" or "ladder". 
            This can also be left blank to give the general definition of a harmonic oscillator.

    Returns:
        The Harmonic Oscillator Hamiltonian. Note that "p" is the linear momentum operator, and so the 
        first term is the kinetic_energy() function mentioned above.

    """
    
    p, k, m, x, alpha, omega, ladder, ground_state, symbol = symbols("p k m x alpha omega ladder ground_state symbol")
    if condition == None:
        return (p**2)/(2*m) + (1/2)*k*(x**2)
    if condition == ladder:
        return h_b*omega*(a_lowering(symbol)*a_raising(symbol)-(1/2))
    if condition == ground_state:
        return ((alpha/pi)**(1/4))*exp(-(1/2)*((m*omega)/h_b)*(x**2))
    

""" The following are for Spherical Harmonics """


def spherical(expr):
    """

    Args:
        expr: The expression of interest to be changed into spherical coordinates
    
    Returns:
        The expression of interest, A, in terms of spherical coordinates
    
    """
    
    x, r, theta, phi, y, z = symbols("x r theta phi y z")
    return expr.replace(Derivative(1, x), sin(theta)*cos(phi)*Derivative(1, r) + (1/r)*cos(phi)*cos(theta)*Derivative(1, theta) - (1/r)* \
                        ((sin(phi)/sin(theta)))*Derivative(1, phi)) \
            .replace(Derivative(1, y), sin(theta)*sin(phi)*Derivative(1, r) + (1/r)*sin(phi)*cos(theta)*Derivative(1, theta) + (1/r)* \
                     (cos(phi)/sin(theta))*Derivative(1, phi)) \
            .replace(Derivative(1, z), cos(theta)*Derivative(1, r) - (1/r)*sin(theta)*Derivative(1, theta)) \
            .replace(x, r*sin(theta)*cos(phi)).replace(y, r*sin(theta)*sin(phi)).replace(z, r*cos(theta)) 



def ang_mom_2(ang_mom, mag):
    """
    Args:    
        ang_mom: The total angular momentum quantum number
        mag: The magnetic quantum number
    
    Returns:    
        The L^2 vector magnitude eigenvalue for spherical harmonics.
    
    Note:
        Bra(str(j), str(","), str(m))*j*(j+1)*h_b**2*Ket(str(j), str(","), str(m))
    
    """
    
    j, m = symbols("j, m")
    if ang_mom == j and mag == m:
        return Bra(str(j), str(","), str(m))*j*(j+1)*h_b**2*Ket(str(j), str(","), str(m))
    if (str(ang_mom), str(","), str(mag)) == (str(ang_mom), str(","), str(mag)):
        return ang_mom*(ang_mom+1)*h_b**2
    else:
        return 0 



def ang_mom_z(ang_mom, mag):
    """
    Args:    
        ang_mom: The total angular momentum quantum number
        mag: The magnetic quantum number
    
    Returns:    
        The L_z projection (in the z direction) eigenvalue for spherical harmonics.
    
    Note:
        Bra(str(j), str(","), str(m))*m*h_b*Ket(str(j), str(","), str(m))
    
    """
    
    j, m = symbols("j, m")
    if ang_mom == j and mag == m:
        return Bra(str(j), str(","), str(m))*m*h_b*Ket(str(j), str(","), str(m))
    if (str(ang_mom), str(","), str(mag)) == (str(ang_mom), str(","), str(mag)):
        return mag*h_b
    else:
        return 0 



def ang_mom_raising(ang_mom = None, mag = None):
    """
    Args:    
        ang_mom: The total angular momentum quantum number
        mag: The magnetic quantum number
    
    Returns:    
        If ang_mom == None and mag == None, the general formula for the raising operator for spherical harmonics is returned.    
        Else, the formula for the raising operator is computed using Dirac notation
    
    Note:
        Bra(str(j), str(","), str(m))*h_b*sqrt(j*(j+1)-m*(m+1))*Ket(str(j), str(','), str(m+1))

        .. math::
          \\langle j, m | \\hbar \\sqrt{ j (j+1) - m (m+1) }  | j, m+1 \\rangle
        
    """
    
    L_x, L_y, j, m = symbols("L_x L_y j m")
    if ang_mom == None and mag == None:
        return Operator(L_x) + I*Operator(L_y)
    if ang_mom == j and mag == m:
        return Bra(str(j), str(","), str(m))*h_b*sqrt(j*(j+1)-m*(m+1))*Ket(str(j), str(','), str(m+1))
    if (str(ang_mom), str(","), str(mag)) == (str(ang_mom), str(','), str(mag+1)):
        return h_b*sqrt(ang_mom*(ang_mom+1)-mag*(mag+1))
    else:
        return 0 

    

def ang_mom_lowering(ang_mom = None, mag = None):
    """
    Args:    
        ang_mom: The total angular momentum quantum number
        mag: The magnetic quantum number
    
    Returns:    
        If ang_mom == None and mag== None, the general formula for the lowering operator for spherical harmonics is returned.    
        Else, the formula for the lowering operator is computed using Dirac notation
    
    Note:
        Bra(str(j), str(","), str(m))*h_b*sqrt(j*(j+1)-m*(m-1))*Ket(str(j), str(','), str(m-1))
    
    """
        
    L_x, L_y= symbols("L_x L_y")
    if ang_mom == None and mag == None:
        return Operator(L_x) - I*Operator(L_y)
    if ang_mom == j and mag == m:
        return  Bra(str(j), str(","), str(m))*h_b*sqrt(j*(j+1)-m*(m-1))*Ket(str(j), str(','), str(m-1))
    if (str(ang_mom), str(","), str(mag)) == (str(ang_mom), str(','), str(mag-1)):
        return h_b*sqrt(ang_mom*(ang_mom+1)-mag*(mag-1))
    else: 
        return 0 


    
def ang_mom_x(ang_mom = None, mag = None):
    """
    Args:    
        ang_mom: The total angular momentum quantum number
        mag: The magnetic quantum number
    
    Returns:    
        If ang_mom == None and mag == None, the general formula for the L_x operator for spherical harmonics is returned.    
        Else, the formula for the L_x operator is computed using Dirac notation
    
    """
        
    if ang_mom == None and mag == None:
        L_R, L_L = symbols("L_+ L_-")
        return (1/2)*(L_R + L_L)
    else:
        return (1/2)*(ang_mom_raising(ang_mom, mag) + ang_mom_lowering(ang_mom, mag))
       
    

