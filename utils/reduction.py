from sage.all import  QQ, PolynomialRing


def sub_x_for_u(f, a, Ring =  PolynomialRing(QQ, 'x,y,u,v,dx,dy')):
    """Substitute x^a with u"""
    x, u = Ring.gens()[0], Ring.gens()[2]                                   #x is the first generator, u is the third
    if f ==0:
        return 0
    
    result = 0
    for term in f.monomials():
        coeff = f.monomial_coefficient(term)                                #gets coefficient of the monomial
        xdeg = term.exponents()[0][0]                                       #gets exponent of x in the monomial
        q, r = xdeg//a, xdeg%a                                              #q is the quotient, r is the remainder; xdeg = q*a + r       
        result += Ring(coeff * term * u**q/(x**(xdeg-r)))                   #x^xdeg -> u^q * x^r  
    return result

def sub_y_for_v(f, b, c, Ring =  PolynomialRing(QQ, 'x,y,u,v,dx,dy')):
    """Substitute y^c with x^b-v"""
    x, y, v = Ring.gens()[0], Ring.gens()[1], Ring.gens()[3]                #x is the first generator, y is the second, v is the fourth

    if f ==0:
        return 0
    result = 0
    for term in f.monomials():
        coeff = f.monomial_coefficient(term)                                #gets coefficient of the monomial
        ydeg = term.exponents()[0][1]                                       #gets exponent of y in the monomial
        q, r = ydeg//c, ydeg%c                                              #q is the quotient, r is the remainder; ydeg = q*c + r
        result += Ring(coeff * term * (x**b-v)**(q)/(y**(ydeg-r)))          #y^ydeg -> (x^b-v)^q * y^r
    return result

def reduce_poly(f,a,b,c):
    """Apply substitution rules until achieving fixed point (which happens when xdeg < a and ydeg < c for all monomials in f)."""
    current = f
    new = sub_y_for_v(sub_x_for_u(current, a), b, c)
    while current!=new:
        current = new
        new = sub_y_for_v(sub_x_for_u(new, a), b, c)
    return new
