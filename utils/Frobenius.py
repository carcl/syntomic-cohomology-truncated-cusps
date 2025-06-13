from sage.all import QQ, PolynomialRing

def Frob_i(f, p, i, a, b, c, Ring = PolynomialRing(QQ, 'x,y,u,v,dx,dy')):
    """Truncated Frobenius divided by p^i, followed by reduction."""
    if f ==0:
            return 0
    x, y, u, v, dx, dy = Ring.gens()
    result = 0
    for term in f.monomials():    
        coeff = f.monomial_coefficient(term)
        (e1, e2, e3, e4, e5, e6) = term.exponents()[0]
        xfrob = x**(p*e1)
        yfrob = y**(p*e2)
        ufrob = u**(p*e3)
        vfrob = (x**(b*p) - y**(c*p))**(e4)
        dxfrob = (dx*p*x**(p-1))**e5 #the ^dxdeg just ensures it returns 1 if dxdeg=0; and the correct Frob if dxdeg=1
        dyfrob = (dy*p*y**(p-1))**e6 #idem for ^dydeg
        result += Ring(coeff/(p**i) *xfrob*yfrob*ufrob*vfrob*dxfrob*dyfrob)
    return result

def degree_of(f,a,b,c):
    """Computes the degree of a monomial in terms of a, b, c."""
    if f ==0:
        return 0
    if len(f.monomials())>1:
        raise ValueError("Input is not a monomial, so this procedure won't work.")
    (e1,e2,e3,e4,e5,e6) = f.exponents()[0]
    return e1*c + e2*b + e3*a*c + e4*b*c + e5*c + e6*b

def cap_deg(f,a,b,c,N): 
    """Function to get rid of terms above F-degree N."""
    result = 0
    for term in f.monomials():
        if degree_of(term, a, b, c) <= N:
            result += term * f.monomial_coefficient(term)
    return result



def faster_Frob(f, p, i, a, b, c, N, Ring = PolynomialRing(QQ, 'x,y,u,v,dx,dy')):
    """Efficient Frobenius: avoids computing terms above F-degree N."""
    if f == 0:
        return 0
    x, y, u, v, dx, dy = Ring.gens()
    result = 0
    for term in f.monomials():

        # Compute new F-degree under Frobenius map; check whether it exceeds N
        # Frobenius map scales the degree of a monomial by p
        # eg v^e4 = (x^b-y^c)^e4 -> (x^{pb}-y^{pc})^e4 is homogeneous of degree p*b*c*e4 as |x^b| = |y^c| = b*c
        new_deg = p*degree_of(term, a, b, c) 
        if new_deg > N:
            continue

        # If we reach this point, the monomial contributes to the result
        coeff = f.monomial_coefficient(term)
        (e1, e2, e3, e4, e5, e6) = term.exponents()[0]
        xfrob = x**(p*e1)
        yfrob = y**(p*e2)
        ufrob = u**(p*e3)
        vfrob = (x**(b*p) - y**(c*p))**e4
        dxfrob = (dx * p * x**(p-1))**e5
        dyfrob = (dy * p * y**(p-1))**e6

        monomial = xfrob * yfrob * ufrob * vfrob * dxfrob * dyfrob
        result += Ring(coeff / (p**i)) * monomial

    return result