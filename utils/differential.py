from sage.all import QQ, PolynomialRing
from utils.basis import exp_from_monomial, build_poly_from_exp

def diff_from_exp_to_exp(exp, a, b, c):
    """
    Logic for computing the result of applying the differential (essntially the exteriror derivative) to a monomial in Prism_R or its Nygaard-filtered version.
    Input: exp = [coeff, e1, e2, e3, e4, e5, e6] where e1,e2,e3,e4 are the exponents of x,y,u,v respectively and e5,e6 are the exponents of dx,dy.
    Output: list of exponent representations of the result of applying the differential to the monomial represented by exp.
    """
    coeff,e1,e2,e3,e4,e5,e6 = exp[0],exp[1], exp[2], exp[3], exp[4],exp[5],exp[6]
    sgn = (-1)**(exp[5])
    potential_diffs = []
    if e1>0:
        dif_x = [coeff*e1,e1-1,e2,e3,e4,e5+1,e6]
        potential_diffs.append(dif_x)
    if e2>0:
        dif_y = [coeff*sgn*e2,e1,e2-1,e3,e4,e5,e6+1]
        potential_diffs.append(dif_y)
    if e3>0:
        dif_u = [coeff*a,e1+a-1,e2,e3-1,e4,e5+1,e6]
        potential_diffs.append(dif_u)
    if e4>0:
        dif_vx = [coeff*b,e1+b-1,e2,e3,e4-1,e5+1,e6]
        potential_diffs.append(dif_vx)

        dif_vy = [coeff*(-c)*sgn,e1,e2+c-1,e3,e4-1,e5,e6+1]
        potential_diffs.append(dif_vy)

    result = []
    for entry in potential_diffs:
        if entry[0] != 0 and entry[5]<2 and entry[6]<2: 
            result.append(entry)
    return result

def diff_function(f,a,b,c, Ring = PolynomialRing(QQ, 'x,y,u,v,dx,dy')):
    """Compute differential of polynomial in Prism_R or Nygaard-filtered version using diff_from_exp_to_exp and reduction."""
    if f ==0:
        return 0
    result = 0
    for term in f.monomials():
        term_as_exp = exp_from_monomial(term*f.monomial_coefficient(term))
        result += sum([build_poly_from_exp(entry, Ring) for entry in diff_from_exp_to_exp(term_as_exp,a,b,c)]) 
    return result
