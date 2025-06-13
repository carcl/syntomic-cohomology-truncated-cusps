from sage.all import factorial, Polyhedron, QQ, PolynomialRing


def build_poly_from_exp(exp, Ring =  PolynomialRing(QQ, 'x,y,u,v,dx,dy')):
    """Builds a polynomial from a list exponent representation."""
    x, y, u, v, dx, dy = Ring.gens()
    coef, e1, e2, e3, e4, e5, e6 = exp
    return coef * x**e1 * y**e2 * (u**e3 / factorial(e3)) * (v**e4 / factorial(e4)) * dx**e5 * dy**e6

def exp_from_monomial(f):
    """Extracts exponent representation from monomial."""
    if len(f.monomials()) > 1:
        raise ValueError("Not a monomial")
    e3 = f.exponents()[0][2]
    e4 = f.exponents()[0][3]
    return [f.coefficients()[0] * factorial(e3) * factorial(e4)] + list(f.exponents()[0])

def compute_basis_exps(a, b, c, N):
    """Builds exponent dictionaries for the bases in 0,1,2-degree."""
    exps_0 = {}

    #nonnegative integer solutions of the equation wt = e_1c+e_2b+e_3ac+e_4bc for 1<=wt<=N
    for wt in range(N+1):
        eq_wt = [(-wt, c, b, a*c, b*c)]                                     #selects for (e1,e2,e3,e4) s.t. F-weight of corresponding monomial = wt
        ieq_1 = [(0,1,0,0,0),(0,0,1,0,0),(0,0,0,1,0),(0,0,0,0,1)]           #e1,e2,e3,e4 are non-negative integers
        ieq_2 =[(a-1,-1,0,0,0),(c-1,0,-1,0,0)]                              #e1<a and e2<c by choice of basis (x^a=u; y^c=v)
        ieq = ieq_1 + ieq_2                                                 #all inequality constraints combined
        convex_poly = Polyhedron(eqns= eq_wt, ieqs = ieq , base_ring=QQ)    #LP problem computing basis for F^{wt}\Prism_R^0
        exps_0[wt] = convex_poly.integral_points()
        exps_0[wt] = [[1]+list(exp)+[0,0] for exp in exps_0[wt]]            #[e1, e2, e3, e4] -> [coef=1, e1, e2, e3, e4, e5=0, e6=0]

    #use deg 0 part of chain complex to compute deg 1, 2 parts of chain cx (ie from \Prism_R^0 to \Prism_R^1 and \Prism_R^2)
    exps_1 = {wt: [] for wt in range(1, N+1)}
    exps_2 = {wt: [] for wt in range(1, N+1)}
    for wt in range(N+1):
        for exp in exps_0[wt]:
            #f -> f.dx, where f is a monomial in the basis of \Prism_R^0
            if wt+c < N+1:
                new_dx = exp[:]
                new_dx[5] = 1
                exps_1[wt+c].append(new_dx)
            #f -> f.dy
            if wt+b < N+1:
                new_dy = exp[:]
                new_dy[6] = 1
                exps_1[wt+b].append(new_dy)
            #f -> f.dx.dy
            if wt+b+c < N+1:
                new_dxdy = exp[:]
                new_dxdy[5], new_dxdy[6] = 1, 1
                exps_2[wt+b+c].append(new_dxdy)
    
    return exps_0, exps_1, exps_2
