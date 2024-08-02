

import time #just to timestamp computations
import copy as copier #to avoid clash with the sagemath copy method
import numpy

from sage.all import *
#from IPython import get_ipython


#input: (p,i,a,b,c);
#output: dictionary with entries 'time' (time of computation), 'Hn' (homology in deg n=0,1,2), ??

def syntomic_package(prime,b,c,a,ii):
    
    #Only consider b>a, otherwise the ring is just Z_p[x,y]/(x**a, y**c) which can be derived from Sulyma's as tensor.
    if b>a:
        return syntomic_package(prime,a,c,a,ii)
    
    start_time = time.time()
    N=a*c*ii+(b-1)*(c-1)-1
    
    #incorporating dx,dy as part of the tuples/exponents
    #input a,b,c, where we are looking at ring R = Z_p[x,y]/(x**a, x**b-y**c)
    #https://doc.sagemath.org/html/en/reference/padics/sage/rings/padics/factory.html
    R = PolynomialRing(QQ, 'x,y,u,v,dx,dy')
    x,y,u,v,dx,dy = R.gen(0), R.gen(1), R.gen(2), R.gen(3), R.gen(4), R.gen(5)
    S = Zp(prime, prec = 30, print_mode= 'terse', show_prec = False) 
    

    #write u = x**a, v=x**b-y**c.
    #our basis is {1,x,..,x**{a-1}}\otimes{1,y,...,y**{c-1}} \otimes {1,u,u**[2],...}\otimes {v,v**[2],...}
    #** instead of ^ throughout for being able to run as .py

    def poly_from_exp(exp): #takes a list as input; exp[0] is really a coefficient, not an exponent
        return exp[0]*x**(exp[1])*y**(exp[2])*(u**(exp[3])/factorial(exp[3]))*(v**(exp[4])/factorial(exp[4]))*(dx**(exp[5]))*(dy**(exp[6]))

    def exp_from_base(monomial): #takes a monomial as input.
        if len(monomial.monomials())>1:
            raise ValueError("Input is not a monomial, so this procedure won't work.")
        e3 = monomial.exponents()[0][2]
        e4 = monomial.exponents()[0][3]
        return [monomial.coefficients()[0]*factorial(e3)*factorial(e4)] + list(monomial.exponents()[0]) 


    exps_0 = {}
    basis_0 = {}

    #the following basically solves the equation dd = e_1c+e_2b+e_3ac+e_4bc
    for dd in range(N+1):
        eq_dd = [(-dd, c, b, a*c, b*c)]                             #exponents s.t. total weight of monomials is dd
        ieq_pos = [(0,1,0,0,0),(0,0,1,0,0),(0,0,0,1,0),(0,0,0,0,1)] #all exponents must be positive
        ieq_nil =[(a-1,-1,0,0,0),(c-1,0,-1,0,0)]                    #exponents of x,y are bounded by the ideal equations
        ieq = ieq_pos + ieq_nil                                     #all inequality constraints combined
        convex_poly = Polyhedron(eqns= eq_dd, ieqs = ieq , base_ring=QQ)  #LP problem computing basis for wt dd part
        exps_0[dd] = convex_poly.integral_points()
        exps_0[dd] = [[1]+list(exp)+[0,0] for exp in exps_0[dd]]
        basis_0[dd] = [poly_from_exp(exp) for exp in exps_0[dd]]

    for dd in range(-N-1,0):    #just ensures that computing exps_1, exps_2 returns no error when index out of range
        exps_0[dd] = []

    #use deg 0 part of chain complex to compute deg 1, 2 parts of chain cx 
    exps_1 = {dd: [] for dd in range(N+1+b+c)}
    exps_2 = {dd: [] for dd in range(N+1+b+c)}
    for dd in range(N+1):
        for exp in exps_0[dd]:
            new_dx, new_dy, new_dxdy = copier.deepcopy(exp), copier.deepcopy(exp), copier.deepcopy(exp)
            new_dx[5] = 1
            new_dy[6] = 1
            new_dxdy[5] = 1
            new_dxdy[6] = 1
            exps_1[dd+c].append(new_dx)
            exps_1[dd+b].append(new_dy)
            exps_2[dd+b+c].append(new_dxdy)

    basis_1 = {dd: [poly_from_exp(exp) for exp in exps_1[dd]] for dd in range(N+1+b+c)}  #dx, dy bases of deg 1, wt dd part of LWOmega
    basis_2 = {dd: [poly_from_exp(exp) for exp in exps_2[dd]] for dd in range(N+1+b+c)}  #dx*dy base of deg 2, wt dd part of LWOmega (NOTE dx,dy DO NOT COMMUTE!!!)


    #computing differentials without a reduction algorithm is ~not great~
    #reduction algo (ie how to express a given poly as a linear combination of basis elements)
    #basically does the next two subsitutions until no longer possible: x**N = u*x**(N-a), y**M = y**(M-c)*(x**b-v)
    #so we really get x**N -> x**(N%a)*u**(N//a), y**M -> y**(M%a)*(x**b-v)**(M//c)
    #takes in polynomial in R and returns a polynomial
    def sub_x_for_u(f):                   #input and output are polys in R
        if f ==0:
            return 0
        result = 0
        for term in f.monomials():    #CAREFUL .monomials() forgets coefficients
            xdeg = term.exponents()[0][0] #gets exponent of x in the monomial
            q = xdeg//a
            r = xdeg%a
            result += R(term*u**q/(x**(xdeg-r)) * f.monomial_coefficient(term))
        return result


    def sub_y_for_v(f):                   #input and output are polys in R
        if f ==0:
            return 0
        result = 0
        for term in f.monomials():
            ydeg = term.exponents()[0][1] #gets exponent of y in the monomial
            q = ydeg//c
            r = ydeg%c
            #(q_poly, r_poly) = term.quo_rem(y**(ydeg-r))  #polynomial division of term by y**(ydeg-r); r should be zero
            #if r_poly!=0:
                #print("error subbing y")
            #result += q_poly*(x**b-v)**(q)
            result += R(term*(x**b-v)**(q)/(y**(ydeg-r)) * f.monomial_coefficient(term) )
        return result

    def reduct(poly):  #input: poly f, output: same poly, but having substituted powers of x,y for u,v as much as ppossible
        current = poly
        new = sub_y_for_v(sub_x_for_u(current))
        while current!=new:
            current = new
            new = sub_y_for_v(sub_x_for_u(new))
        return new




    #function computing the differential in terms of the exponent representation, and to the exp rep
    #should handle any monomial in x,y,u,v,dx,dy
    #input is the exp rep of a basis elt; spits out a list of exp reps of the basis elts corresponding to the summands of the image
    def diff_from_exp_to_exp(exp):  
        e0,e1,e2,e3,e4,e5,e6 = exp[0],exp[1], exp[2], exp[3], exp[4],exp[5],exp[6]
        sgn = (-1)**(exp[5])
        potential_diffs = []
        if e1>0:
            dif_x = [e0*e1,e1-1,e2,e3,e4,e5+1,e6]
            potential_diffs.append(dif_x)
        if e2>0:
            dif_y = [e0*sgn*e2,e1,e2-1,e3,e4,e5,e6+1]
            potential_diffs.append(dif_y)
        if e3>0:
            dif_u = [e0*a,e1+a-1,e2,e3-1,e4,e5+1,e6]
            potential_diffs.append(dif_u)
        if e4>0:
            dif_vx = [e0*b,e1+b-1,e2,e3,e4-1,e5+1,e6]
            potential_diffs.append(dif_vx)

            dif_vy = [e0*(-c)*sgn,e1,e2+c-1,e3,e4-1,e5,e6+1]
            potential_diffs.append(dif_vy)
        result = []
        for entry in potential_diffs:
            if entry[0] != 0 and entry[5]<2 and entry[6]<2: 
                result.append(entry)
        return result

    def diff_from_exp_to_poly(exp): #input an exp, takes differential, returns differential as a poly
        return sum([poly_from_exp(entry) for entry in diff_from_exp_to_exp(exp)])

    def diff_from_exp_to_basis(exp):
        p = diff_from_exp_to_poly(exp)
        return reduct(p)      #using .monomials() would forget coefficients

    #dictionary of bases for the crystalline/prismatic cohomology complexes in deg=0,1,2; F-degree dd; ie Delta^deg_dd
    prismatic_exps = {0:{}, 1:{}, 2:{}}
    for dd in range(N+1):
        prismatic_exps[0][dd] = exps_0[dd]
        prismatic_exps[1][dd] = exps_1[dd]
        prismatic_exps[2][dd] = exps_2[dd]

    prismatic_bases = {0:{}, 1:{}, 2:{}}
    for dd in range(N+1):
        prismatic_bases[0][dd] = basis_0[dd]
        prismatic_bases[1][dd] = basis_1[dd]
        prismatic_bases[2][dd] = basis_2[dd]

    #dictionary of transition matrices for the crystalline/prismatic cohomology complexes in deg=0,1; F-degree dd
    #computes up to degree N

    prismatic_matrices = {0:{}, 1:{}}
    for deg in [0,1]:
        for dd in range(N+1):
            M=Matrix(ZZ,len(prismatic_bases[deg+1][dd]),len(prismatic_exps[deg][dd]));M

            for e in range(len(prismatic_exps[deg][dd])):
                #next lines until def of diffs compute differential and split it into a list of bases monomials with corresponding coeffs 
                #print(prismatic_exps[deg][dd][e])
                diff_as_poly = R(diff_from_exp_to_basis(prismatic_exps[deg][dd][e]))
                #print(diff_as_poly)
                diffs= [term*diff_as_poly.monomial_coefficient(term) for term in diff_as_poly.monomials()] 
                for t in range(len(diffs)):
                    #print(prismatic_exps[deg][dd][e])
                    #print(prismatic_bases[deg+1][dd])
                    #print(deg,dd,e,t, diffs)
                    prismatic_bases_without_coeff = [term.monomials()[0] for term in prismatic_bases[deg+1][dd]]
                    i=prismatic_bases_without_coeff.index(diffs[t].monomials()[0])
                    coeff=diffs[t].coefficients()[0] / prismatic_bases[deg+1][dd][i].coefficients()[0]
                    M[i,e]+=coeff
                    #print(coeff)
            prismatic_matrices[deg][dd] = copier.deepcopy(M)


    #method to print a cryst/prism complex: [Delta**0_dd, matrix, Delta**1_dd, matrix, Delta**2_dd]
    def prismatic_complex_in_weight(dd):
        return(prismatic_bases[0][dd], prismatic_matrices[0][dd],prismatic_bases[1][dd],prismatic_matrices[1][dd],prismatic_bases[2][dd])

    #double parentheses function from Sulyma
    def nonneg(a):
        return max(a,0)

    #dictionary/list of bases for the Nygaard-filtered complexes in filtration ii; deg=0,1,2; weight dd


    Nyg_exps = {0:{}, 1:{}, 2:{}} #the to-be three terms of the ii-th filtered part of the Nyg chain complex
    for dd in range(N+1):    
        Nyg_exps[0][dd] = [(e[0]*prime**(nonneg(ii-e[3]-e[4])),e[1],e[2],e[3],e[4],e[5],e[6]) for e in exps_0[dd]]
        Nyg_exps[1][dd] = [(e0*prime**(nonneg(ii-e3-e4-1)),e1,e2,e3,e4,e5,e6) for (e0,e1,e2,e3,e4,e5,e6) in exps_1[dd]]
        Nyg_exps[2][dd] = [(e0*prime**(nonneg(ii-e3-e4-2)),e1,e2,e3,e4,e5,e6) for (e0,e1,e2,e3,e4,e5,e6) in exps_2[dd]]
        #REVISE EXPONENT OF p IN THE LINES ABOVE. I think it's good now


    Nyg_bases = {0:{}, 1:{}, 2:{}}
    for dd in range(N+1): 
        Nyg_bases[0][dd] = [poly_from_exp(term) for term in Nyg_exps[0][dd]]
        Nyg_bases[1][dd] = [poly_from_exp(term) for term in Nyg_exps[1][dd]]
        Nyg_bases[2][dd] = [poly_from_exp(term) for term in Nyg_exps[2][dd]]        



    #dictionary/list of transition matrices for the Nygaard filtered complexes in deg=0,1; F-degree dd; filtration ii

    Nyg_matrices = {0:{}, 1:{}}
    for deg in [0,1]:
        for dd in range(N+1):
            M=Matrix(ZZ,len(Nyg_bases[deg+1][dd]),len(Nyg_exps[deg][dd]));M #OK to define over ZZ rather than QQ?

            for e in range(len(Nyg_exps[deg][dd])):
                #next lines until def of diffs compute differential and split it into a list of bases monomials with corresponding coeffs 
                #print(prismatic_exps[deg][dd][e])
                diff_as_poly = R(diff_from_exp_to_basis(Nyg_exps[deg][dd][e])) #CHECK WORKS IF COEFFE != 1
                #print(diff_as_poly)
                diffs= [term*diff_as_poly.monomial_coefficient(term) for term in diff_as_poly.monomials()] 
                for t in range(len(diffs)):
                    #print(prismatic_exps[deg][dd][e])
                    #print(prismatic_bases[deg+1][dd])
                    #print(deg,dd,e,t, diffs)
                    Nyg_bases_without_coeff = [term.monomials()[0] for term in Nyg_bases[deg+1][dd]]
                    i=Nyg_bases_without_coeff.index(diffs[t].monomials()[0])
                    coeff=diffs[t].coefficients()[0] / Nyg_bases[deg+1][dd][i].coefficients()[0]
                    M[i,e]+=coeff
                    #print(coeff)
            Nyg_matrices[deg][dd] = copier.deepcopy(M)

    #---------------------


    #method to print a Nygaard-filtered complex: [N^ii,0_dd; matrix; N^ii,1_dd; matrix; N^ii,2_dd]
    def Nyg_ii_complex_in_weight(dd):
        return(Nyg_bases[0][dd], Nyg_matrices[0][dd],Nyg_bases[1][dd],Nyg_matrices[1][dd],Nyg_bases[2][dd])

    #method to compute phi_i = phi/p^i; input: poly in R[x,y,u,v,dx,dy]; output: reduced poly in R[x,y,u,v,dx,dy]
    #POTENTIAL SPEEDUP: we only use the 'truncated Frobenius', so could kill terms of degree>ideal_truncation before reducing
    def Frob_ii(f):
        if f ==0:
            return 0
        result = 0
        for term in f.monomials():    #CAREFUL .monomials() forgets coefficients
            xdeg = term.exponents()[0][0] #gets exponent of x in the monomial
            ydeg = term.exponents()[0][1]
            udeg = term.exponents()[0][2]
            vdeg = term.exponents()[0][3]
            dxdeg = term.exponents()[0][4]
            dydeg = term.exponents()[0][5]
            xfrob = x**(prime*xdeg)
            yfrob = y**(prime*ydeg)
            ufrob = u**(prime*udeg)
            vfrob = (x**(b*prime) - y**(c*prime))**(vdeg)
            dxfrob = (dx*prime*x**(prime-1))**dxdeg #the ^dxdeg just ensures it returns 1 if dxdeg=0; and the correct Frob if dxdeg=1
            dyfrob = (dy*prime*y**(prime-1))**dydeg #idem for ^dydeg
            result += R(f.monomial_coefficient(term)/(prime**ii) *xfrob*yfrob*ufrob*vfrob*dxfrob*dyfrob)
        return reduct(result)

    #better method of computing the differential: input: poly in R[x,y,u,v,du,dv]; output: reduced poly in idem
    def diff_function(f):
        if f ==0:
            return 0
        result = 0

        for term in f.monomials():    #CAREFUL .monomials() forgets coefficients
            result += diff_from_exp_to_poly( exp_from_base( term*f.monomial_coefficient(term) ))
        return reduct(result)


    #TODO: could write a function matrix_of_f(basis_1, basis_2, func) that takes a linear function going from basis_1 to basis_2
    #and spits out the associated matrix (where all these are polys/diff forms in Q_p and the matrix entries are in Z_p)
    #could use that to replace some code in computing prismatic, Nygaard, and Frob matrices!

    def degree_of(f):
        if f ==0:
            return 0
        if len(f.monomials())>1:
            raise ValueError("Input is not a monomial, so this procedure won't work.")
        (e1,e2,e3,e4,e5,e6) = f.exponents()[0]
        return e1*c + e2*b + e3*a*c + e4*b*c + e5*c + e6*b

    def matrix_of(basis_1, basis_2, func):

        M=Matrix(ZZ,len(basis_2),len(basis_1));M

        for e in range(len(basis_1)):
            #next lines until def of diffs compute differential and split it into a list of bases monomials with corresponding coeffs 
            #print(prismatic_exps[deg][dd][e])
            func_as_poly = R(func(basis_1[e]))
            #print(diff_as_poly)
            func_terms= [term*func_as_poly.monomial_coefficient(term) for term in func_as_poly.monomials()] 
            for t in range(len(func_terms)):
                #print(prismatic_exps[deg][dd][e])
                #print(prismatic_bases[deg+1][dd])
                #print(deg,dd,e,t, diffs)
                basis2_without_coeff = [term.monomials()[0] for term in basis_2]
                i=basis2_without_coeff.index(func_terms[t].monomials()[0])
                coeff=func_terms[t].coefficients()[0] / basis_2[i].coefficients()[0]
                M[i,e]+=coeff
                #print(coeff)        
        return M

    #helper function; input: elementary divisors list; output: list of nonzero v_prime valuations
    def torsion_exps_from_elem_divs(ed): 
        result = []
        for term in ed:
            if term ==0 or term == 1:
                pass
            else:
                result.append(valuation(term,prime))
        return result

    
    #trunc=N+1 is the ideal cutoff; 

    #-----------------------

    #now that we have ideal cutoff, recompute matrices one last time
    #list of numbers coprime to p and less than N+1
    I_p = [n for n in range(N+1) if gcd(prime, n)==1]
    rays = {}
    for n in I_p:
        k=0
        rays[n] = []
        while n*prime**k < N+1:
            rays[n].append(n*prime**k)
            k+=1

    #terms in the total chain complex: N^0 -> N^1+Prism^0 -> N^2+Prism^1 -> Prism^2; whose cohom is the syntomic cohom
    #ii is Nyg filtered degree; deg =0,1,2,3; dd\leq trunc
    total_complex_bases = {0:{}, 1:{}, 2:{}, 3:{}}
    total_Nyg_bases = {0:{}, 1:{}, 2:{}}
    total_prism_bases = {0:{}, 1:{}, 2:{}}

    for n in I_p:
        for i in range(3):
            total_Nyg_bases[i][n] = []
            total_prism_bases[i][n] = []

        for dd in rays[n]:
            #N^0, N^1, N^2
            total_Nyg_bases[0][n] += Nyg_bases[0][dd]
            total_Nyg_bases[1][n] += Nyg_bases[1][dd] 
            total_Nyg_bases[2][n] += Nyg_bases[2][dd]

            #P^1,P^2,P^3
            total_prism_bases[0][n] += prismatic_bases[0][dd]
            total_prism_bases[1][n] += prismatic_bases[1][dd]
            total_prism_bases[2][n] += prismatic_bases[2][dd]

        #N^0, N^1+P^0, N^2+P^1, P^2
        total_complex_bases[0][n] = total_Nyg_bases[0][n]
        total_complex_bases[1][n] = total_Nyg_bases[1][n] + total_prism_bases[0][n]
        total_complex_bases[2][n] = total_Nyg_bases[2][n] + total_prism_bases[1][n]
        total_complex_bases[3][n] = total_prism_bases[2][n]    



    def cap_deg(f): #helper function for when image of frob_ii is above the degree we are computing, it just forces it to be zero
        if degree_of(f.monomials()[0]) >= N+1:
            return 0
        else:
            return f

    total_complex_matrices = {0:{}, 1:{}, 2:{}}    

    for n in I_p:
        for i in range(3):
            total_complex_matrices[i][n] = None

        #TODO: this seems what takes longest, not the SNF procedures in next  cell
        #transition functions in total complex chain complex
        #from N^0->P^0+N^1
        f0 = lambda t : diff_function(t) + t - cap_deg(Frob_ii(t))    #this is d+id-frob_ii
        #print(ii, total_complex_bases0], total_complex_bases[1])
        #from P^0+N^1 -> P^1+N^2
        def f1(t):
            if t in total_Nyg_bases[1][n]:
                #d+id-frob_ii on the Nygaard summand
                return diff_function(t) + t - cap_deg(Frob_ii(t))    
            else:
                #-d on the prism summand
                return -diff_function(t)                  
        #from P^1+N^2 -> P^2
        def f2(t):
            if t in total_Nyg_bases[2][n]:
                #d+id-frob_ii on the Nygaard summand
                return t - cap_deg(Frob_ii(t))
            else:
                #-d on the prism summand
                return -diff_function(t)

        total_complex_matrices[0][n] = matrix_of(total_complex_bases[0][n], total_complex_bases[1][n], f0)
        total_complex_matrices[1][n] = matrix_of(total_complex_bases[1][n], total_complex_bases[2][n], f1)
        total_complex_matrices[2][n] = matrix_of(total_complex_bases[2][n], total_complex_bases[3][n], f2)
        
      #----------------

    #We have the transition matrices of total complex 
    #the upcoming processes them to give the syntomic cohomology data: 
    #1) computes their Smith Normal Form; or rather their elementary divisors
    #2) uses this to give the homology
    #https://www.matem.unam.mx/~omar/mathX27/smith-form.html is a concise reference
    #we compute over ZZ, then throw out anything coprime to p=prime

    elem_divisors = {0:{}, 1:{}, 2:{}}
    dims_bases = {0:{}, 1:{}, 2:{}, 3:{}}
    ranks = {0:{}, 1:{}, 2:{}}
    homology_dims = {0:{}, 1:{}, 2:{}, 3:{}}
    torsion_exps = {1:{}, 2:{}, 3:{}}

    for n in I_p:
        #elementary divisors of each transition matrix
        elem_divisors[0][n] = total_complex_matrices[0][n].elementary_divisors()
        elem_divisors[1][n] = total_complex_matrices[1][n].elementary_divisors()
        elem_divisors[2][n] = total_complex_matrices[2][n].elementary_divisors()

        #dimensions of bases of the terms in total complex
        dims_bases[0][n] = len(total_complex_bases[0][n])
        dims_bases[1][n] = len(total_complex_bases[1][n])
        dims_bases[2][n] = len(total_complex_bases[2][n])
        dims_bases[3][n] = len(total_complex_bases[3][n])

        #ranks of the transition matrices (r in the link above)
        ranks[0][n] = numpy.count_nonzero(elem_divisors[0][n])
        ranks[1][n] = numpy.count_nonzero(elem_divisors[1][n])
        ranks[2][n] = numpy.count_nonzero(elem_divisors[2][n])

        #ranks of H^0, H^1, H^2 of the synotmic complex (m-r-s in link above)
        homology_dims[0][n] = dims_bases[0][n] - ranks[0][n]
        homology_dims[1][n] = dims_bases[1][n] - ranks[1][n] - ranks[0][n]
        homology_dims[2][n] = dims_bases[2][n] - ranks[2][n] - ranks[1][n]
        homology_dims[3][n] = dims_bases[3][n]               - ranks[2][n] 

        #torsion summands (a_i in the link above)
        torsion_exps[1][n] = torsion_exps_from_elem_divs(elem_divisors[0][n])
        torsion_exps[2][n] = torsion_exps_from_elem_divs(elem_divisors[1][n])
        torsion_exps[3][n] = torsion_exps_from_elem_divs(elem_divisors[2][n])

    #for nn in I_p:
        #print("complex for j="+str(nn))
        #print(total_complex_bases[0][nn], total_complex_matrices[0][nn], total_complex_bases[1][nn], total_complex_matrices[1][nn], total_complex_bases[2][nn], total_complex_matrices[2][nn], total_complex_bases[3][nn]) 
        #print(torsion_exps[1][nn], torsion_exps[2][nn], torsion_exps[3][nn])
        
        
    #combine the info from different rays
    elem_divs_0 = []
    elem_divs_1 = []
    elem_divs_2 = []
    dims_bases_0 = 0
    dims_bases_1 = 0
    dims_bases_2 = 0
    dims_bases_3 = 0
    rank_0 = 0
    rank_1 = 0
    rank_2 = 0
    homology_dim_0 = 0
    homology_dim_1 = 0
    homology_dim_2 = 0
    homology_dim_3 = 0
    torsion_1 = []
    torsion_2 = []
    torsion_3 = []

    for n in I_p:
        elem_divs_0 += elem_divisors[0][n]
        elem_divs_1 += elem_divisors[1][n]
        elem_divs_2 += elem_divisors[2][n]
        dims_bases_0 += dims_bases[0][n]
        dims_bases_1 += dims_bases[1][n]
        dims_bases_2 += dims_bases[2][n]
        dims_bases_3 += dims_bases[3][n]
        rank_0 += ranks[0][n]
        rank_1 += ranks[1][n]
        rank_2 += ranks[2][n]
        #homology_dims should be zero; they are the ranks of the homology (number of free generators)
        homology_dim_0 += homology_dims[0][n]
        homology_dim_1 += homology_dims[1][n]
        homology_dim_2 += homology_dims[2][n]
        homology_dim_3 += homology_dims[3][n]
        #this is the interesting bit; it seems torsion_3 always vanishes?
        torsion_1 += torsion_exps[1][n]
        torsion_2 += torsion_exps[2][n]
        torsion_3 += torsion_exps[3][n]
        
    #sort the torsion exponents for readability
    torsion_1.sort()
    torsion_2.sort()
    torsion_3.sort()
        
    computation_time = time.time()-start_time
    
    synt_data = {'computed': True, 'time':computation_time, 'p':prime, 'i':ii, 'a':a, 'b':b, 'c':c, 'input_N': N, 'trunc':N+1, 'length':a*c, 
                'rkH_0':homology_dim_0, 'rkH_1':homology_dim_1, 'rkH_2':homology_dim_2, 'rkH_3':homology_dim_3,
                'H1':torsion_1, 'H2':torsion_2, 'H3':torsion_3}
    
    #raise warning if we ever get non-torsion stuff (we shouldn't!) or H3 doesn't vanish (conjecture it always does)
    if homology_dim_0+homology_dim_1+homology_dim_2+homology_dim_3 >0 or len(torsion_3)>0:
        with open('log.py', 'w') as f:
            f.write(str(synt_data))
            f.close()
        raise ValueError("look at log!")
    
    return synt_data

