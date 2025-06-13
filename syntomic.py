import time
import argparse
from sage.all import PolynomialRing, QQ, gcd, is_prime
from IPython.display import display
from utils.basis import build_poly_from_exp, compute_basis_exps
from utils.reduction import reduce_poly
from utils.differential import diff_function
from utils.matrix_and_cohom import matrix_of, p_power_count, combine_power_counts, format_group_str
from utils.Frobenius import Frob_i, cap_deg, faster_Frob


def syntomic_complex(p, i, a, b, c, verbose=False, fast_frob=True):
    """
    Main function encapsulating the central algorithm described in the dissertation. 
    Computes a cochain complex model of Z_p(i)(R) for a given prime p, syntomic weight i, and R=F_p[x,y]/(x^a, x^b-y^c).
    The complex is constructed as the homotopy fiber of the map can-Frob_i: N^* -> P^*.
    Here, N^* is a concrete model for the i-th Nygaard-filtered part of the prismatic complex and P^* for the prismatic complex.

    Parameters:
    - p: prime number
    - i: syntomic weight (non-negative integer)
    - a, b, c: positive integers
    - verbose: if True, prints detailed information about the computation
    - fast_frob: if True, uses a faster version of the Frobenius map that discards terms of sufficiently high degree during the computation.

    Returns:
    
    - I_p: a list of integers coprime to p and less than N+1, which index the direct summands of the total complex
    - syntomic_complex_bases: a dictionary, indexed by F-degree and j, containing bases for the modules of the total complex equivalent to Z_p(i)(R)
    - syntomic_complex_matrices: a dictionary, indexed by F-degree and j, containing the transition matrices of the maps in the total complex
    """
    if b > a:
        return syntomic_complex(p, i, a, c, a)

    # N is the maximum F-degree of the polynomials in the bases to be computed
    # Crucially, F^{>N}Z_p(i)(R) = 0, as explained in the thesis
    N = a * c * i + (b - 1) * (c - 1) - 1

    # Polynomial Ring and generators
    base_ring = PolynomialRing(QQ, 'x,y,u,v,dx,dy')
    x, y, u, v, dx, dy = base_ring.gens()

    """Step 1: initialize linear bases for the F-graded prismatic complex (F^w\\Prism_R^* with 1<= w <= N; *=0,1,2)"""
    # See basis.py
    exps_0, exps_1, exps_2 = compute_basis_exps(a, b, c, N)
    prismatic_bases = {0:{}, 1:{}, 2:{}}
    prismatic_bases[0] = {wt: [build_poly_from_exp(exp, base_ring) for exp in exps_0[wt]] for wt in range(1, N+1)}  #0-forms
    prismatic_bases[1] = {wt: [build_poly_from_exp(exp, base_ring) for exp in exps_1[wt]] for wt in range(1, N+1)}  #1-forms (dx, dy)
    prismatic_bases[2] = {wt: [build_poly_from_exp(exp, base_ring) for exp in exps_2[wt]] for wt in range(1, N+1)}  #2-forms (dxdy)

    """Step 2: initialize linear bases for the F-graded Nygaard-filtered prismatic complex (F^w\\mathcal{N}\\Prism_R^* with 1 <= w <= N; *=0,1,2)"""
    Nygaard_exps = {0:{}, 1:{}, 2:{}} 
    for wt in range(1, N+1):    
        Nygaard_exps[0][wt] = [(e0*p**(max(i-e3-e4, 0)),e1,e2,e3,e4,e5,e6) for (e0,e1,e2,e3,e4,e5,e6) in exps_0[wt]]
        Nygaard_exps[1][wt] = [(e0*p**(max(i-e3-e4-1, 0)),e1,e2,e3,e4,e5,e6) for (e0,e1,e2,e3,e4,e5,e6) in exps_1[wt]]
        Nygaard_exps[2][wt] = [(e0*p**(max(i-e3-e4-2,0)),e1,e2,e3,e4,e5,e6) for (e0,e1,e2,e3,e4,e5,e6) in exps_2[wt]]

    Nygaard_bases = {0:{}, 1:{}, 2:{}}
    for wt in range(1, N+1): 
        Nygaard_bases[0][wt] = [build_poly_from_exp(term, base_ring) for term in Nygaard_exps[0][wt]]
        Nygaard_bases[1][wt] = [build_poly_from_exp(term, base_ring) for term in Nygaard_exps[1][wt]]
        Nygaard_bases[2][wt] = [build_poly_from_exp(term, base_ring) for term in Nygaard_exps[2][wt]]

    """Step 3: collate linear bases for the prsimatic and Nygaard-filtered prismatic complexes"""
    # Use Rem 5.1.1 to split computation of total complex into rays indexed by integers I_p coprime to p and less than N+1
    I_p = [j for j in range(1, N+1) if gcd(p, j)==1]
    rays = {}
    for j in I_p:
        k=0
        rays[j] = []
        while j*p**k < N+1:
            rays[j].append(j*p**k)
            k+=1

    # Z_p-modules in the total chain complex: N^0 -> N^1+Prism^0 -> N^2+Prism^1 -> Prism^2; which is quasiisomrphic to Z_p(i)(R)
    total_Nygaard_bases = {0:{}, 1:{}, 2:{}}
    total_prism_bases = {0:{}, 1:{}, 2:{}}

    for j in I_p:
        for deg in range(3):
            total_Nygaard_bases[deg][j] = []
            total_prism_bases[deg][j] = []

        for j_mult in rays[j]:
            #N^0, N^1, N^2
            total_Nygaard_bases[0][j] += Nygaard_bases[0][j_mult]
            total_Nygaard_bases[1][j] += Nygaard_bases[1][j_mult] 
            total_Nygaard_bases[2][j] += Nygaard_bases[2][j_mult]

            #P^1,P^2,P^3
            total_prism_bases[0][j] += prismatic_bases[0][j_mult]
            total_prism_bases[1][j] += prismatic_bases[1][j_mult]
            total_prism_bases[2][j] += prismatic_bases[2][j_mult]


    """Subalgorithm: reduce"""
    # See reduction.py
    reduce = lambda f: reduce_poly(f, a, b, c) 


    """Step 4: implement internal differentials, canonical map, and Frobenius map"""
    # For differentials, see differential.py
    # The canonical map is the inclusion of the Nygaard-filtered prismatic complex into the prismatic complex
    # For the Frobenius map, see frobenius.py


    """Step 5: implement above maps for the truncated cochain complexes (1<= F-degree <= N) in terms of chosen bases"""
    diff = lambda f: reduce(diff_function(f,a,b,c, base_ring))
    can = lambda f: reduce(f)
    if fast_frob:
        frob = lambda f: reduce(faster_Frob(f, p, i, a, b, c, N, Ring=base_ring))
    else:
        frob = lambda f: reduce(cap_deg(Frob_i(f, p, i, a, b, c, Ring=base_ring), a,b,c,N))

    """Step 6: compute the matrices of the above maps in terms of chosen bases, and thus the syntomic complex Z_p(i)(R)"""
    # The total complex computing the homotopy fiber of can-phi_i^[1,N]:N^* -> P^* is a model for the syntomic complex Z_p(i)(R).
    syntomic_complex_bases = {0:{}, 1:{}, 2:{}, 3:{}}
    syntomic_complex_matrices = {0:{}, 1:{}, 2:{}}  

    # Z_p modules of the homotopy fiber: N^0, N^1+P^0, N^2+P^1, P^2
    for j in I_p:
        syntomic_complex_bases[0][j] = total_Nygaard_bases[0][j]
        syntomic_complex_bases[1][j] = total_Nygaard_bases[1][j] + total_prism_bases[0][j]
        syntomic_complex_bases[2][j] = total_Nygaard_bases[2][j] + total_prism_bases[1][j]
        syntomic_complex_bases[3][j] =                             total_prism_bases[2][j]
    
    # transition functions; if: elt in N^*; else: elt in P^*
        syn_0 = lambda t: diff(t) + can(t) - frob(t)  
        syn_1 = lambda t: diff(t) + can(t) - frob(t) if t in total_Nygaard_bases[1][j] else -diff(t)
        syn_2 = lambda t:           can(t) - frob(t) if t in total_Nygaard_bases[2][j] else -diff(t)  

    # transition matrices
        syntomic_complex_matrices[0][j] = matrix_of(syntomic_complex_bases[0][j], syntomic_complex_bases[1][j], syn_0)
        syntomic_complex_matrices[1][j] = matrix_of(syntomic_complex_bases[1][j], syntomic_complex_bases[2][j], syn_1)
        syntomic_complex_matrices[2][j] = matrix_of(syntomic_complex_bases[2][j], syntomic_complex_bases[3][j], syn_2)


    if verbose:
        print("Syntomic complex computed successfully.")
        print(f"Parameters: p={p}, i={i}, a={a}, b={b}, c={c}.")
        print(f"F-degrees: |x|=|dx|={c}, |y|=|dy|={b}, |u|={a*c}, |v|={b*c}.")
        print(f"Suffices to truncate computation at N={N}.")
        print(f"The direct summands of the total complex are indexed by the integers I_p = {I_p}.")
        for j in I_p:
            print(f"\nFor j={j}, the direct summand contributed to the total complex is the cochain complex:")
            print(syntomic_complex_bases[0][j])
            print(syntomic_complex_matrices[0][j])
            print(syntomic_complex_bases[1][j])
            print(syntomic_complex_matrices[1][j])
            print(syntomic_complex_bases[2][j])
            print(syntomic_complex_matrices[2][j])
            print(syntomic_complex_bases[3][j])

    return I_p, syntomic_complex_bases, syntomic_complex_matrices


def syntomic_cohomology(I_p, syntomic_bases, syntomic_matrices, p, verbose=1):
    """
    Computes the cohomology of the cochain complex model of Z_p(i)(R) produced by the main algorithm.
    Indeed, the cohomology of a complex can be computed from the elementary divisors of the transition matrices.
    See https://www.matem.unam.mx/~omar/mathX27/smith-form.html for a handy reference.
    
    Input:
    - I_p: a list of integers coprime to p and less than N+1, which index the direct summands of the total complex
    - syntomic_bases: a dictionary containing the bases for the modules of the syntomic complex, indexed by F-degree and j
    - syntomic_matrices: a dictionary containing the transition matrices of the maps in the syntomic complex, indexed by F-degree and j
    - p: prime number
    - verbose: 0 -> no info printed; 1 -> readable result of the computation printed; 2 -> additional details on the contribution of each j in I_p

    Returns:
    - A dictionary containing the computed syntomic complex and related information.
    """

    elem_divisors = {0:{}, 1:{}, 2:{}}
    nt_elem_divisors = {0:{}, 1:{}, 2:{}}
    dims_bases = {0:{}, 1:{}, 2:{}, 3:{}}
    ranks = {0:{}, 1:{}, 2:{}}
    homology_dims = {0:{}, 1:{}, 2:{}, 3:{}}

    for j in I_p:
        #elementary divisors of each transition matrix (a_i in the reference link; SageMath method already discards zero ones)
        elem_divisors[0][j] = syntomic_matrices[0][j].elementary_divisors()
        elem_divisors[1][j] = syntomic_matrices[1][j].elementary_divisors()
        elem_divisors[2][j] = syntomic_matrices[2][j].elementary_divisors()

        #elementary divisors that are >1; i.e. those that contribute to cohomology
        nt_elem_divisors[0][j] = sorted([ed for ed in elem_divisors[0][j] if ed>1])
        nt_elem_divisors[1][j] = sorted([ed for ed in elem_divisors[1][j] if ed>1])
        nt_elem_divisors[2][j] = sorted([ed for ed in elem_divisors[2][j] if ed>1])

        #dimensions of bases of the Z_p-modules in cochain complex modeling Z_p(i)(R) (namely, N^\bullet+P^\bullet[1]; \bullet=0,1,2)
        dims_bases[0][j] = len(syntomic_bases[0][j])
        dims_bases[1][j] = len(syntomic_bases[1][j])
        dims_bases[2][j] = len(syntomic_bases[2][j])
        dims_bases[3][j] = len(syntomic_bases[3][j])

        #ranks of the transition matrices (r in the link above)
        ranks[0][j] = len(elem_divisors[0][j])
        ranks[1][j] = len(elem_divisors[1][j])
        ranks[2][j] = len(elem_divisors[2][j])

        #ranks of H^0, H^1, H^2 of the synotmic complex (m-r-s in link above)
        homology_dims[0][j] = dims_bases[0][j] - ranks[0][j]
        homology_dims[1][j] = dims_bases[1][j] - ranks[1][j] - ranks[0][j]
        homology_dims[2][j] = dims_bases[2][j] - ranks[2][j] - ranks[1][j]
        homology_dims[3][j] = dims_bases[3][j]               - ranks[2][j]

        
    #combine the info from different j in I_p
    nt_elem_divs_0 = []
    nt_elem_divs_1 = []
    nt_elem_divs_2 = []
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

    for j in I_p:
        nt_elem_divs_0 += nt_elem_divisors[0][j]
        nt_elem_divs_1 += nt_elem_divisors[1][j]
        nt_elem_divs_2 += nt_elem_divisors[2][j]
        dims_bases_0 += dims_bases[0][j]
        dims_bases_1 += dims_bases[1][j]
        dims_bases_2 += dims_bases[2][j]
        dims_bases_3 += dims_bases[3][j]
        rank_0 += ranks[0][j]
        rank_1 += ranks[1][j]
        rank_2 += ranks[2][j]
        #homology_dims are the ranks of the homology (number of free generators); they should be zero
        homology_dim_0 += homology_dims[0][j]
        homology_dim_1 += homology_dims[1][j]
        homology_dim_2 += homology_dims[2][j]
        homology_dim_3 += homology_dims[3][j]

    cohomology = {1:sorted(nt_elem_divs_0), 2:sorted(nt_elem_divs_1), 3:sorted(nt_elem_divs_2)}
    
    # raise warning if we ever get non-torsion stuff (for theoretical reasons, we shouldn't!) 
    # or if H3 doesn't vanish (we conjecture it always does)
    if homology_dim_0+homology_dim_1+homology_dim_2+homology_dim_3 >0 or len(nt_elem_divs_2)>0:
        synt_data = {
        'p': p,
        'I_p': I_p,
        'syntomic_bases': syntomic_bases,
        'syntomic_matrices': syntomic_matrices,
        'dims_bases': dims_bases,
        'ranks': ranks,
        'homology_dims': homology_dims,
        'elem_divisors': elem_divisors,
        'nt_elem_divisors': nt_elem_divisors,
        'cohomology': cohomology
        }

        with open('log.py', 'w') as f:
            f.write(str(synt_data))
            f.close()
        raise ValueError("look at log!")
    
    # readable output of the computation
    if verbose >= 1:
        print("Syntomic cohomology computed successfully.")
        print("The syntomic cohomology groups are as follows:")
        print(f"H^0(Z_p(i)(R)) = 0")
        for d in range(1,4):
            print(f"H^{d}(Z_p(i)(R)) = {format_group_str(p_power_count(elem_divs_list=cohomology[d], p=p), p=p)}")
        
    # optional details about the computation
    if verbose >= 2:
        for j in I_p:
            
            print(f"\nFor j={j}, the contribution to syntomic cohomology (i.e. the non-trivial elementary divisors) are cyclic groups of order:")
            print("H^0:", nt_elem_divisors[0][j],  "\nH^1:", nt_elem_divisors[1][j], "\nH^2:", nt_elem_divisors[2][j])

    return cohomology


def syntomic_data(p, i, a, b, c, fast_frob=True):
    """
    Computes the syntomic complex and its cohomology for given parameters.
    
    Parameters:
    - p: prime number
    - i: syntomic weight (non-negative integer)
    - a, b, c: positive integers
    - verbose: if True, prints detailed information about the computation
    - fast_frob: if True, uses a faster version of the Frobenius map

    Returns:
    - A dictionary containing information about the syntomic cohomology computation.
    """
    
    if b > a:
        return syntomic_data(p, i, a, c, a)

    # N is the maximum F-degree of the polynomials in the bases to be computed
    # Crucially, F^{>N}Z_p(i)(R) = 0, as explained in the thesis
    N = a * c * i + (b - 1) * (c - 1) - 1

    start_time = time.time()
    
    I_p, syntomic_bases, syntomic_matrices = syntomic_complex(p=p, i=i, a=a, b=b, c=c, fast_frob=fast_frob)
    
    syntomic_cx_computation_time = time.time() - start_time
    
    cohomology = syntomic_cohomology(I_p, syntomic_bases, syntomic_matrices, p, verbose = 0)

    cohomology_computation_time = time.time() - start_time - syntomic_cx_computation_time

    # As the groups get large, it is more convenient to encode them as follows
    torsion_exps_dict = {0:{}, 1:{}, 2:{}}
    for d in range(1,4):
        torsion_exps_dict[d] = p_power_count(cohomology[d], p)

    total_time = time.time() - start_time

    synt_data = {
                'computed': True, 
                'total_time':total_time,
                'synt_cx_time': syntomic_cx_computation_time, 
                'cohom_time':cohomology_computation_time,
                'p':p, 
                'i':i,
                'a':a, 
                'b':b, 
                'c':c,  
                'F_deg_cutoff':N, 
                'I_p':I_p,
                'length':a*c,                   # length of R as a F_p-module 
                'H^0':{}, 
                'H^1':torsion_exps_dict[1], 
                'H^2':torsion_exps_dict[2], 
                'H^3':torsion_exps_dict[3]
                }

    return synt_data

def p_adic_K_thry_of_cusp(p, i, a, b, c, fast_frob=True):
    """
    Computes the (2*i-2)-th and (2*i-1)-th p-adic K-group of the cusp R=F_p[x,y]/(x^a, x^b-y^c).
    
    Parameters:
    - p: prime number
    - i: nonnegative integer (syntomic weight)
    - a, b, c: positive integers
    - verbose: if True, prints detailed information about the computation
    - fast_frob: if True, uses a faster version of the Frobenius map

    Output:
    - Prints readable strings detailing the structure of the K-groups.
    """

    if not is_prime(p):
        raise ValueError(f"p={p} is not a prime.")
    

    # i-th syntomic cohomology group
    synt_data = syntomic_data(p=p, i=i, a=a, b=b, c=c, fast_frob=fast_frob)
    print(f"Syntomic cohomology of R=F_p[x,y]/(x^{a}, x^{b}-y^{c}) computed successfully.")
    print("The syntomic cohomology groups are as follows:")
    for d in range(4):
        print(f"H^{d}(Z_{p}({i})(R)) = {format_group_str(synt_data[f'H^{d}'], p=p)}")
    

    # Passage to K-theory
    synt_data_next = syntomic_data(p=p, i=i+1, a=a, b=b, c=c, fast_frob=fast_frob)
    K_2i_2 = format_group_str(synt_data['H^2'], p=p)
    print(f"\nThe p-adic K-group K_{2*i-2}(R; Z_{p}) is: {K_2i_2}.")

    if len(synt_data_next['H^3']) == 0:
        # Presumably, this 'if' clause always triggers; see discussion after Corollary 1.0.2 or Conjecture 4.3.3.
        K_2i_1 = format_group_str(synt_data['H^1'], p=p)
        print(f"The p-adic K-group K_{2*i-1}(R; Z_{p}) is: {K_2i_1}.")
    
    else:
        # Presumably, this clause never triggers.
        if p == 2:
            print(f"The p-adic K-group K_{2*i-1}(R; Z_{p}) fits in the short exact sequence: ")
            print(f"0 -> {format_group_str(synt_data_next['H^3'], p=p)}  -> K_{2*i-1}(R; Z_{p}) -> {format_group_str(synt_data['H^1'], p=p)} -> 0.")
            
        elif p > 2:
            print(f"The p-adic K-group K_{2*i-1}(R; Z_{p}) is: {format_group_str(combine_power_counts(synt_data['H^1'], synt_data_next['H^3']), p=p)}.")


def main():
    parser = argparse.ArgumentParser(description="Compute p-adic K-theory in degrees 2i-2 and 2i-1 for a ring R=F_p[x,y]/(x^a, x^b-y^c).")
    parser.add_argument('p', type=int, help="A prime number p.")
    parser.add_argument('i', type=int, help="Syntomic weight i.")
    parser.add_argument('a', type=int, help="Exponent x^a.")
    parser.add_argument('b', type=int, help="Exponent x^b.")
    parser.add_argument('c', type=int, help="Exponent y^c.")

    args = parser.parse_args()

    p = args.p
    a, b, c, i = [args.a, args.b, args.c, args.i]

    p_adic_K_thry_of_cusp(p=p, i=i, a=a, b=b, c=c)

if __name__ == "__main__":
    main()

