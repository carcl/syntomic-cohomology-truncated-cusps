from sage.all import PolynomialRing, QQ, ZZ, Matrix, valuation
from collections import Counter

def matrix_of(basis_from, basis_to, func, Ring = PolynomialRing(QQ, 'x,y,u,v,dx,dy')):
    """
    Compute the matrix of a linear transformation defined by `func`
    from the basis `basis_from` to `basis_to` over a polynomial ring.
    """
    rows = len(basis_to)
    cols = len(basis_from)
    M = Matrix(ZZ, rows, cols)

    basis_to_monomials = [term.monomials()[0] for term in basis_to]

    for col_index, basis_elem in enumerate(basis_from):
        transformed = Ring(func(basis_elem))
        for monomial in transformed.monomials():
            full_term = monomial * transformed.monomial_coefficient(monomial)
            if monomial in basis_to_monomials:
                row_index = basis_to_monomials.index(monomial)
                coeff_in_func = full_term.coefficients()[0]
                coeff_in_basis = basis_to[row_index].coefficients()[0]
                M[row_index, col_index] += coeff_in_func / coeff_in_basis

    return M

def p_power_count(elem_divs_list, p):
    """
    Encode the direct sum decomposition of a finite abelian p-group.
    Ex: (Z/p)^5 + (Z/p^3) is encoded by {1: 5, 3:1} 
    """
    count = Counter()
    for elem_div in elem_divs_list:
        e = elem_div.valuation(p)
        if p**e != elem_div:
            raise ValueError(f"{elem_div} is not a pure power of {p}")
        count[e] += 1
    
    return dict(count)


def combine_power_counts(p_power_dict_1, p_power_dict_2):
    """
    Outputs the dictionary encoding the direct sum of two abelian groups encoded as by p_power_count.
    """
    return dict(Counter(p_power_dict_1) + Counter(p_power_dict_2))

def format_group_str(p_power_dict, p):
    """
    Convert a dictionary encoding a group to a string representation of the group.
    Eg, input: {1: 1, 3: 2, 4: 1}; output: 'Z/5 + (Z/5^3)^2 + Z/5^4'
    """
    if len(p_power_dict) == 0:
        return 0
    
    terms = []
    for e in sorted(p_power_dict):
        count = p_power_dict[e]
        group = f"Z/{p}^{e}" if e > 1 else f"Z/{p}"
        if count == 1:
            terms.append(group)
        else:
            terms.append(f"({group})^{count}")
    return " + ".join(terms)
