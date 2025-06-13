# Syntomic Cohomology of Truncated Cuspidal Curves

## Overview

Implements algorithm to compute the $p$-adic syntomic cohomology $\mathbb{Z}_p(i)(R)$ of rings of the form $R=F_p[x,y]/(x^a, x^b-y^c)$ at arbitrary levels $i\geq 0$. As a corollary, one also obtains individual $p$-adic algebraic K-groups of $R$.

The code was implemented with ``SageMath 10.4``.
## Example usage

To compute the syntomic cohomology and $K$-theory corresponding to the input $(p,i,a,b,c) = (3,2,3,2,3)$ one can run

```python syntomic.py 3 2 3 2 3```

from the CLI. It outputs
```Syntomic cohomology of R=F_p[x,y]/(x^3, x^2-y^3) computed successfully.
The syntomic cohomology groups are as follows:
H^0(Z_3(2)(R)) = 0
H^1(Z_3(2)(R)) = (Z/3)^13 + (Z/3^2)^3
H^2(Z_3(2)(R)) = (Z/3)^3
H^3(Z_3(2)(R)) = 0

The p-adic K-group K_2(R; Z_3) is: (Z/3)^3.
The p-adic K-group K_3(R; Z_3) is: (Z/3)^13 + (Z/3^2)^3.
```
as the result of the computation.

The notebook ``example.ipynb`` illustrates this computation in greater detail as well as the usage of the functions implemented in ``syntomic.py``.

## References

Cortez Lemos, Carlos. _Syntomic Cohomology of Truncated Cuspidal Curves_. Northwestern University, 2024. ProQuest Dissertations & Theses, [https://www.proquest.com/docview/3097700872](https://www.proquest.com/docview/3097700872).

The Sage Developers. (2024). *SageMath, the Sage Mathematics Software System (Version 10.4)*. [https://www.sagemath.org](https://www.sagemath.org).
