Lecture 3 hw problems:
1. Given that $a_1 + a_2 + \cdots + a_n = 0$ and $a_1^2 + a_2^2 + \cdots + a_n^2 = 1$, find the largest possible value of $a_1a_2 + a_2a_3 + \cdots + a_na_1$. Or, writing with inner product form, given $\langle a, 1\rangle = 0$ and $\langle a, a\rangle = 1$, what is the largest possible value of $\langle a, Ra\rangle$ where $R$ is the "rotate right" matrix?

For lecture 3 examples:
- Find roots of Legendre polynomials using QR algorithm.

## Finite Element Methods

Suppose we have some function $f$, and a bunch of trial functions $\Phi_1, \Phi_2, \dots, \Phi_n$. We approximate $f$ as a sum of these trial functions: $$f(x)\approx \sum_{i=1}^n c_i\Phi_i$$for some coefficients $c_i$. Then, we test how good the approximation is using other test functions, $\Psi_1, \Psi_2, \dots, \Psi_m$. Usually the test functions are the same as the trial functions (so $m=n$ and $\Psi_i=\Phi_i$), but they do not have to be. The test is whether the two integrals match up: $$\int_{-1}^1\sum_{i=1}^nc_i\Phi_i(x)\Psi_\text{test}(x)\text{d}x = \int_{-1}^1 f(x)\Psi_\text{test}(x)\text{d}x.$$In order to save room, most mathematicians use the notation $\langle f, \Psi\rangle$ to refer to $$\int_{-1}^1 f(x)\Psi(x)\text{d}x.$$In general the notation $\langle \cdot, \cdot\rangle$ refers to an [inner product](https://mathworld.wolfram.com/InnerProduct.html), and an integral happens to be one inner product. Anyway, we can pull the summation outside the integral and write it in terms of inner products to get the test $$\sum_{i=1}^nc_i\langle \Phi_i, \Psi_\text{test}\rangle = \langle f, \Psi_\text{test}\rangle.$$Remember, we have $m$ different test functions, so if $m = n$ we have $n$ equations and $n$ unknown $c_i$, meaning we can find a unique solution. It becomes a lot easier if we make two additional assumptions:
- The test and trial functions are the same.
- $\langle \Phi_i, \Phi_j\rangle = \delta_{ij}$ (one when $i=j$ and zero otherwise, known as the [Kronecker delta](https://en.wikipedia.org/wiki/Kronecker_delta)).
Then our equations simplify to $$\sum_{i=1}^n c_i\langle\Phi_i, \Psi_\text{test}\rangle = c_i,$$and we're left with $c_i = \langle f, \Phi_i\rangle$.

To be continued...