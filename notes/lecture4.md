---
title: "Lecture 4"
author: James Camacho
geometry: margin=2cm
output: pdf_document
colorlinks: true
---

Covering:

- Finite element methods,
- Various bases,
- Fourier/DCT transforms,
- Fourier analysis (stability/convergence)
- Fourier analysis (faster solvers)

-----

# Finite Element Methods

In finite element methods, you approximate a function $f$ as a sum of trial functions, $$f(x)\approx \sum_{j=1}^n c_j\Phi_j,$$then test the approximation using test functions $\{\Psi_i\}_{i=1}^m$. The test is whether $$\int_{0}^1\sum_{j=1}^nc_j\Phi_j(x)\Psi_\text{test}(x)\text{d}x = \int_{0}^1 f(x)\Psi_\text{test}(x)\text{d}x$$for every test function. To save room, mathematicians use the notation $\langle u, \Psi\rangle$ to refer to the  [inner product](https://mathworld.wolfram.com/InnerProduct.html) $$\int_{0}^1 f(x)\Psi(x)\text{d}x.$$Sometimes it's useful to insert a nonnegative weight $w$, so $$\langle f,\Psi\rangle_w = \int_0^1 f(x)\Psi(x)w(x)\text{d}x.$$Writing the test equations in matrix form gives $$\begin{bmatrix}\langle \Psi_1, \Phi_1\rangle&\langle \Psi_1, \Phi_2\rangle&\cdots&\langle \Psi_1, \Phi_n\rangle\\\langle \Psi_2, \Phi_1\rangle&\langle \Psi_2, \Phi_2\rangle&\cdots&\langle \Psi_2, \Phi_n\rangle\\\vdots&\vdots&\ddots&\vdots\\\langle \Psi_m, \Phi_1\rangle&\langle \Psi_m, \Phi_2\rangle&\cdots&\langle \Psi_m, \Phi_n\rangle\end{bmatrix}\begin{bmatrix}c_1\\c_2\\\vdots\\c_n\end{bmatrix} = \begin{bmatrix}\langle \Psi_1, f\rangle\\\langle \Psi_2, f\rangle\\\vdots\\\langle \Psi_m, f\rangle\end{bmatrix},$$or in [bra-ket notation](https://en.wikipedia.org/wiki/Bra%E2%80%93ket_notation), $$\langle \Psi|\langle\Phi|\textbf{c}\rangle|w\rangle = \langle\Psi |f| w\rangle.$$To solve the Poisson equation ($\nabla^2 u = f$), we can approximate $u$ as a sum of basis functions, and solve $$\langle \Psi|\langle\nabla^2\Phi|\textbf{c}\rangle|w\rangle = \langle\Psi |f| w\rangle.$$
Then $u\approx \langle \Phi | \textbf{c}\rangle$.

Usually the test and trial bases are the same. Also, we can solve the system faster if the matrix has a small bandwidth; in particular orthogonal functions ($\langle \Phi_i, \Phi_j\rangle = 0$ for $i\ne j$) will give a diagonal system. In addition, it's good if they are easy to integrate and differentiate. Some common bases include hat functions, Legendre polynomials, and the Chebyshev polynomials.

For large systems, you divide the structure into a mesh and choose a basis for each cell (usually the same basis but shifted over). Outside its cell each function is defined to be zero so you only need to take inner products between functions within a cell. To keep it smooth across cells, you create boundary conditions, evaluating the function and its derivatives at the boundary. The Bernstein bases are particularly good for this, as although they are not orthogonal most of them vanish on the boundary ([source](https://boundaryvalueproblems.springeropen.com/articles/10.1155/2011/829543)).

-----

# Common Bases
## Hat Functions

![Hat Basis](images/hat_basis.png){ width=500px }

The hat functions have a peak at the center of their cell. An equivalent basis would put a line with positive slope and a line with negative slope in each cell, so the hat functions are also called the "piecewise linear basis". We can extend them to higher dimensions using a tensor product, $$\Phi^{2D} = \Phi^{1D}\otimes\Phi^{1D}\Longleftrightarrow \Phi_{ij}^{2D} = \Phi_i^{1D}\Phi_j^{1D}.$$Here I've shown the top half of some two dimensional hat functions (the pyramids should continue down after intersection):

![Two Dimensional Hat Basis](hat2d_basis.png){ width=500px }

There is no weight for these bases, i.e. $w\equiv 1$. For solving the Poisson equation, the second derivative is zero, so $$\langle \Psi|\langle\nabla^2\Phi|\textbf{c}\rangle|1\rangle = \langle\Psi |f| 1\rangle$$looks like it'll become $$0 = \langle\Psi|f|1\rangle,$$which doesn't really work. Luckily, the product rule comes to the rescue. We have $$\nabla\cdot (\Psi\nabla\Phi) = \nabla\Psi\cdot \nabla\Phi + \Psi\nabla^2\Phi\implies \int_{0}^{1}\Psi\nabla^2\Phi\text{d}x = \Psi\nabla\Phi\biggr\rvert_0^1 - \int_0^1\nabla\Psi\cdot\nabla\Phi\text{d}x.$$As $$\Psi(0) = \Psi(1) = 0,$$this reduces to $$\langle \Psi, \nabla^2\Phi\rangle = -\langle\nabla\Psi, \nabla\Phi\rangle.$$This trick is often called integration by parts. We need to solve $$-\langle\nabla \Psi|\langle\nabla\Phi|\textbf{c}\rangle|w\rangle = \langle\Psi |f| 1\rangle.$$In one dimension we find $$-\langle \Psi_i', \Phi_j'\rangle = \begin{cases}-1&i=j\\0.5&|i-j|=1\\0&\text{otherwise}.\end{cases}$$The finite element matrix ends up being the exact same as the finite difference one, times a constant. This happens in general for all dimensions, e.g. the two dimensional hat functions will give rise to the nine-point stencil finite difference matrix. The only difference from the finite difference equation is the right hand side, where we convolve $f$ with the hat function rather than taking a single point.

## Legendre Polynomials

![](legendre.png){ width=500px }

Again, no weight for this basis, although the bounds on the integral go from minus one to one: $$\langle P_i, P_j\rangle = \int_{-1}^{1} P_i(x)P_j(x)\text{d}x.$$The Legendre polynomials are orthogonal, and can be computed via Bonnet's recursion formula: $$P_0=1;\quad P_1=x;\quad P_{n+1} = \frac{(2n+1)xP_n - nP_{n-1}}{n+1}.$$Due to orthogonality, the approximation formula $$\langle \Psi|\langle\Phi|\textbf{c}\rangle|w\rangle = \langle\Psi |f| w\rangle,$$reduces to $$c_i = \frac{\langle P_i, f\rangle}{\langle P_i, P_i\rangle}.$$Solving the Poisson equation would give a much denser matrix.

## Chebyshev Polynomials

![](chebyshev.png){ width=500px }

The Chebyshev polynomials are defined by $$T_n(\cos(x)) = \cos(nx),$$and satisfy the recurrence $$T_0=1;\quad T_1=x;\quad T_{n+1} = 2xT_n - T_{n-1}.$$Under the inner product $$\langle T_i, T_j\rangle = \int_{-1}^1 \frac{T_i(x)T_j(x)}{\sqrt{1-x^2}}\text{d}x$$(i.e. $w(x) = 1/\sqrt{1-x^2}$) we get $$\langle T_i, T_j\rangle = \begin{cases}\pi&i=j=0\\\pi/2&i=j\ne0\\0&\text{otherwise},\end{cases}$$so the approximation formula would reduce to $$c_0 = \frac1\pi\langle f, 1\rangle;\quad c_{i} = \frac{2}{\pi}\langle f, T_i\rangle.$$Integrating and multiplying together Chebyshev polynomials is also pretty quick: $$\begin{aligned}\int T_n \text{d}x = \frac12\left(\frac{T_{n+1}}{n+1}-\frac{T_{n-1}}{n-1}\right),\\T_aT_b = \frac12(T_{a+b}+T_{a-b}).\end{aligned}$$Unfortunately taking derivatives is not so nice, and so the Poisson formula ends up with a pretty dense matrix.
