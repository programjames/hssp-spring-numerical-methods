# Numerical Methods
## An HSSP Spring Course

**Difficulty:** High.

**Required Knowledge:**
- Middle school algebra.
- *If you wish to do homework* you should know Python. You would probably be fine with MATLAB, C++, or C#, but all code here will be written in Python with the packages `numpy, scipy, matplotlib,` and `jupyter`.

-----

# Schedule (tentative)

1. The derivative, fixed point iteration, Newton's, Euler & Runge-Kutta methods.
2. Gauss-Legendre quadrature, solving PDE's (e.g. heat equation), finite difference methods & correctors.
3. Matrices: finding eigenvalues/vectors, Broyden's method, condition number & stiffness. Common bases and their condition numbers (i.e. why polynomials are bad).
4. Finite element methods w/ Chebyshev polynomials or sines/cosines. Fourier/DCT transforms, Fourier analysis (faster solvers) & Fourier analysis (determining convergence rate).
5. Optimization: golden section, simplex, gradient descent, conjugate gradient method, Adam.
6. Different ideas, maybe a combination of them: Image recognition w/ Bayes' theorem, the Lagrangian w/ the simplex method, image compression w/ principal component analysis.

-----

# Lecture Notes (homeworks at the end)

1. [Week 1](notes/lecture1.pdf)
2. [Week 2](notes/lecture2.pdf)
3. [Week 3](notes/lecture3.pdf)
4. [Week 4](notes/lecture4.pdf)

-----

# Installation
Open up the terminal (cmd line on Windows) and run
```shell
git clone https://github.com/programjames/hssp-spring-numerical-methods
cd hssp-spring-numerical-methods
pip install -r requirements.txt
```

To view a `.ipynb` file run `jupyter notebook` and open the corresponding file.
