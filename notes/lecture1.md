## By James Camacho

Covering:
- The derivative + trig basics,
- Fixed point iteration,
- Newton's method,
- Euler & Runge-Kutta methods.
-----
### The Derivative
The equation for a line is $y = mx + b$, where $m$ is the "slope" of the line—the change in $y$ over the change in $x$, or in math terms, $\frac{\text{d}y}{\text{d}x}$. The "d" comes from "delta", and in fact people also use $\frac{\Delta y}{\Delta x}$ and $\frac{\delta y}{\delta x}$ "(D/d)elta y over (D/d)elta x". There's also a curly "d" $\frac{\partial y}{\partial x}$ that we'll get into later. But, all of these refer to a rate of change, just in slightly different contexts.

Calculus is generalizing the idea of slopes to more functions. For example, what is the slope of the line $y = x^2$ at the point $(x, y) = (1, 1)$?

![[slope_x2.excalidraw]]

We could draw a line tangent to the function at that point and compute its slope. We call this the "derivative". How do you compute this? Well, we just take add a small $\text{d}x$ and find how much that changes $y$, like so:
$$\begin{align*}y &= x^2,\\y + \text{d}y &= (x + \text{d}x)^2\\&\implies \\\text{d}y &= 2x\cdot \text{d}x + (\text{d}x)^2\\&\implies\\\frac{\text{d}y}{\text{d}x} &= 2x + \text{d}x.\end{align*}$$
Of course, when $\text{d}x$ gets really small this is infinitesimally close to $2x$, so as we get better and better approximations to a tangent line, the slope gets closer and closer to the tangent line's. The derivative of $y = x^2$ at $x=1$ is $2x = 2$, so the slope of that line must be two.

Now, what about those other variants, $\Delta, \delta, \partial$? The first two are used pretty interchangeably to mean a *finite difference*, basically the same idea but we don't use arbitrarily small $\text{d}x$, we fix it at some finite value. For example, with $\Delta x = 1$ we would get $$\frac{\Delta y}{\Delta x} = \frac{(x + \Delta x)^2 - x^2}{\Delta x} = 2x + \Delta x = 2x + 1.$$Sometimes people use centered differences, because they can be more accurate, like so:
$$\frac{\Delta y}{\Delta x} = \frac{(x + \Delta x)^2 - (x - \Delta x)^2}{(\Delta x + \Delta x)} = 2x.$$
The last variant, $\partial$ (called the "partial (derivative)"), means we pretend all variables are constants except the one in the denominator. For example, if we made an additional equation $z = x^2 + y$, then $$\frac{\partial z}{\partial x} = 2x,$$but $$\frac{dz}{dx} = 2x + \frac{dy}{dx} = 2x + 2x = 4x$$as we need to take into account that $y$ changes when $x$ changes.

-----
### Derivatives of Common Functions

Some of the most common functions you see are monomials and exponentials. Consider the function $f(x) = x^n$ (a monomial). What is $\frac{\text{d}f}{\text{d}x}$? Well, we have
$$f(x + \text{d}x) = \underbrace{(x+\text{d}x)(x+\text{d}x)\cdots (x+\text{d}x)}_{n\text{ times}}.$$
If we multiply it all out, we can choose either an $x$ or a $\text{d}x$ from each binomial. There is one way to choose all $x$'s, $n$ ways to choose a single $\text{d}x$ and have all the rest be $x$'s, $\frac{n(n-1)}{2}$ ways to choose two $\text{d}x$'s, and so forth:
$$f(x+\text{d}x) = x^n + nx^{n-1}\text{d}x + \frac{n(n-1)}{2}x^{n-2}\text{d}x^2 + \cdots.$$
These coefficients are known as [binomial coefficients](https://en.wikipedia.org/wiki/Binomial_coefficient). We only really care about the first few because $$\frac{\text{d}f}{\text{d}x} = \frac{f(x+\text{d}x)-f(x)}{\text{d}x} = nx^{n-1} + \frac{n(n-1)}{2}x^{n-2}\text{d}x + \cdots,$$
but $\text{d}x$ is infinitesimally small, so everything except $nx^{n-1}$ approaches zero. Therefore, the derivative of $f(x) = x^n$ is $nx^{n-1}$.

Alright, what about exponentials? For example, what is the derivative of $f(x) = 2^x$?
![[exponential.excalidraw]]
The intuition behind compounding interest tells us that the slope should grow at the same pace as the function, e.g. $\frac{\text{d}f}{\text{d}x} = c2^x$ for some constant $c$. Working out the math gives us $$\frac{\text{d}f}{\text{d}x} = \frac{2^{x+\text{d}x} - 2^x}{\text{d}x} = \frac{2^{\text{d}x} - 1}{\text{d}x}2^x.$$As $2^{0} - 1 = 0$, when $\text{d}x = 0$ both the numerator and denominator of $\frac{2^{\text{d}x}-1}{\text{d}x}$ are zero, which doesn't really make sense. However, if we take successively smaller $\text{d}x$, and compute the limit as $\text{d}x$ goes to zero, we find $$\lim_{\text{d}x\to 0}\frac{2^{\text{d}x}-1}{\text{d}x}\approx 0.69.$$So, $$\frac{\text{d}f}{\text{d}x}\approx 0.69\cdot 2^x = 0.69f(x).$$For $f(x) = 3^x$ we get $$\frac{\text{d}f}{\text{d}x}\approx 1.1\cdot 3^x = 1.1f(x).$$It'd be more natural to find derivatives if $\frac{\text{d}f}{\text{d}x} = f,$ so mathematicians defined a base that satisfies this—Euler's number, $e$, which is approximately $2.7183$. Two ways to compute it are with the formulae $$\begin{align*}e &= \lim_{n\to\infty}\left(1+\frac{1}{n}\right)^n\\e&=\frac{1}{0!}+\frac{1}{1!}+\frac{1}{2!}+\frac{1}{3!}+\cdots.\end{align*}$$Call the natural logarithm (ln) the logarithm with this base (so $\ln(7.4)\approx 2$), then for any base $a$, $$a^x = (e^{\ln a})^x = e^{x\ln a}.$$Taking a derivative gives us $$\frac{\text{d} (a^x)}{\text{d}x} = \ln a\cdot\frac{\text{d} (e^{x\ln a})}{\text{d} (x\ln a)} = \ln a\cdot e^{x\ln a} = \ln a\cdot a^x.$$So, throwing back to $f(x) = 3^x$ above, $$\frac{\text{d}(3^x)}{\text{d}x} = \ln 3 \cdot 3^{x},$$and in fact $\ln 3$ is approximately $1.1$.

-----
### Other Calculus
1. People often use prime notation, writing $y'$ instead of $\text{d}y/\text{d}x$. This should only ever be used for the total derivative.
2. The product rule is $(fg)' = f'g + g'f$. Try plugging in the deltas if you wish to prove it.
3. The chain rule is $$\frac{\text{d}(f(g(x))}{\text{d}x} = \frac{\text{d}f(g)}{\text{d}g}\cdot \frac{\text{d}g}{\text{d}x},$$for example with $f(g)=e^{g}$ and $g(x) = x\ln 3$ we get $$\frac{\text{d} e^{x\ln 3}}{\text{d}x} = e^{x\ln 3}\cdot \ln 3.$$I used this a couple paragraphs ago.
4. The *integral* is the area under a function (if the function is negative, that area is considered negative). It can be computed by making a bunch of small boxes of width $\text{d}x$ and height $f(x)$, then adding them up: ![[integral.excalidraw]]
5. The derivative of the integral is the original function. Sometimes the integral is called the derivative for this reason. (Example: The integral of velocity is how far you've moved, i.e. your position, while the derivative of your position is how fast you're moving, i.e. your velocity.)
6. In general integrals are much harder to compute than derivatives, and many functions do not have simple integrals. We'll show how to numerically approximate integrals later on. But first...
-----
### Fixed Point Iteration
![[fixed_point.excalidraw]]

Say you take a function and apply it over and over again. If it's like $f(x) = x^2$ and you start with $x=2$, then the number will get infinitely large. However, if the derivative is between $-1$ and $1$, we get a "contraction mapping" where the range shrinks after every iteration until we're left with just a few fixed points. For example, the only fixed point of $\cos(x)$ is $x\approx 0.7391$.

We can use fixed point iteration to find roots. For example, with the equation $$x^5 - 7x^3 + 3x^2 + x + 2 = 0$$ we can rearrange it to $$x = \sqrt[5]{7x^3 - 3x^2 - 2}$$then apply fixed point iteration on the right-hand side (RHS) to get equality when $x\approx 2.362$.

-----
### Newton's Method
