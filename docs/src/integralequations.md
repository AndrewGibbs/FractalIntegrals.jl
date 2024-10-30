```@setup tutorial
using FractalIntegrals
using Plots; gr()
Plots.reset_defaults()
```

# Integral Equations

In [the previous section](integrals.md) we saw how to approximate integrals on fractals. Now, we turn our attention to solving integral equations on fractals. Many constant coefficient PDEs can be reformulated as integral equations, and in some circumstances, the integral equation is easier to solve than the original PDE.

## A brief introduction to integral equations on fractals

When solving problems on fractals, a common approach is to approximate the fractal $\Gamma$ by a shape that $\tilde{\Gamma}\approx\Gamma$ that is sufficiently detailed, but non-fractal, and instead solve the same problem on $\tilde{\Gamma}$. This is appealing because it enables the use of standard methods and/or software. However, this geometric approximation introduces an error that is difficult to quantify.

Recently, in the context of integral equations [ProcRSCaetal:23, HausdorffBEM](@cite), techniques have been developed for solving problems on $\Gamma$, without the need for a geometrical approximation $\tilde{\Gamma}$.

Earlier versions of `FractalIntegrals.jl` were used in two 

For our purposes, we consider integral equations of the form:

```math
(\lambda~\mathcal{I} + \mathcal{K})~v= f,\quad\text{on }\Gamma,
```
where $v$ is the unknown solution, $\lambda\in\mathbb{C}$, $\mathcal{I}$ is the identity operator, $f$ is known data and $\mathcal{K}$ is an integral operator of the form
```math
\mathcal{K}~\phi~(x)= \int_\Gamma K(x,y)~\phi(y)~ \mathrm{d}\mu(y),\quad x\in\Gamma,
```
where $K$ is a known integral kernel.

For our purposes, we make interpret $\Gamma$ as an attractor, or a union of attractors. Typically, when one has solved for $v$, a *potential operator* is used to obtain a solution to the equivalent PDE
```math
\mathcal{V} v(x) = \int_\Gamma ~V(x,y)~v(y)~ \mathrm{d}\mu(y), \quad x \in \mathbb{R}^N.
```

## Classification

There are two key classes of integral equations, and this distinction is particularly significant in the fractal setting:

* **1st kind**: When $\lambda=0$, equivalently, the identity term is removed.
* **2nd kind**: When $\lambda\neq0$.

As mentioned in [](geometry.md#measures), for many attractors $\Gamma$, we do not know $\mathcal{H}^d(\Gamma)$. In this case, we cannot work with $\mu=\mathcal{H}^d$ directly, instead we must work with $\tilde\mu(x) = \mathcal{H}^d(x)/\mathcal{H}^d(\Gamma)$. Hence, if the integral equation is defined with respect to Hausdorff measure, we can (usually) only solve a *related* integral equation with respect to $\tilde\mu$.

In the case of first kind integral equations, the solution to the *related* problem $\tilde v=v\cdot \mathcal{H}^d(\Gamma)$. If the *related* potential operator $\tilde{\mathcal{V}}$ is defined analogously with respect to $\tilde\mu$, then $\mathcal{V}v = \tilde{\mathcal{V}} \tilde v$; in other words, the solution of the underlying PDE is unaffected by our scaling of the measure, thus **we do not need to know** $\mathcal{H}^d(\Gamma)$.

In the case of second kind integral equations, we are not so fortunate: $\tilde v \neq v\cdot \mathcal{H}^d(\Gamma)$. Therefore, when posing a second-kind integral equation, for meaningful results, it is essential that the `suppmeasure` parameter, mathematically defined as $\mu(\Gamma)$, is known and included in the definition of the measure. We note that for the case when $\Gamma$ is open in $\mathbb{R}^N$ with a fractal boundary, for example the Koch Snowflake, the $N$-dimensional Lebesgue measure is usually known.

## Fractal Operators

The type hierarchy for `FractalOpeartor` is as follows:
```
FractalOpeartor
├── IntegralOperator
    ├── AbstractSingularIntegralOperator
        ├── SingularIntegralOperator
        ├── SeparableIntegralOperator
    ├── SmoothIntegralOperator
├── IdentityOperator
├── ScaledOperator
├── SumOperator
├── BlockOperatorx

```