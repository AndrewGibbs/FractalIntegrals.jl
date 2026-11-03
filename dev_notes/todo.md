# To-do list

*Converted from `todo.odt` on 2026-11-03. This is now the shared copy; the original `.odt` is gitignored.*

Can you turn off macros? If so, it'd be worth having the option to turn off parallelism.

Variable names:

- measure / \mu
- selfmeasure / measure / etc
- DihedralGroup shouldn’t be capitalised

Decide on suitable name for:

- Nonrotating attractors
- measure
- selfmeasure
- weights (context of probabilities)
- d (Hausdorff dimension)
- Similarity/similitude
- symmetries / invariant maps / etc
- AffineMap – need something which encompasses similarities and symmetries
- Quadrature rules, a consistent convention is needed

Docs

- diameter
- translation
- Bounding balls

Here’s a list of all the easy and hard things I need to do, hopefully in some kind of order.

Bugs

- Discretisation of identity matrix errors when L_h={\[0\]}

General

- Global variable names

  - Cq (barycentric mindset!)
  - Gmres accuracy / threshold

- Check precompilation time doesn’t explode after modifying plots, or anything else

- Turn off macros via a flag? Would be useful to turn off parallelism

- Is it weird that the basis element has ‘measure’ as a field? Should this be ‘supp’ instead, and the measure is deduced as the restriction of the measure of the operator?

- Consider adding HomogenousHausdorffMeasure type. This would simplify sytax in a few places, and avoid the use of lots of redundant types in function defs.

- Nonrotating attractors – are they worth including as their own type?

- Should we have a Lebesgue measure type? If so, I think it makes sense to have a ‘non-empty-interior’ attractor, and Lebesgue must have this as its attractor. A key difference also is that Lebesgue has Int64 dimension, it is not a type parameter.

- Error catching for s-energy based on (s,d)

- Warn when someone tries to solve a problem on a fractal with d ≤ n-2 in ambient dimension n, e.g. a 1D fractal in 3D. The set has zero capacity, so the solution is always zero

- Consider removing barycentres from measure properties? It’s only fundamental in our opinion, as the creators of the barycentre rule…

- The bounding ball problem. Is not true in general, but I think the constant can be improved from a half to something much bigger.

- N-dimensional Laplacian, why not define the operator too?

- Should grading alg output a mesh? Can be misleading as if you then apply Hausdorff measure to each, they will all have measure one

- Laplace potential doesn’t work, at least in case where ambient dimensions don’t agree

- Composing similarities of different type, e.g. integer translation

- Fractal Unions need some consistency, in terms of numeric type and dimension. Adjust to Fractal{N,T}

- add collocation Laplace test on \[0,1\] to test suite

- Identity matrix should be sparse / diagonal

- Need a ‘getmeasure’ feature, which defaults to the right suppmeasure value for known cases

- Check for type stability of Galerkin & Collocation discretisations, and anything else time consuming

Vector indexing

- s-energy, and specifically the bit of it that gets called in the Galerkin discretisation, are still using trad indexing. This is OK, as it doesn’t use up memory, but conversion is inelegant and perhaps slightly innefficient.
- Overlfow of integer happens when you go like 60 levels deep

Inner product and discretisation

- Make new type ‘NonRotatingAttractor’, instead of current awkward implicit typing def
- Define inner product for inhomogenous attractors (Attractor / AbstractAttractor)
- Define generic routines which can take in a simple quadrature rule.
- Choose defaults to be different in 1D / 2D+
- Exploit symmetry of attractor
- Exploit the homogenous attractors, in terms of repeated values only, in discretise
- Currently the bounding ball is created O(N^2) times in the inner_product constructor
- Constructing the basis actually takes a little while? After subdivide_indices, which is still fast

Affine maps

- Add all the symmetry stuff from the previous version
- Add 3D stuff too, look into ways of taking cartesian products of groups

Simple (in principle) ports

- Fractal presets

  - Agree on a more sensible naming convention
  - Ignore weights, do at attractor level
  - Consider own submodule

- Gauss quadrature

- Operator presets (also consider own submodule)

- diameter function

- use the ‘cis’ function for exp(i\*x) computations in 3d kernel

Sugar

- Attractor constructor

- Measure constructor

- Basic operators on fractals

  - addition / translation
  - avoid ambiguous scaling
  - rotation

- IO for attractors, measures, bases

- Permit Attractor to take in n-tuple of similarities

Plotting

- Decide on a scheme to choose marker size
- For the boundary plot, look into importing a C polygonal clipping library: <https://github.com/rickbrew/GeneralPolygonClipper> using the instructions at <https://docs.julialang.org/en/v1/manual/calling-c-and-fortran-code/> having looked at this in more detail – there appears to be bugs and it does appear to be prone to errors :/ This website looks decent? <https://rosettacode.org/wiki/Sutherland-Hodgman_polygon_clipping#Julia>

Documentation

- Start with docstrings… when variable names are set
- Split into another layer of subsections

Testing

- Use Andrea’s Matlab
- Could import IFSintegrals into Test environment
- Test with/without fractal matrix exploitation

Mose

- The similarity composition issue/
- How do sub-environments work, if they are a thing? Is it just based on folder structure?
