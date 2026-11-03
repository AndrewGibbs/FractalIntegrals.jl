# Next steps

This document is written by Andrew, for Josh, explaining future directions we can get Codex and Claude to work on.

## How we work

I think first we should both share our ideas, and check we aren't stepping on each others' toes with these. Then each create a new branch on GitHub (Codex can talk you through this) and pull to our local machine. These can be `Andrew-dev` and `Josh-dev` if you like. We can both merge into main periodically, when new features are complete. 

### Testing

This is something I am happy to hand over to agents fully. Every change must have associated unit tests added to the `test` folder - but the AIs can decide what is best here - they are much more thorough than I am with this.

More generally, I will get Claude to create some new unit tests on the main branch, so these are also in our new branches.

AI is usually very good at checking that all tests pass before finalising any changes.

### Some ground rules for the agents

I've worked on several projects with AI now, and I'm wary of a couple of pitfalls.

Something that can easily happen is the AI runs away with the project, churning out ideas faster than the humans have time to understand them, and it is too tempting to push on to see how effective the end result will be *before* understanding how we got here. This can happen with mathematical derivations and code.

I like understanding everything that is in the codebase right now. AI has a tendancy to write far more code than is needed, and I am concerned it will change my code significantly for a few performance gains, to the extent where I won't understand it at all. I accept that we want AI to write new stuff, and I can accept that I won't understand or necessarily even read the code for that, but we should be mindful of keeping the AI contributions separate.

**Ground rules for AIs**
* New features live in new subfolders inside of `src`.
* Changes to existing code must only be made if they are essential for compatibility with new features, and must be kept minimal. The AI should document these.
* Separately, existing code will be tested and optimised for performance, but changes to the code must be kept minimal.
* New features should be thoroughly documented, and the documentation should be read by a human. The mathematical derivation should also be stored somewhere, although this needn't be public.
* Official new version release should be authorised by humans.

## Summary of my suggested changes

### Sexy new features
* FFT-Toeplitz solver - need to think about how much information is needed for this. We want the user to provide the simplest possible representation of the fractal tiling, and the code to do the rest. And possibly, should we generalise the existing [Julia Toeplitz package](https://github.com/JuliaLinearAlgebra/ToeplitzMatrices.jl), which can't do block-Toeplitz? It doesn't make sense to have our own in-repo Toeplitz algebra package.
* A polynomial integrand calculator. Jeevon's MSc thesis has an exact formulat for this for a $\mathbf{p}$-degree polynomial.
* High order quadrature - Jeevon's method and Zois' code.
* H-matrix solver. Can use the bounding balls to determine the tree, certainly possible. If we go down the separable-kernel route, this could use the calculator above, although I think it is actually fewer FLOPs to use the cubature rule, which can be made exact if we choose enough points.
* Reptile type. Need to find a clean way of encoding adjacency though. If Josh works on this, it would be worth respecting the existing type hierarchy.
* Exploit symmetry fractals to repeat blocks of matrices.
* Modified Helmholtz and other radial kernels

### Things which I never got around to implementing
* The block operator logic for collocation. Similarly, of new features, which are compatible with collocation? For e.g. can we force translation invariance?
* Some other little stuff and bug fixes in `todo.odt`.
* Most of the documentation - although I'm undecided if I'm happy with AI writing this. This includes the stuff in `docs/` and also the docstrings in the code. Again, let's discuss.

### Things which should be optimised
* The $d=2$ plotting algorithm, is a bit random and needs a sensible restart condition. Relies on a dodgy C port which is poorly documented and we've had to trial and error the arguments, AI would be great at ironing this out. Also if it can be made faster that would be great.
* Sometimes can take a long time to import. I'd like to know why this is - I've never gotten to the bottom of it.
* The near and far field potential operators, these should be parallelised.
* All functions should be reviewed by AI, checked for efficiency, and rewritten *for efficiency purposes only*, keeping changes  minimal otherwise.

### Some questions about the direction of the project
* Should we include more general basis functions?
* Can we generalise to non-radial kernels?

## Responses to Josh's ideas
* It looks like many of our ideas, high order quadrature, new kernels, H-matrices, are shared. 
* Variable contrast - have you worked this through mathematically? I think we'd need a more general quadrature for this, graded, which we implemented already but only in FractalIntegrals. There may be other things I'm missing.
* All of the internal operations assume an IFS attractor, or a union of these. Other types of fractals such as Julia sets do not fit within this framework - and it might make sense for them to be in a separate code base.
* A lot of what you're asking can already be done - but my documentation is the weakness. If you want to learn about large scale codebases then it is worth rebuilding from scratch, but if you want to add these new features quickly, it makes sense to tweak what I've done and add new features.