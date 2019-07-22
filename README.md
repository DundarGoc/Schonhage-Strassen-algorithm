# A benchmark of the Schönhage-Strassen algorithm.

This is my Master's Thesis. It is a benchmarking of polynomial multiplication modulo an odd number using the [Schönage-Strassen algorithm](https://en.wikipedia.org/wiki/Sch%C3%B6nhage%E2%80%93Strassen_algorithm). The benchmarking is done against the multiplication function `_nmod_poly_mul` in the [FLINT library](https://github.com/wbhart/flint2). The code also requires FLINT to run. The code was developed by using the [zn_poly](https://gitlab.com/sagemath/zn_poly/) library as the starting point.

Note that this isn't an implementation of a "pure" Schönhage-Strassen. I wouldn't recommend using this to try to learn it. The code also uses Nussbaumer multiplication, Truncated Fourier Transform, Kronecker substitution, Montgomery reductions and much more.

## Can I use it for my computations?

I advice against it. There are already libraries such as the aforementioned [FLINT](https://github.com/wbhart/flint2) and [GMP](https://gmplib.org/) that does this. The performance of this code is largely dependent on many tuning parameters. It's possible to achieve better results than FLINT if you use the right tuning parameters, but knowing what they are is impossible meaning that you'd have to guess what those values would be. The optimal tuning parameters are also entirely dependent on architecture; two different computers would need different parameters. The next step in the thesis would be to create a more consistent tuning but I unfortunately didn't have enough time.

If you still wish to give it a try the multiplication algorithm is in SSA.c.

## Questions?

If you have questions about the code then I'll do my best to answer them. If you have questions about the algorithms you're probably better off asking someone more knowledgeable since my understanding of the algorithms are shaky at best.
