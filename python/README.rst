
General Recognition Theory
==========================

This is a Python implementation of a subset of the tools used in General
Recognition Theory. It is based on the GRT toolboox for Matlab by Leola
Alfonso-Reese, hosted on `her website
<http://www-rohan.sdsu.edu/~leola/toolbox.html>`_. I was working on porting it
over to Python for my own research, and realized I should make my efforts
public. 

So far I have only been working on the General Linear Classifier. Alfonso-Reese
has also implemented the General Quadratic Classifier.

A note on changes in this version
---------------------------------

First, forgive the mess, I should be able to clean it up a bit in the next few
days.

It turns out the method used in the GRT toolbox does not always find the right
answer: using the built-in optimization toolkit  optimizer with a single
starting point, a variance of 10 and paramters derived from the Fisher
coefficient, is not always adequate to find the maximum likelihood solution.
I have found much better answers by first using the one-parameter optimization
suggested by [Ashby1992]_ (referencing earlier methods summarized by
[Fukunaga1972]_). Because minimization tends to get close to the best solution
but often doesn't find it, the algorithm here then does Nelder-Meade gradient
descent search on all parameters.

This optimization method is now much better than the Alfonso-Reese
method, but I can't guarantee that it always finds the best solution. I'm
still trying to work on better guarantees that the optimization will work
correctly.


.. [Ashby1992] Ashby, F. G. (1992). "Multidimensional Models of Categorization."
    In F. G. Ashby (Ed.), *Multidimensional Models of Perception and Cognition*
    (pp. 449-483). Hillsdale, NJ: Lawrence Erlbaum Associates.

.. [Fukunaga1972] Fukunaga, K. (1972).  *Introduction to Statistical Pattern
    Recognition*. San Diego, CA: Academic Press.

System requirements
-------------------
The code here requires Python 2.x (availible `here <http://python.org>`_) with
the numpy and scipy libraries (availible `here <http://numpy.org>`_).

Plotting also requires Matplotlib (available `here
<http://matplotlib.sourceforge.net/>`_).


Copyright Information
---------------------
All files are marked internally with copyright information. I have made
minimal changes to the Matlab files, focusing instead on porting the code over
to Python. Work by me (all the Python files) is copyrighted by me, but I assent
to copying for personal or academic reasons. I do not assent to any copying or
publication for any commercial reason without my express written consent. I
also make no claim that the work provide here is useful for any particular
purpose, and am not responsible for any ill that may befall you as a result of
using it.

If you are a copyright holder of any material on this site and are displeased
with the attribution provided here or wish that any material be taken down,
please contact me and I will fulfill your request as quickly as possible.



