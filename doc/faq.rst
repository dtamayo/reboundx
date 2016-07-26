.. _faq:

Frequently Asked Questions
==========================

Here are some answers to frequently asked questions, as well as not so frequently asked questions that are known issues.
If you have a question that is not answered here (likely while this is built up!), please open an issue on github.  
In principle this makes it possible for others to find it with a google search.
To do this, go to https://github.com/dtamayo/reboundx, go to the `Issues` tab along the top, and click the green `New Issue` button.

.. contents::
    :local:
    :depth: 2

.. _kerneldied:

I get a Kernel has died message in iPython notebook
---------------------------------------------------

These problems are typically due to mismatched installations of REBOUND and REBOUNDx.  
Version conflicts among python modules are nicely resolved by virtual environments (or conda environments if you use the Anaconda distribution of python).  
If you run into this issue, please follow the steps in 5.1.1 of http://rebound.readthedocs.org/en/latest/python_quickstart.html.  
Then activate your newly created environment (as described in 5.1.1) and::
    
    pip install numpy matplotlib

or::

    conda install numpy matplotlib

if you're using Anaconda (and any other packages you like).
Then follow 5.1.3 to install REBOUND if you'd like the source code and examples, or 5.1.2 for the quick installation.
Finally, follow the installation instructions for REBOUNDx (with your environment activated) at :ref:`python_quickstart`.

If this does not resolve your problem, please open an issue on github including the output from the terminal you used to launch the ipython notebook.

When I change one particle's parameter, it also changes another particle's value
--------------------------------------------------------------------------------

This is likely because you assigned one particle to another, e.g., ::
    
    sim = rebound.Simulation()
    rebx = reboundx.Extras(sim)

    sim.add(m=1.)
    sim.particles[0].params["tau_a"] = 3.

    sim.add(sim.particles[0])

Here we added ``particles[0]`` to the simulation again, so the simulation now has two copies, i.e., ``particles[1] == particles[0]``.
However, by default copies are shallow, i.e., they don't copy memory pointed to by pointers, which is how the particle parameters are stored.
This means that the two particles share the same memory, so when we change one, it changes the other too.
