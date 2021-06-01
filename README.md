# LIGA
Code that implements the LIGA algorithm.

This is used to solve the structure of a material given only its atomic pair distribution function (PDF) as described here:

P. Juhas, D. M. Cherba, P. M. Duxbury, W. F. Punch, and S. J. L. Billinge. \Ab initio
determination of solid-state nanostructure". In: Nature 440.7084 (2006), pp. 655{658. doi:
10.1038/nature04556. url:
http://www.nature.com/nature/journal/v440/n7084/abs/nature04556.html

and in more boring detail, here:

P. Juhas, L. Granlund, P. M. Duxbury, W. F. Punch, and S. J. L. Billinge. \The Liga algorithm
for ab initio determination of nanostructure". In: Acta Crystallogr. A 64.6 (2008), pp. 631{640.
doi: 10.1107/S0108767308027591. url: http://dx.doi.org/10.1107/S0108767308027591

and tested on some real structures here:

P. Juhas, L. Granlund, S. R. Gujarathi, P. M. Duxbury, and S. J. L. Billinge. \Crystal structure
solution from experimentally determined atomic pair distribution functions". In: J. Appl.
Crystallogr. 42.3 (2010), pp. 623{629. doi: 10.1107/S002188981000988X. url:
http://dx.doi.org/10.1107/S002188981000988X

The initial Liga version was tested on Bucky ball, generated Platonic solids,
and Lennard-Jones clusters - some of these are available in a supplement
to the algorithm paper in Acta Cryst A, https://doi.org/10.1107/S0108767308027591.
The small organic molecules from various sources the authors have tested with a mixed success can be accessed through u24. Those structures and their distance
data should be under /u24/juhas/research/LIGA_runs/organics.