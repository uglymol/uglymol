---
title: 'UglyMol: a WebGL macromolecular viewer focused on the electron density'
tags:
 - protein structure visualization
 - electron density
 - WebGL
authors:
 - name: Marcin Wojdyr
   orcid: 0000-0003-3980-4092
   affiliation: 1,2
affiliations:
 - name: Diamond Light Source Ltd, Harwell Campus, Didcot, OX11 0DE, UK
   index: 1
 - name: Global Phasing Ltd, Sheraton House, Cambridge, CB3 0AX, UK
   index: 2
date: 18 July 2017
bibliography: paper.bib
---

# Summary

UglyMol [@uglymol] is a macromolecular viewer specialized in presenting
macromolecular models together with the electron density.
It uses web technologies (JavaScript and WebGL) and is suitable for embedding
in web applications. The project was started as a fork of xtal.js [@xtaljs].

UglyMol is aimed at crystallographers who inspect electron density
at various stages of structure solution and model completion.
It can be also used by researchers who before using a model from
the Protein Data Bank want to check how well the model is supported by the experimental data.

To make UglyMol easy to use by its audience,
the user interface is closely resembling Coot [@coot],
a desktop program that is popular among crystallographers.

UglyMol is currently used in a few web applications.
In SynchWeb [@synchweb] in Diamond Light Source it presents results from
the refinement pipeline Dimple [@dimple].
In EXI [@exi] in European Synchrotron Radiation Facility it shows
results from a phasing pipeline.
And in CCP4 web services [@ccp4web] -- results from multiple programs
and pipelines.

# References
