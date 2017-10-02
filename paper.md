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

Three-dimensional structural models of macromolecules are used to gain
insights into biological processes. Most of the macromolecular structures
are determined using X-ray crystallography, which provides information
about electron density in a crystal. The electron density map is used to
build a model and can be later used to check the local quality of the model.

UglyMol is aimed at crystallographers who inspect electron density
at various stages of structure solution and model completion.
It can be also used by researchers who before using a model from
the Protein Data Bank want to check how well the model is supported
by the experimental data.

To make UglyMol easy to use by its audience,
the user interface is closely resembling Coot [@coot],
a desktop program popular among crystallographers.

Originally, UglyMol was developed to present results from the refinement
pipeline Dimple [@dimple]. Currently, it has also other uses.
It is included in at least five web applications:

* SynchWeb [@synchweb] in Diamond Light Source,
* EXI [@exi] in European Synchrotron Radiation Facility,
* CCP4 web services [@ccp4web],
* molstack [@molstack]
* and ContaMiner [@contaminer].

# References
