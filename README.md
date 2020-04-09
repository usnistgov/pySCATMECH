# pySCATMECH Package

pySCATMECH is a Python interface to SCATMECH: Polarized Light Scattering C++
Class Library.

## SCATMECH

SCATMECH is an object-oriented C++ class library developed to
distribute models for light scattering applications. Included in the
library are models for diffuse surface scattering that predict the
bidirectional reflectance distribution function (BRDF), codes for
calculating scattering by isolated particles, and codes for
reflection, transmission, and diffraction from gratings. Emphasis has
been given to those diffuse scatter models that are physics-based and
which predict the polarization properties of the scattered light. The
library also includes a number of classes that are useful for working
with polarized light or the optics of thin films. The library is
designed to enable expansion of new models.

See [https://pages.nist.gov/SCATMECH/index.htm](https://pages.nist.gov/SCATMECH/index.htm) for full SCATMECH documentation.

## pySCATMECH Modules

PySCATMECH contains nine modules:

* mueller - Tools for handling Mueller matrices, Stokes vectors, Jones matrices, and Jones vectors
* model - Tools for handling the SCATMECH::Model class, which handles generic models
* fresnel - Tools for handling optical functions, thin films, and reflection and transmission coefficients
* brdf - Tools for handling bidirectional reflectance distribution function (BRDF) models (SCATMECH::BRDF_Model)
* local - Tools for handling differential scattering cross-section (DSC) models of local defects on surfaces (SCATMECH::Local_BRDF_Model)
* rcw - Tools for handling rigorous couple wave (RCW) analysis of 1D periodic gratings (SCATMECH::RCW_Model)
* crossrcw - Tools for handling RCW analysis of 2D periodic gratings (SCATMECH::CrossRCW_Model)
* scatterer - Tools for handling free-space scattering functions (SCATMECH::FreeSpaceScatterer)
* integrate - Tools for integrating BRDF or DSC models

