.. pySCATMECH documentation master file, created by
   sphinx-quickstart on Tue Mar 17 07:29:04 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. raw:: html
	 
   <table width="100%" border="0" cellpadding="0" cellspacing="0" align="left">
     <tr>
       <td>
          <img src="_static/head_otd.jpg" width="750" height="68" border="0" USEMAP="#nav_otd" alt="Sensor Science Division">
       </td>
     </tr>
   </table>
   <br clear=all>
   <MAP NAME="nav_otd">
       <AREA SHAPE="RECT" COORDS="15, 5, 150, 20" HREF="http://physics.nist.gov" TARGET="TOP" ALT="Physical Measurement Laboratory Home">
       <AREA SHAPE="RECT" COORDS="15, 20, 290, 50" HREF="http://www.nist.gov/pml/div685/index.cfm" TARGET="TOP" ALT="Sensor Science Division Home">
       <AREA SHAPE="RECT" COORDS="590, 5, 725, 50" HREF="http://www.nist.gov" TARGET="TOP" ALT="NIST Home">
    </MAP>


pySCATMECH : A Python Interface to the SCATMECH Library
=======================================================

pySCATMECH is a Python interface to the SCATMECH library of scattering codes.
The documentation for the 
specific models are available at `SCATMECH Polarized Light
Scattering C++ Class Library <https://pages.nist.gov/SCATMECH/docs/index.htm>`_.


:Author: Thomas A. Germer
:E-mail: thomas.germer@nist.gov
:Telephone: 301-975-2876
:Address:
   | Sensor Science Division
   | National Institute of Standards Technology
   | 100 Bureau Drive Stop 8443
   | Gaitherburg, MD 20899-8443
:Current Version: 0.0.2.1 (8 July 2020)
:Requires: Python 3.6 or greater
:To install: Use ``pip install pySCATMECH`` at a command prompt.
	     
.. toctree::
   :maxdepth: 2
   :caption: Examples

   UsingMueller
   UsingFresnel
   UsingBRDF
   UsingLocal
   UsingIntegrate
   UsingScatterer
   UsingRCW
   UsingCrossRCW
	   
.. toctree::
   :maxdepth: 2
   :caption: Modules
	     
   mueller
   fresnel
   model
   brdf
   local
   integrate
   scatterer
   rcw
   crossrcw

   
Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
