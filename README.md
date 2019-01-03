# MABC
<!-- Created by Wim Delva, 27 November 2018 -->
Model calibration as a missing data problem: Approximate Bayesian Computation
with MICE (Multivariate Imputation by Chained Equations)


MABC is an R package for calibrating stochastic simulation models to data. It
implements a sequential Approximate Bayesian Computation method that employs
MICE (Multivariate Imputation by Chained Equations) as an emulator to link model
inputs to model outputs.

## CONTENTS

This repository holds the source code for the MABC package:

* [Code files](#code-files)
* [System and software requirements](#system-and-software-requirements)
* [Copyright and licensing information](#copyright-and-licensing-information)
* [Contact information](#contact-information)


## CODE FILES 

All code is written in R. R is a statistical programming language and software
package that is distributed under a GNU General Public License. R documentation
and software is available for free download through the R Project for
Statistical Computing website at http://www.r-project.org. The software is
available as an executable file for a wide variety of operating systems and
computer architectures, and a compilable binary is also available should you
need it.

MABC.R -- The core function to run MABC.

mice.parallel, model.mpi.run.R, model.parallelrun.R -- Helper functions used by
the above two functions.
  

## SYSTEM AND SOFTWARE REQUIREMENTS

### Operating system

  We have only tested this code on personal computers (OS X Version 10.11.6) and
  the golett cluster of the Flemish Supercomputer Centre (VSC).

### Required software

  Open MPI or Microsoft MPI Redistributable Package, and Rmpi. For installation instructions, see <http://fisher.stats.uwo.ca/faculty/yu/Rmpi/>.
  
  R version >= 3.4.0
  
  In addition, a long list of auxiliary R packages is required to run the
  functions in this package. All of these dependencies will get installed
  automatically when installing the software.
  
  install.packages("devtools")
  
  library(devtools)
  
  devtools::install_github(wdelva/MABC) 

## COPYRIGHT AND LICENSING INFORMATION

All files are copyright protected and are made available under the GPL 3.0
License <https://www.gnu.org/licenses/gpl-3.0.en.html>. This means that this
work is suitable for commercial use, that licensees can modify the work, that
they must release the source alongside with Derivative Work, and that Derivative
Work must be released under the same terms.

## CONTACT INFORMATION

Please contact Prof Wim Delva with questions regarding the files provided and
their authorized usage:

Wim Delva Email: <DELVAW@sun.ac.za>
