## Test environments
* windows-latest (on github, with INLA), R 4.1.2
* windows-latest (local, with INLA), R 4.1.0
* ubuntu 20.04 (on github, with INLA), R 4.1.2
* ubuntu 20.04 (on github, with INLA), R 4.1.2, R devel
* macOS-latest (local, with INLA), R 4.1.0

## Submission notes
The package contains several vignettes which are available at:

* <https://emi-sstcdapp.unavarra.es/bigDM/bigDM-1-fitting-spatial-models.html>
* <https://emi-sstcdapp.unavarra.es/bigDM/bigDM-2-parallel-and-distributed-modelling.html>
* <https://emi-sstcdapp.unavarra.es/bigDM/bigDM-3-fitting-spatio-temporal-models.html>

## R CMD check results
Comments:

* Almost all check result are 0 errors, 0 warnings and 0 notes.
* The last test environment (macOS) generates a code note because, the package holds 2 high-dimensional datasets with 5.1Mb.
