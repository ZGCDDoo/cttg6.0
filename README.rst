.. image:: https://travis-ci.com/ZGCDDoo/cttg6.0.svg?token=qU5DPfuxBM537pQaUFgc&branch=master
    :target: https://travis-ci.com/ZGCDDoo/cttg6.0
    
==========================================================================
 CTTG6.0 : Continuous time Tremblay Group
==========================================================================
<<<<<<< HEAD

=======
    
>>>>>>> master
:Authors: Charles-David Hébert, Maxime Charlebois, Patrick Sémon 
:Date: $Date: 2018-09-07 $
:Revision: $Revision: 6.0.0 $
:Description: Description

.. image:: https://travis-ci.com/ZGCDDoo/cttg6.0.svg?token=qU5DPfuxBM537pQaUFgc&branch=master
    :target: https://travis-ci.com/ZGCDDoo/cttg6.0
    
    
Graham
-------

g++
^^^^^^

* module reset 
* module load nixpkgs/16.09  gcc/5.4.0 armadillo boost-mpi
* mkdir build && cd build && cmake -DHOME=OFF -DGRAHAM=ON -DMPI_BUILD=ON .. && make

icpc (mpic++)
^^^^^^^^^^^^^^
* module reset
* module load intel/2017.5 armadillo boost-mpi
* mkdir build && cd build && cmake -DHOME=OFF -DGRAHAM=ON -DMPI_BUILD=ON .. && make

Mp2
------

g++ and icpc (mpic++)
^^^^^^^^^^^^^^^^^^^^^^
* module reset
* module load cmake/3.6.1  gcc/6.1.0  intel64/17.4  boost64/1.65.1_intel17 openmpi/1.8.4_intel17  armadillo/8.300.0
* mkdir build && cd build && cmake -DHOME=OFF -DMP2=ON -DMPI_BUILD=ON .. && make



    
