# Inverse Kinematic for Serial 6R/P robots

This is a C/C++ dynamic library for computing the inverse kinematics of serial 6R/P robots. It is based on the original 
executables created by 
[Arthur Angerer](https://www.umit-tirol.at/page.cfm?vpath=personendetails&genericpageid=100151&switchLocale=de_AT). The algorithm used is 
based on 
the [Dissertation of Martin Pfurner](http://repository.uibk.ac.at/alo?objid=1015078) and the accompanying paper [Husty 
M. et. al., A new and efficient algorithm for the inverse kinematics of a general serial 6R 
manipulator](https://doi.org/10.1016/j.mechmachtheory.2006.02.001). My former PhD student, S.M. Manongsong, and I 
extended the algorithm to be used for serial robots having also prismatic joints. The generalized form of the 
algorithm is published in [Capco J. et. al., Implementing HuPf algorithm for the inverse kinematics of general 6R/P manipulators](https://doi.org/10.1007/978-3-030-26831-2_6) and is implemented in this library.

This library has C-exports for use with Python (See the Python section for more information).

## Compiling and Building

**Disclaimer:** All MSVC project files are for VS2008 and setup for Release. You have to manually setup other build configurations if you want them. The project files should be forward-compatible to newer MSVC. 

**Compiling with GCC** 

Sample makefile are available in the gcc folder. One for a 32-bit build and one for a 64-bit build. This should work for both Windows and Linux operating systems.

## Prebuilt binaries (Windows)

1. Dynamically shared libraries are available for Windows and can be found in the `./bin` folder
1. For these binaries, the [Microsoft Visual C++ 2008 Redistributable Package](https://www.microsoft.com/en-us/download/details.aspx?id=26368) is needed.


## Python

A Python package that calls the `libhupf` library is available. This is available for both 32-bit and 64-bit versions of Python 2.x and 3.x. Users can install this in python using `pip`:

```
pip install libhupf
```

For Windows users `pip` will get the prebuilt binary wheel for Python 2.x and Python 3.x. Both 32bit and 64bit Python installations are supported. Windows users will need the [Microsoft Visual C++ 2008 Redistributable Package](https://www.microsoft.com/en-us/download/details.aspx?id=26368) .

For Linux/macOS users `pip` will get the source and compile using `gcc`.

Here is a sample Python code using the `libhupf` package
```
from libhupf import iksolver

#DH-parameters of a UR5-Robot
a = [0.0,-0.425,-0.39225,0.0,0.0,0.0]
#offset
d=[0.089159,0.0,0.0,0.10915,0.09465,0.0]
#rotation offset
theta = [0.0,0.0,0.0,0.0,0.0,0.0]
#twists
alpha = [1.5708,0.0,0.0,1.5708,-1.5708,0.0]
#R or P
rots = [1,1,1,1,1,1]   
 
iks = iksolver(a,d,theta, alpha, rots)

#flat end-effector
ee = [-0.582564,0.694272,0.422618,-0.333226,-0.271654,0.323744,-0.906308,-0.275819,-0.766044,-0.642788,0,-0.489578,0,0,0,1]

#solutions
iks.solve(ee)

```
## Authors
A. Angerer, J. Capco

## Acknowledgement

Many thanks to:
Arthur Angerer, Manfred Husty, Saraleen Mae Manongsong, Martin Pfurner and Hans-Peter Schröcker.

## Citing
 Link to Zenodo DOI: Coming soon.
 