# NumSim
A flow solver for the Navier-Stokes equations with some nice features such as arbitrary geometries, heat transport and interaction with solids from scratch.
The project consists of four parts:
1. basic flow solver in a unit square with simple boundary conditions (2D driven cavity)
2. parallelization of the basic flow solver using MPI for distributed memory system
3. more complex simulations with an open source simulation software
4. project work with free choice of topics

## Installation & Prerequisites

You need
- gcc
- libvtk7.1(p) and
- cmake
- paraview

Executing ``` cd build && cmake .. & make -j && make install ``` in the main project folder compiles the sources and copys the executable to the build directory.
The first time the command above will fail, because there is no build  directory. Simply enter ``` mkdir build ``` to create the directory.

## Usage

- Run the programm: ``` ./numsim ../ini/your_szenario.txt ```

## Visualization

The results of the simulation (velocity u, velocity v, pressure p) will be saved for every timestep in the folder ```/out``` and can be visualized with paraview.

## Documentation
Run ```doxygen ./docs/doxyfile``` in the project folder to create the html and latex documentation. You can access the main page of the html documentation by opening ```docs/html/index.html``` with a browser of your choice.

## Running the tests

``` cd build && cmake .. & make -j && ctest -r ```

## Styleguide

This project follows the [C++ core guidelines](https://github.com/isocpp/CppCoreGuidelines).

## Authors

* **Author 1** - *tba* - [janiswissinger](https://github.com/janiswissinger)

* **Author 2** - *tba* - [kimkroener](https://github.com/kimkroener)

* **Author 3** - *tba* - [magnusostertag](https://github.com/magnusostertag)
