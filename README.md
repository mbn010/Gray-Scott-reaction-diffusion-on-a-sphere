# Project Title

 Gray-Scott reaction diffusion on a sphere in Matlab

## Description

Reaction-diffusion models are notable for producing many patterns found in nature. In this short project I put together a quick and simple surface finite element method which I've used to solve the Gray-Scott reaction-diffusion equations on the surface of a sphere given by gamma.

<p align="center"><img src="https://user-images.githubusercontent.com/17126595/50368113-549ae900-057d-11e9-9510-933ffef61c6f.png" width="250" alt="domain"/></p>


The model itself is described in the following figure.

<p align="center"><img src="https://user-images.githubusercontent.com/17126595/50368112-549ae900-057d-11e9-85db-464ab4caaae8.png" width="250" alt="gray scott model"/></p>

And the parameters I've selected to model are given in:

<p align="center"><img src="https://user-images.githubusercontent.com/17126595/50368114-55337f80-057d-11e9-880e-b23dedfa18fa.png" width="130" alt="params"/></p>

This code makes use of the Matlab vectorization notation and iterative solvers to produce a high performance, and colourful, output given in the animation.

<p align="center">
<img src="https://user-images.githubusercontent.com/17126595/50352623-b20f4580-053d-11e9-8a04-0c6b055fc0c5.gif" alt="animation" />
  </p>

These animations correspond to two separate initial conditions which are available from the code.

## Deployment

To demo this application, execute the .m script file in a Matlab environment.

## Built With

* [Matlab](https://www.mathworks.com/products/matlab.html) - Matlab
* [Gray Scott](https://groups.csail.mit.edu/mac/projects/amorphous/GrayScott/) - Further details of the Gray-Scott model

## Versioning

Version 1.0 only.

## Authors

* **Michael N** - *Initial work* - [mbn010](https://github.com/mbn010)




