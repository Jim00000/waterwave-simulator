# waterwave-simulator

## Description 

This is the term project for parallel programming course and our team implements the waterwave simulator

## Theory and implementation

In this project, we use this equation to implement waterwave
simulation

![$\frac{\partial^2{u}}{\partial{t^2}} + k\frac{\partial{u}}{\partial{t}} = c^2(\frac{\partial^2{u}}{\partial{x^2}} + \frac{\partial^2{u}}{\partial{y^2}})$](imgs/formula.png)

where parameters are

- u is the abbreviation of u(x, y, t)
- k is the damping constant
- c is the wave speed

and then we use a few parallel techniques to accelerate this program

For more details, please refer [Wiki-The-theory-and-wave-equation](https://github.com/Jim00000/waterwave-simulator/wiki/The-theory-and-wave-equation)

## Build 

Make sure that you have installed python2 and pip2

Carefully, python3 is incompatible

At first, install vtk, qt4, Boost library

`$ sudo apt-get install python-vtk python-qt4 python-qt4-dev libboost-dev-all`

`$ pip2 install -r requirements.txt` to install necessary libraries

`$ sh build.sh` to build the shared library for wave_equation_api

## Demo

![sequential_gif](https://i.imgur.com/jKV2ShV.gif)

In fact, it will work faster and more smoothly than the gif presented

## Authors

- [Jim00000](https://github.com/Jim00000)
- [ActKz](https://github.com/ActKz)
- [pityYo](https://github.com/pityYo)

## References

- [Y. Nishidate and G. P. Nikishkov, "Fast Water Animation Using the Wave Equation with Damping"](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.98.3126&rep=rep1&type=pdf)

## License

Copyright (C) 2017 the team of Jim00000, ActKz and pityYo

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.