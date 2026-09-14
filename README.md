# CADET Projects

## Motivation

This repo is used to archive my files for [CADET](https://github.com/cadet) simulation (chromatography simulations).

It will be used for both :

- [CADET-Python](https://cadet.github.io/v5.1.X/developer_guide/cadet_python.html)
- [CADET-Process](https://cadet-process.readthedocs.io/en/latest/index.html)

Both are front-ends for [CADET-Core](https://cadet.github.io/v5.1.X/index.html) :

Cadet-python is more oriented for developers or users wanting latest models/isotherms, setting everything is a must.

Cadet-process is a more guided one, selecting some parameters by default and making the code more readable overall (it is object oriented).

One should need both even if Cadet-python is enough for experimented users.

____

## Structure

In this repo, you can find :

- both existing tutorial for the front-ends,
- guide to install adequate python env ([anaconda](https://www.anaconda.com/)) and install latest CADET-Core
- Code for the projects
- Some chromatograms

```
|
|
|-- cadet-python ---|
                    |
                    |-- Install           : .yml + .sh
                    |-- Tuto              : cloned
|
|
|---- Projects ----|
                   |
                   |-- Graphs      : .yml + .sh
                   |-- H5 files    : cloned
                   |-- Code        : python codes for the projects
|
|
|------ Data ------|
                   |
                   |-- Chromatogram
                   |-- BT measurements
|
|
|- cadet-process --|
                   |
                   |-- Install  : .yml + .sh
                   |-- Tuto     : cloned
```
