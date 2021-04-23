[![DOI](https://zenodo.org/badge/357641950.svg)](https://zenodo.org/badge/latestdoi/357641950)

# Asymptotic_metacavity

## Reference

This is the code associated with the article:

- C. Carvalho and Z. Moitier, _Asymptotics for metamaterial cavities and their effect on scattering_ [[arXiv](https://arxiv.org/abs/2010.07583), [HAL](https://hal.archives-ouvertes.fr/hal-02965993)]

## Requirements

- Python version:

  - Tested on Python 3.8;
  - Should works on Python 3.7 but not tested.

- Require the following libraries:

  - [Matplotlib](https://github.com/matplotlib/matplotlib),
  - [Numba](https://github.com/numba/numba),
  - [NumPy](https://github.com/numpy/numpy),
  - [SciPy](https://github.com/scipy/scipy),
  - [cxroots](https://rparini.github.io/cxroots/).

## Install

Clone from GitHub repository:

```bash
  git clone https://github.com/zmoitier/Asymptotic_metacavity.git
```

## Instructions for usage

### Compute resonances

The folder `Compute resonances` contains codes necessary to compute resonances using contour integration techniques. For convenience some cases have been pre-computed and stored in the `data` folder.

Run

```bash
     python3 calc_cxroots eps
```

to compute the 64 resonances for the circular cavity with permittivity eps.

### data

The folder `data` contains computed resonances and modes for the circular metamaterial cavity. The data can be obtained by excuting the codes in the folder `Compute resonances`.

### Resonances

Run the Jupyter notebook `Interactive resonant modes for circular negative metamaterial cavity .ipynb` to get the resonances sets and associated resonant modes for the circular cavity. This interactive notebook allows you to recover Figures 3, 8, 9 from the paper.

### Scattering

Run the Jupyter notebook `Interactive scattering for circular negative metamaterial cavity .ipynb` to get the scattering solution and stability constant for the circular cavity. This interactive notebook allows you to recover Figures 1, 2 from the paper.

### Asymptotics

Run the Jupyter notebook `Quasi resonances.ipynb` to compute the quasi-resonances for a given metamaterial cavity. This notebook allows you to recover quasi-resonances used in Figure 7 from the paper.
