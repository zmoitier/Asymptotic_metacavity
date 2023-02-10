[![DOI](https://zenodo.org/badge/357641950.svg)](https://zenodo.org/badge/latestdoi/357641950)

# Asymptotic_metacavity

## Reference

This is the code associated with the article:

- C. Carvalho and Z. Moitier, _Scattering resonances in unbounded transmission problems with sign-changing coefficient_, IMA Journal of Applied Mathematics. [[arXiv](https://arxiv.org/abs/2010.07583), [HAL](https://hal.science/hal-02965993), [DOI](https://doi.org/10.1093/imamat/hxad005)]

## Requirements

- Python version:

  - Tested on Python 3.10;
  - Should work on Python 3.7 but not tested.

- Require the following libraries:

  - For numerical computations: [cxroots](https://github.com/rparini/cxroots), [Numba](https://github.com/numba/numba), [NumPy](https://github.com/numpy/numpy), and [SciPy](https://github.com/scipy/scipy);

  - For symbolic computations: [SymPy](https://github.com/sympy/sympy);

  - For visualization: [Matplotlib](https://github.com/matplotlib/matplotlib);

  - For Jupyter notebook: [IPython](https://github.com/ipython/ipython), [ipywidgets](https://github.com/jupyter-widgets/ipywidgets), and [JupyterLab](https://github.com/jupyterlab/jupyterlab).

Might works with previous versions of the libraries but if it does not works try to update the libraries for example through pip

```bash
python3 -m pip install --user --upgrade -r requirements.txt
```

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

### Data

The folder `data` contains computed resonances and modes for the circular metamaterial cavity. The data can be obtained by executing the codes in the folder `Compute resonances`.

### Resonances

Run the Jupyter notebook `Interactive resonant modes for circular negative metamaterial cavity .ipynb` to get the resonances sets and associated resonant modes for the circular cavity. This interactive notebook allows you to recover Figures 3, 8, 9 from the paper.

### Scattering

Run the Jupyter notebook `Interactive scattering for circular negative metamaterial cavity .ipynb` to get the scattering solution and stability constant for the circular cavity. This interactive notebook allows you to recover Figures 1, 2 from the paper.

### Asymptotic

Run the Jupyter notebook `Quasi resonances.ipynb` to compute the quasi-resonances for a given metamaterial cavity. This notebook allows you to recover quasi-resonances used in Figure 7 from the paper.

## Contact

If you have any questions or suggestions please feel free to create [an issue in this repository](https://github.com/zmoitier/Asymptotic_metacavity/issues/new).
