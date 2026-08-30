# README

[![Python Version](https://img.shields.io/badge/python-3.10%2B-blue)](https://python.org)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE.md)

The goal of this repository is to provide a validation of the analytical solution for for the Geometric Brownian Motion Stochastic Differential Equation. Then, a from first principles discretization will be applied to derive the Euler-Maruyama method, which will be implemented, and its convergence tested. Finally, a practical application will be implemented (population growth under uncertainty).

## Quick start / Setup

To explore the project, clone or download the repository and open the notebooks in Jupyter or VS Code. The analysis is carried out in Python and makes use of the standard scientific stack, including NumPy, matplotlib, and Jupyter.

```bash
git clone https://github.com/pablo-rtz/sde_analysis.git
cd sde_analysis
```

Then open either of the notebooks in the [src/sde_analysis](src/sde_analysis) folder:

- [src/sde_analysis/Convergence_Test.ipynb](src/sde_analysis/Convergence_Test.ipynb)
- [src/sde_analysis/Modelling_test.ipynb](src/sde_analysis/Modelling_test.ipynb)

The full mathematical writeup is available in the [documentation](Docs/Documentation.md).

### Prerequisites

- Python 3.10 or later
- Jupyter Notebook / VS Code with Python support
- Numpy and matplotlib

## Project overview

This repository is structured around three connected goals:

1. Validate the analytical solution to the GBM SDE
2. Discretize it using Euler-Maruyama and validate the convergence orders numerically
3. Apply the model to real population dynamics with stochastic uncertainty

The core stochastic differential equation studied is

$$
dS_t = \mu S_t \, dt + \sigma S_t \, dW_t
$$

with exact solution

$$
S_t = S_0 \exp\left(\left(\mu - \frac{1}{2}\sigma^2\right)t + \sigma W_t\right).
$$

The repo then shows how this solution can be approximated numerically and how the error behaves under refinement of the time step.

## Data

The project includes historical population data for Ethiopia, obtained from Our World in Data / UN World Population Prospects 2024, stored in the [Data](Data) folder.

This data is used to estimate the drift and volatility parameters for the stochastic population model and to produce multi-path forecasts.

## License

This project is licensed under the GNU General Public License v3.0 or later (GPL-3.0-or-later).

See the [LICENSE](LICENSE.md) file for details.

If you use this repository, please cite it following the [citation](citation.cff).
