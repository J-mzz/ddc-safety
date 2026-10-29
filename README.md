# Robust Data-Driven Safe Control Using Density Functions

MATLAB code for the paper

> J. Zheng, T. Dai, J. Miller, and M. Sznaier, "Robust data-driven safe control using density functions," *IEEE Control Systems Letters*, vol. 7, pp. 2611–2616, 2023.
> [doi:10.1109/LCSYS.2023.3287801](https://doi.org/10.1109/LCSYS.2023.3287801) · [arXiv:2303.09004](https://arxiv.org/abs/2303.09004)

The method starts from noisy measurements of an unknown polynomial system and solves a sum-of-squares program for a rational state feedback `u = psi/rho` and a density function `rho`. Together they certify that trajectories starting in the initial set `X0` avoid the unsafe set `Xu`. The certificate holds for every system consistent with the data and under bounded process noise.

## Requirements

- MATLAB (tested with R2024a) with the Symbolic Math, Optimization, and Statistics and Machine Learning Toolboxes
- [YALMIP](https://yalmip.github.io/)
- [MOSEK](https://www.mosek.com/) (free academic licenses are available); tested with MOSEK 9.1

## Usage

To draw the panels of the paper's figures from the saved results in `data/`, run:

```matlab
plot_example1('fig1b')   % Example 1: 'fig1a', 'fig1b', 'fig2a', 'fig2b'
plot_example2('fig3b')   % Example 2: 'fig3a', 'fig3b'
```

To recompute the results, set `panel` at the top of a script and run it:

```matlab
example1_flow    % Example 1, flow system (Figs. 1 and 2)
example2_twist   % Example 2, twist system (Fig. 3)
```

| `panel` | Paper | Controller | Simulation |
|---|---|---|---|
| `fig1a` | Fig. 1(a) | non-robust | open loop, no process noise |
| `fig1b` | Fig. 1(b) | robust, designed for `‖w‖∞ ≤ 2` | closed loop, `‖w‖∞ ≤ 2` |
| `fig2a` | Fig. 2(a) | non-robust (designed for `w = 0`) | closed loop, no process noise |
| `fig2b` | Fig. 2(b) | non-robust | closed loop, `‖w‖∞ ≤ 2` (100 trajectories, zoomed in) |
| `fig3a` | Fig. 3(a) | robust, designed for `‖w‖∞ ≤ 1` | open loop, no process noise |
| `fig3b` | Fig. 3(b) | robust, designed for `‖w‖∞ ≤ 1` | closed loop, `‖w‖∞ ≤ 1` |

The scripts plot the results directly. To draw a new run in the style of the paper, set `save_file` in the script and call `plot_example1(panel, save_file)` or `plot_example2(panel, save_file)`.

The fixed random seeds reproduce the data exactly, and the controllers match the saved ones to within 0.5%. The simulated trajectories use randomly switching noise and `ode15s`, so they can differ slightly from the saved results, which were computed with MATLAB R2020b. In Fig. 2(b) the trajectories stop where |u| reaches 10^4, as stated in the paper; the saved data of `fig2b` are cut there, while a new run keeps the full trajectories.

## Repository layout

```
example1_flow.m, example2_twist.m   scripts for Examples 1 and 2
plot_example1.m, plot_example2.m    draw the panels from saved results
data/                               saved results of each panel (fig1a_*.mat, ...)
figures/                            PNG images of the panels
utils/
  ddc_safety_robust.m      SOS program for rho and psi (Example 1)
  simulateTrajectories.m   closed/open-loop simulation of Example 1
  getDataRND.m, getCons.m  noisy data and the consistency set
  nontrivial_constraints.m removes redundant constraints of a polytope
  getSystem.m, setCons.m   system of Example 1 and problem dimensions
  third_party/             external code, see utils/third_party/README.md
```

## Citation

If you use this code, please cite:

```bibtex
@article{zheng2023robust,
  title={Robust data-driven safe control using density functions},
  author={Zheng, Jian and Dai, Tianyu and Miller, Jared and Sznaier, Mario},
  journal={IEEE Control Systems Letters},
  volume={7},
  pages={2611--2616},
  year={2023},
  publisher={IEEE},
  doi={10.1109/LCSYS.2023.3287801}
}
```

## License

The code in this repository is released under the [MIT License](LICENSE), except the files in `utils/third_party/`, which keep their original licenses.
