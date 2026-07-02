# Metasurface Simulations (MATLAB)

MATLAB code supporting my review and reproduction of Lin, Fan, Hasman, and Brongersma, "Dielectric gradient metasurface optical elements," *Science* 345(6194), 2014. These scripts independently reproduce the paper's optical results and analyze the meta-atom design.

The Huygens-Fresnel propagation and all analysis code here are mine. The full-wave meta-atom sweeps that produce the input `.mat` data (`Tx_phix_data.mat`, `Ty_phiy_data.mat`) were run in Ansys Lumerical in collaboration with a colleague, as credited in the paper.

## Scripts

- **axicon.m** — Simulates the Bessel beam from an axicon by numerically propagating an incident Gaussian through the axicon's radial phase profile with the Huygens-Fresnel integral. Produces the xz and xy intensity profiles and integrates the inner spot against the outer rings to confirm the ~6:1 ratio that is the Bessel-beam signature.

- **metalens.m** — Simulates the metalens by applying the hyperboloidal lens phase profile to a plane wave over a circular aperture and propagating to the focal region. Recovers the near-diffraction-limited focal spot (~0.6 um FWHM) and the small focal shift.

- **numerical_simulation_analysis.m** — Processes the full-wave meta-atom sweep (nano-beam width 80 to 160 nm, period 160 to 240 nm at 550 nm). Computes co- and cross-polarization efficiency (CoPE, CrPE) and diffraction efficiency over the 2D design space, locates the minimum-CoPE design point, and compares it against the paper's 120 nm / 200 nm design.

- **numerical_transmission_phase.m** — Companion analysis of transmission and relative phase across the same width/period sweep.

- **transmission.m** — Plots the digitized experimental transmission coefficients and relative phase from the paper's supplementary data and computes CoPE/CrPE from them, showing the polarization-conversion optimum sits nearer 520 nm than the 550 nm design wavelength.

- **refractive_index.m** — Plots the real and imaginary refractive index of the deposited poly-Si film, digitized from the paper's supplementary figures.

## Results/

Generated figures (Bessel beam, metalens focus, axicon cross-section, CoPE/CrPE and diffraction-efficiency plots) and the digitized CSV / `.mat` data the scripts consume.

## Running

Run each script from inside this folder in MATLAB so the relative paths to the CSV and `.mat` files resolve. No toolboxes beyond base MATLAB are required.
