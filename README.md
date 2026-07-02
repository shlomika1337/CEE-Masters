# CEE Masters

Graduate technical work from my M.Sc. in Computer and Electrical Engineering at Ben-Gurion University of the Negev. My specialization is RF integrated circuits, nanotechnology, and VLSI; these works sit on the semiconductor-device, photonics, and applied-electromagnetics side of the program.

These are not summaries. Each builds the physics from first principles, and the metasurface work in particular reproduces and critiques a published result with independent simulation. The transfer-matrix and Huygens-Fresnel optical modeling, semiconductor transport, and thin-film physics here overlap directly with integrated-circuit and RF work.

## Reports

### Metasurface Optics
*Folder: `Metasurface Optics/`*

An independent reproduction and critique of Lin, Fan, Hasman, and Brongersma, "Dielectric gradient metasurface optical elements," *Science* 345(6194), 2014, a landmark flat-optics paper on silicon-nanobeam metasurfaces that use the geometric Pancharatnam-Berry phase for wavefront control.

Rather than summarizing, this paper verifies and challenges the original:

- Reproduced the axicon, metalens, and beam-deflector results with my own MATLAB Huygens-Fresnel propagation simulations, recovering the 6:1 Bessel-beam signature, a near-diffraction-limited metalens spot (~0.6 um FWHM, ~0.05 um focal shift), and the deflection behavior.
- Found that the nanobeam half-wave-plate is non-ideal at the paper's 550 nm design wavelength: by interpolating the authors' own transmission and phase data and computing the co- and cross-polarization efficiencies (CoPE/CrPE), I showed the true polarization-conversion optimum sits closer to 520 nm, and quantified that up to ~86% of incident light is lost to reflection and absorption.
- Identified an improved meta-atom design point (136 nm width, 230 nm period) with sub-1% phase error and near-unity transmission ratio, and then honestly characterized the CoPE/CrPE trade-off that prevents it from being a free improvement.
- Corrected an inaccuracy in the paper's TE/TM polarization terminology.
- Worked around real hardware limits: the far-field data ran to 17000x17000 Float64 matrices (~50 GB), handled by downsampling, with the limits of what my hardware could integrate stated plainly.

The MATLAB reproduction and analysis code is included under `Metasurface Optics/MATLAB/` (see its own README) (Huygens-Fresnel propagation for the axicon and metalens, and the CoPE/CrPE and design-sweep analysis). The full-wave Ansys Lumerical simulations that produce the meta-atom sweep data were built and run in collaboration with a colleague, as credited in the paper; the MATLAB propagation reproductions, the CoPE/CrPE analysis, and the critique are my own.

### Organic Semiconductors Applications
*Folder: `Organic Semiconductors Applications/`*

A 26-page review of the device physics and engineering of series-connected tandem organic photovoltaic cells (OSCs), argued from first principles.

It builds the case for tandems from the Shockley-Queisser detailed-balance limit and its four intrinsic loss channels, extends it through de Vos's multi-junction treatment (and why the one-to-two-junction step dominates the gain), and then develops the organic-specific device physics: why the low dielectric constant of organic semiconductors produces tightly bound Frenkel excitons, why this forces the donor-acceptor bulk-heterojunction architecture, and how Marcus electron-transfer theory and charge-transfer-state energetics set the voltage. It covers series vs. parallel interconnection with the current- and voltage-matching conditions, the strict electrical, optical, and chemical requirements on the interconnecting layer (ICL), and the optical design problem, treating the sub-wavelength stack as a coherent cavity and computing field-driven generation via the transfer-matrix method to current-match the subcells. It closes on the state of the art (single-junction OSCs past 20.82%, tandems past a certified 21.2%) and the efficiency-stability conundrum that now dominates the field, citing 2025-2026 literature throughout.

### Nanoelectronics
*Folder: `Nanoelectronics/`. Bandgap-Engineered Solar Cells, in Hebrew.*

A 39-page survey, in Hebrew, of solar-cell physics and bandgap engineering, built from first principles and supervised by Dr. Ilan Shalish. It opens with a terminology glossary and the global energy motivation, then develops the physics from the photoelectric effect and band structure through the ideal- and double-diode solar-cell models, the Shockley-Queisser limit, and the Varshni temperature dependence of the bandgap. It then covers advanced structures and bandgap engineering in depth: heterojunctions and the Anderson model (with its correction), alloys and solid solutions, quantum wells (via the confined Schrodinger problem), quantum dots and the blue shift (Brus equation), Bragg reflectors, bandgap grading, perovskite composition engineering, and strain engineering through the deformation potential.

It closes with two original numerical problems, each with a full worked solution. The first derives an approximate fill-factor expression directly from the ideal-diode model and uses it to compute the fill factor and PCE of perovskite cells at two cesium concentrations, validating against published values (~1% error at the well-behaved concentration) and reasoning about why the model breaks down at the phase-segregating one. The second is an intermediate-band / quantum-dot solar-cell comparison that derives the open-circuit-voltage and efficiency ratio between a QD and non-QD GaAs cell from Kirchhoff's laws and the diode model, and deliberately embeds distractor data to reward critical reading. Writing solvable original problems, deriving the results, and checking them against the literature is the part that best demonstrates command of the material.

## Topics

Semiconductor device physics, nanophotonics and metasurfaces, geometric (Pancharatnam-Berry) phase, Huygens-Fresnel and transfer-matrix optical modeling, photovoltaics and bandgap engineering, and charge transport in disordered semiconductors.

## About

M.Sc. student in Computer and Electrical Engineering at Ben-Gurion University of the Negev, specializing in RF integrated circuits, nanotechnology, and VLSI. First place, Iron Codes 2024 CTF (Tel Aviv University CyberWeek).

My other repositories cover embedded hardware, RF, and cybersecurity work.
