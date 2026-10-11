# Desktop-Speaker

**Empirical study of how port length and diameter shape the bass response of a 3D-printed bass-reflex speaker.**

![MATLAB](https://img.shields.io/badge/analysis-MATLAB-0076A8?style=flat-square)
![CAD](https://img.shields.io/badge/CAD-Fusion-555555?style=flat-square)
![Course](https://img.shields.io/badge/MIT-2.671-8A1538?style=flat-square)

<p align="center">
  <img src="images/Isometric%20View.png" alt="Isometric view of the speaker enclosure" height="46%">
  <img src="images/Section%20View.png" alt="Cross-section view of the speaker enclosure" height="46%">
</p>

Carter Tucker, MIT 2.671 Measurement and Instrumentation (Spring 2026)
[Paper (PDF)](report/Tucker_C_2671_Final_Paper_2col.pdf) | [Poster (PDF)](report/Th_07_Tucker_C_PortedSpeakerPoster.pdf)

---

## Overview

Ported (bass-reflex) enclosures extend low-frequency output by adding a resonant air mass, a Helmholtz resonator, to the box. In practice, port dimensions are usually picked from heuristics or lumped-parameter models that ignore real effects such as turbulence, port end geometry, non-uniform enclosure pressure, and losses.

This project measures what port geometry actually does. A 3D-printed desktop speaker with interchangeable ports was driven with a 20 Hz to 20 kHz log sweep and recorded at 1 m. Two metrics were extracted over 126 to 500 Hz:

- **Average SPL**, the overall bass output.
- **RMS deviation**, the flatness of the response (lower is flatter).

These were measured across four port diameters plus a sealed baseline, and five port lengths, for **21 configurations** with **5 trials each**.

## Headline Findings

| Parameter | Metric | Slope | 95% CI | Significant |
|---|---|---|---|---|
| Diameter (17.5 to 70 mm) | Average volume | +0.042 dB/mm | [0.033, 0.052] | Yes |
| Diameter (17.5 to 70 mm) | RMS deviation | -0.012 dB/mm | [-0.015, -0.009] | Yes |
| Length (0 to 180 mm) | Average volume | -0.0053 dB/mm | [-0.0091, -0.0014] | Yes |
| Length (0 to 180 mm) | RMS deviation | +0.0012 dB/mm | [-0.0002, 0.0023] | No |

Sealed baseline: 23.25 ± 0.43 dB average level, 1.86 ± 0.03 dB RMS deviation.

- Going from the smallest to the largest port diameter raised average output by about **2.2 dB** and cut RMS deviation by about **0.6 dB**.
- Going from the shortest to the longest port lowered average output by about **0.9 dB**. Its effect on flatness was not statistically significant.
- Both significant trends agree in direction with lumped-parameter predictions (output rises with port area and falls with port length). However, the predicted Helmholtz resonances (about 31 to 641 Hz across the tested ports) fall largely outside the 126 to 500 Hz analysis band. This limits how directly the theory can be tested here; see [Limitations](#limitations).

> **TODO:** Drop in 2 to 4 result plots here (volume vs. diameter, RMS deviation vs. diameter, and the same for length). Export them from the poster or regenerate them with the analysis code, save to `images/`, and embed with `![caption](images/your_plot.png)`.

## Experimental Setup

<p align="center">
  <img src="images/Modified%20Speaker%20Isometric.png" alt="Modified speaker isometric view" width="46%">
  <img src="images/Modified%20Speaker%20Section.png" alt="Modified speaker cross-section" width="46%">
</p>

> **TODO:** One-line caption explaining the "Modified" renders (what changed, e.g. the lid and port interface).

```
Computer / USB (sweep_sine.wav)
        |
   40 W amplifier
        |
 Passive crossover (~2 kHz)
     |             |
  Tweeter     Midrange driver  ->  Enclosure + interchangeable port
                                          |
                          Far-field mic @ 1 m  ->  LabQuest Mini  ->  LoggerPro  ->  CSV
```

| Item | Detail |
|---|---|
| Stimulus | 12 s logarithmic sine sweep, 20 Hz to 20 kHz ([`sweep_sine.wav`](sweep_sine.wav)) |
| Microphones | Vernier MCA (rated 100 Hz to 15 kHz), 1 m far-field used for analysis, 2 cm near-field recorded |
| DAQ | Vernier LabQuest Mini, 12-bit, 10 kHz sample rate, logged in LoggerPro |
| Amplifier | 40 W |
| Crossover | Passive, ~2 kHz (only the midrange driver matters below 500 Hz) |
| Repetitions | 5 trials per configuration, sweep played twice per trial and averaged |
| Ports | Diameters 17.5, 35.0, 52.5, 70.0 mm, plus sealed (0 mm); lengths 0, 45, 90, 135, 180 mm |
| Port seal | PTFE tape on the mating surface, friction-fit |
| Midrange driver / tweeter | **TODO** |
| Amplifier model | **TODO** |
| Enclosure internal volume | **TODO** |
| Enclosure outer dimensions | **TODO** |

### Fabrication

- **CAD:** Fusion. A STEP export is included: [`Optimized Speaker Box.step`](Optimized%20Speaker%20Box.step).
- **Printer:** Bambu Lab P1S.
- **Enclosure:** 4 walls, high infill (~50%).
- **Ports (resonators):** 100% infill.

> **TODO:** Confirm the infill figures above, add filament type, and add the STL files (the `.stl` files are not in the current folder tree).

### Test matrix

The sealed case (0 mm diameter) has no physically meaningful length, giving 21 geometries in total.

| Length \ Diameter | 0 | 17.5 | 35.0 | 52.5 | 70.0 |
|---|:-:|:-:|:-:|:-:|:-:|
| 0 | x | x | x | x | x |
| 45 | | x | x | x | x |
| 90 | | x | x | x | x |
| 135 | | x | x | x | x |
| 180 | | x | x | x | x |

## Method

Recordings are transformed with an FFT, divided by the FFT of the input sweep to form a transfer function, and reduced to two metrics over 126 to 500 Hz:

$$\bar{L} = \frac{1}{f_f - f_i}\int_{f_i}^{f_f} \mathrm{SPL}(f)\,df$$

$$\mathrm{RMS\ Deviation} = \sqrt{\int_{f_i}^{f_f} \left(\mathrm{SPL}(f) - \bar{L}\right)^2 df}$$

The Helmholtz resonance used for the model comparison is

$$f_p = \frac{c}{2\pi}\sqrt{\frac{S}{VL}}$$

where $S$ is the port cross-section, $L$ the port length, $V$ the enclosure volume, and $c$ the speed of sound.

```mermaid
flowchart LR
    A[LoggerPro CSV] --> B[FFT of recording and sweep]
    B --> C[Transfer function]
    C --> D[Average over 5 trials]
    D --> E[Average volume and RMS deviation, 126-500 Hz]
    E --> F[Linear fits with 95% CI]
```

## Repository Structure

```
Desktop-Speaker/
├── analysis_code/            MATLAB analysis (iterative versions, see below)
├── experimental_data_csv/    Raw LoggerPro exports
├── experimental_data_cmbl/   LoggerPro project files
├── images/                   CAD renders
├── report/
│   ├── Tucker_C_2671_Final_Paper_2col.pdf
│   └── Th_07_Tucker_C_PortedSpeakerPoster.pdf
├── Optimized Speaker Box.step
├── sweep_sine.wav
└── README.md
```

> **TODO:** Fill in the file naming convention for the CSV and `.cmbl` files (for example, how port length, diameter, and trial number are encoded).

## Reproducing the Analysis

**Requirements:** MATLAB (**TODO:** version and toolboxes, likely Signal Processing and Curve Fitting for `cftool`).

1. Clone the repo.
2. Open MATLAB and set the working directory to `analysis_code/`.
3. **TODO:** Name the main script (probably the highest-numbered `analyze_ports_*.m`) and state whether paths to `experimental_data_csv/` need editing.
4. **TODO:** State the run order and outputs, e.g. which script produces the per-configuration metrics and which one, [`averaged_metrics_cftool.m`](analysis_code/averaged_metrics_cftool.m), produces the linear fits and plots.

The `analyze_ports_*.m` files are successive iterations of the same analysis, kept for history. **TODO:** Mark which is final, or move older versions into `analysis_code/archive/`.

## Limitations

- **Resonance placement.** Predicted Helmholtz resonances (about 31 to 641 Hz) mostly do not sit inside the 126 to 500 Hz window, so port behavior depends heavily on where each resonance falls relative to the band. This limits what the experiment can say about the quantitative accuracy of the theory.
- **Microphone range.** The Vernier MCA is rated from about 100 Hz, close to the 126 Hz lower analysis bound.
- **Idealized model.** The lumped-parameter model assumes lossless acoustics, uniform enclosure pressure, and ideal cylindrical ports. Real ports have an effective length that differs from the physical length (end effects, flanges, clearance), and losses lower the effective Q.
- **Materials and rig.** Printed-part permeability and amplifier variability add uncontrolled error. Tests were done in an ordinary room, not an anechoic chamber.

**Future work:** design enclosures whose resonances sit consistently below the test band, use less permeable materials, and minimize edge effects, so the measurements reflect port-enclosure interaction rather than rig artifacts.

## Why I Built This

> **TODO:** A short paragraph in your own voice: what drew you to speakers and acoustics, and why this question.

## What I Learned

> **TODO:** Two or three sentences on what surprised you (for example, how the sealed baseline compared to the ported cases) and what you would change next time.

## Citation

```bibtex
@misc{tucker2026ported,
  author = {Tucker, Carter},
  title  = {The Effect of Port Geometry on the Frequency Response of Desktop Speakers},
  year   = {2026},
  note   = {MIT 2.671 Measurement and Instrumentation},
  url    = {https://github.com/TODO-USERNAME/Desktop-Speaker}
}
```

## Acknowledgments

Thanks to Prof. Sili Deng, Jeanne Wildman, and Kevin DiGenova for guidance and feedback, and to Dr. Kevin Cedrone for flexibility and support during equipment issues.

## License

> **TODO:** Choose a license (e.g. MIT for code, CC BY 4.0 for the report and CAD) and add a `LICENSE` file.

## Contact

Carter Tucker: [LinkedIn](https://www.linkedin.com/in/TODO) | **TODO:** email or personal site
