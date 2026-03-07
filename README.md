# Desktop Speaker System

## Overview
This project involved designing and building a compact desktop speaker system using a custom enclosure, off-the-shelf drivers, and a Bluetooth amplifier. The goal was to improve low-frequency performance of a mid-range driver using a port.

The enclosure geometry was optimized to flatten the bass response and reduce **port chuffing** at higher output levels.

---

## Design Goals

- Compact desktop form factor
- Rated power of approximately **40 W**
- Improved low-frequency response compared to sealed enclosures
- Minimal airflow noise through the bass port

---

## Engineering Work

- Designed a **ported enclosure** using 3D-printed structural components based on analysis in WinISD
- Selected and integrated a **tweeter and midrange woofer driver pair**
- Tuned the enclosure as a **Helmholtz resonator** to reinforce low-frequency response
- Optimized port geometry to reduce turbulence and chuffing
- Assembled amplifier, drivers, and enclosure into a fully functional speaker system

---

## Fabrication

- 3D printed enclosure components
- Mechanical assembly of drivers and enclosure
- Electrical integration of Bluetooth amplifier module

---

## Tools Used

- **CAD:** Fusion 360
- **Fabrication:** 3D printing
- **Analysis:** Acoustic modeling and Helmholtz resonance calculations

---

## Images

### Speaker Assembly

### Enclosure Design
![CAD]("Isometric View.png")

---

## Future Improvements

- Measurement of frequency response using two calibrated microphones
- Additional tuning of port length and enclosure volume according to current Helmholtz resonator theory, and experimental data
- Comparison between simulated and measured acoustic performance

## Future Experiments

- Modify the enclosure to have a modular port system.
- Run trials, changing port length and diameter, to improve the amplitude and flatness of the response over a 20 Hz - 20 kHz sine sweep and compare the results of the experiment with current theories about Helmholtz resonators.
- Current Plan: 4 diameter trials + control, 4 length trials + control, 3+ runs per trial for 120+ total trials.
