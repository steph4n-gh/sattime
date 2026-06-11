# Sattime Antenna Design and Tuning Guide

This document provides the mathematical models, engineering specifications, and construction guidelines for fabricating and tuning Software Defined Radio (SDR) antennas optimized for the Low Earth Orbit (LEO) satellite swarms tracked by the `sattime` receiver.

---

## 1. Principles of LEO Satellite RF Reception

Receiving high-speed digital telemetry or carrier leakages from satellites orbiting between $400\text{ km}$ and $1500\text{ km}$ altitude requires consideration of the following physical factors:
- **Polarization Mismatch**: Most LEO satellites transmit using Right-Hand Circular Polarization (RHCP) to prevent fading as the satellite rotates and tumbles relative to the ground station. Using a linearly polarized receiver antenna (like a simple whip or dipole) introduces a constant $3\text{ dB}$ polarization mismatch loss and periodic signal nulls.
- **Free-Space Path Loss (FSPL)**: FSPL increases with the square of the frequency and distance:
  $$\text{FSPL} = \left(\frac{4\pi d f}{c}\right)^2$$
  At VHF ($137\text{–}150\text{ MHz}$), path loss is significantly lower than at L-band ($1626\text{ MHz}$), allowing weak carrier locks even from inside buildings using sub-optimal antennas.
- **Doppler Shift**: The satellite’s relative velocity shifts the received frequency:
  $$f_{\text{received}} = f_{\text{nominal}} \left(1 - \frac{v_{\text{range\_rate}}}{c}\right)$$
  The antenna bandwidth must be wide enough to accommodate the full Doppler sweep (up to $\pm 3.7\text{ kHz}$ at VHF, and $\pm 40\text{ kHz}$ at L-band).

---

## 2. Satellite Swarms and RF Profiles

| Satellite Swarm | Center Downlink Frequency ($f$) | Polarization | Wavelength ($\lambda$) | Recommended Antenna Design |
| :--- | :--- | :--- | :--- | :--- |
| **NOAA Weather** | $137.50\text{ MHz}$ | RHCP | $2.180\text{ m}$ | Quadrifilar Helix (QFH) or V-Dipole |
| **Orbcomm M2M** | $137.50\text{ MHz}$ | RHCP | $2.180\text{ m}$ | Quadrifilar Helix (QFH) or Turnstile |
| **Amateur Satellites** | $145.90\text{ MHz}$ | RHCP/Linear | $2.055\text{ m}$ | Eggbeater or Turnstile |
| **Starlink VHF** | $150.80\text{ MHz}$ | Vertical Linear | $1.988\text{ m}$ | Quarter-Wave Whip on Ground Plane |
| **Iridium L-Band** | $1626.27\text{ MHz}$ | RHCP | $0.184\text{ m}$ | Active RHCP Patch Antenna |

---

## 3. VHF Antenna Tuning Equations

To calculate the physical length of half-wave ($\lambda/2$) dipoles or quarter-wave ($\lambda/4$) monopole elements, we must account for the **Velocity Factor** ($V_f$) of the metal conductor. The speed of electromagnetic propagation in copper or aluminum is slower than in a vacuum:
$$v = V_f \cdot c$$
For typical antenna conductors (e.g., $2.0\text{–}4.0\text{ mm}$ copper wire or aluminum rods), the velocity factor is $V_f \approx 0.95\text{–}0.97$.

### 3.1 Quarter-Wave Element Length Formula
$$L_{\lambda/4} = \frac{c \cdot V_f}{4 \cdot f} = \frac{299.792458 \cdot 0.95}{4 \cdot f_{\text{MHz}}} \approx \frac{71.2}{f_{\text{MHz}}}\text{ meters}$$

### 3.2 Half-Wave Element Length Formula
$$L_{\lambda/2} = \frac{c \cdot V_f}{2 \cdot f} = \frac{299.792458 \cdot 0.95}{2 \cdot f_{\text{MHz}}} \approx \frac{142.4}{f_{\text{MHz}}}\text{ meters}$$

---

## 4. Swarm-Specific Antenna Designs

### 4.1 Starlink VHF ($150.80\text{ MHz}$)
Starlink telemetry leakages are vertically polarized. The most effective omnidirectional antenna is a quarter-wave monopole (whip) mounted vertically over a conducting ground plane.

#### Dimensions at $150.80\text{ MHz}$:
- **Vertical Whip Length**:
  $$L = \frac{299.792458 \cdot 0.97}{4 \cdot 150.80} = 0.482\text{ m} \approx 48.2\text{ cm}$$
- **Radial Elements (Ground Plane)**: Cut 4 radials at $\approx 51.0\text{ cm}$ (angled downwards at $45^\circ$ to match the feedpoint impedance to $50\,\Omega$).

#### Basement/Indoor Test Setup:
For software and EKF verification, a simple magnetic-mount whip cut to $49.7\text{ cm}$ (tuned to account for base matching coils) placed on a large metal baking sheet (cookie sheet) in a basement can lock onto overhead Starlink passes. The metal sheet acts as an image plane, forming the missing half of the dipole.

---

### 4.2 NOAA Weather & Orbcomm M2M ($137.50\text{ MHz}$)
These swarms broadcast circular polarization. Linear antennas will suffer periodic deep fades (polarization rotation nulls) as the satellite passes.

#### 1. The V-Dipole (Linear Approximation)
For a low-cost entry, build a V-dipole using two copper wire legs:
- **Leg Length**: $52.0\text{ cm}$ per dipole element.
- **Apex Angle**: Mount elements at a $120^\circ$ angle in a horizontal plane.
- **Orientation**: Point the apex of the "V" North, with the legs extending South-East and South-West.

#### 2. The Quadrifilar Helix (QFH)
The QFH consists of two orthogonal helical loops fed $90^\circ$ out of phase. It provides an RHCP hemispherical pattern with high gain at the zenith and near the horizon.
- **Tall Loop Height**: $L_{\text{tall\_height}} \approx 81.3\text{ cm}$, Width $\approx 35.6\text{ cm}$
- **Short Loop Height**: $L_{\text{short\_height}} \approx 77.2\text{ cm}$, Width $\approx 33.8\text{ cm}$
- **Phasing**: The loop lengths are sized differently so that one loop is inductive and the other is capacitive at $137.5\text{ MHz}$, creating the required $90^\circ$ phase shift when fed in parallel.

---

### 4.3 Amateur Satellites ($145.90\text{ MHz}$)

#### The Turnstile Antenna (Crossed Dipoles)
A turnstile antenna consists of two identical half-wave dipoles mounted orthogonally ($90^\circ$ mechanical displacement) and fed with a $90^\circ$ electrical phase delay.

#### Dimensions at $145.90\text{ MHz}$:
- **Dipole Leg Length**: $49.0\text{ cm}$ ($98.0\text{ cm}$ tip-to-tip per dipole).
- **Phasing Harness**: Connect the two dipoles using a quarter-wavelength delay line of $75\,\Omega$ coaxial cable (RG-59 or RG-11).
  The physical length of the $75\,\Omega$ delay line is:
  $$L_{\text{delay}} = \frac{c \cdot V_{\text{coax}}}{4 \cdot f} = \frac{299.792458 \cdot 0.66}{4 \cdot 145.90} \approx 33.9\text{ cm}$$
  *(where $V_{\text{coax}} = 0.66$ for solid polyethylene dielectric coax like RG-59).*

---

### 4.4 Iridium L-Band ($1626.27\text{ MHz}$)
At L-band, atmospheric attenuation and building penetration losses are severe. Passive whip antennas or indoor setups are insufficient.

#### Requirements:
- **Antenna Type**: Active Ceramic Patch Antenna or Helical Patch.
- **Built-in LNA**: Must feature a low-noise amplifier ($20\text{–}30\text{ dB}$ gain) integrated directly at the feedpoint.
- **Bias-Tee Powering**: The SDR or daemon must supply $+3.3\text{ V}$ or $+5.0\text{ V}$ DC up the coaxial cable to power the LNA.
- **Element Size**: At $1626\text{ MHz}$, $\lambda \approx 18.4\text{ cm}$. A patch element is approximately $4.6\text{ cm}$ square.
- **Placement**: Outdoor installation with an unobstructed line-of-sight to the horizon is mandatory.

---

## 5. Impedance Matching and Baluns

- **Common-Mode Currents**: When feeding a balanced antenna (like a dipole) with an unbalanced transmission line (coaxial cable), RF current will flow down the outer shield of the coax. This distorts the antenna pattern and introduces local electromagnetic interference.
- **Choke Balun (Ugly Balun)**: Wind 5 to 8 turns of the feedline coaxial cable into a tight solenoid (diameter $\approx 10\text{ cm}$) directly at the antenna feedpoint. This creates a high common-mode impedance that blocks currents from flowing down the shield.
- **Impedance Mismatch**: While SDRs are designed for $50\,\Omega$ systems, using standard $75\,\Omega$ television coax (RG-6) introduces a Voltage Standing Wave Ratio (VSWR) of $1.5:1$, representing a minor reflection loss of only $4\%$. This loss is negligible compared to the benefits of using low-loss double-shielded RG-6 cable.

---

## 6. Tuning with a Vector Network Analyzer (VNA)

To verify the resonance frequency and impedance matching of your constructed antenna:
1. **Calibrate the VNA**: Perform a standard Short-Open-Load-Thru (SOLT) calibration on the VNA up to the coaxial connector interface.
2. **Measure Return Loss ($S_{11}$)**: Connect the antenna and locate the frequency of the deepest dip in the $S_{11}$ trace.
   - For optimal performance, $S_{11}$ should be below $-15\text{ dB}$ at the center frequency (corresponding to a VSWR of $< 1.43:1$).
3. **Trimming Elements**:
   - If the dip is at a lower frequency than targeted, the antenna elements are too long. Trim the conductors in small increments ($2\text{–}5\text{ mm}$).
   - If the dip is at a higher frequency, the elements are too short and must be replaced.
