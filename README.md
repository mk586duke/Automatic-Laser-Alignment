# Automatic-Laser-Alignment
## Introduction
This project implements an automated laser alignment system designed to optimize optical coupling efficiency for experiments at the Duke Quantum Center. Traditionally, aligning laser beams into fibers or detectors is a labor-intensive manual process that can take hours and suffer from inconsistency.

Our solution automates this 4-axis alignment task using a Nelder-Mead optimization algorithm running on an Arduino Mega. By controlling motorized mirror mounts, the system achieves high-precision alignment in minutes with superior repeatability compared to manual methods.

## System Specifications
### Targeting accuracy:
-  Lateral error ≤ ± 250 μm
-  Angular error ≤ ± 5 mrad
### Performance goals:
-  System improves output power relative to misaligned start
-  Manual override available at any time
-  Typical alignment cycle completes within a few to tens of minutes
- System is easy to install onto the prexisting hardware
## Design Description
### Software Architecture

The software is divided into two main components: the embedded firmware running on the Arduino Mega and a Python-based simulation environment.

#### Firmware (Arduino)
The core logic resides in `arduino_full.ino`, which implements the control loop and optimization algorithm.
- **Optimization Algorithm**: A custom C++ implementation of the **Nelder-Mead simplex method** (`nelderMead4D`). It optimizes a 4-dimensional space (X/Y axes for 2 mirrors) to maximize light intensity.
- **Motor Control**: The `setMirrorPositions` function manages four `Stepper` objects (using the standard Arduino `Stepper.h` library). It handles coordinate tracking and sequential motor movement.
- **Sensor Interface**: The `measureIntensity` function reads the analog voltage from the photodetector (pin A0), averaging multiple samples (`ADC_SAMPLES`) to reduce noise.
- **Simulation Mode**: A `SIMULATED` flag allows testing the optimization logic without physical hardware by using a mathematical model of the laser intensity profile.

#### Python Simulation
The `laser_sim.py` script provides a visual simulation environment using **Panda3D**.
- **Visualization**: Renders the laser beam, mirrors, and target in 3D space.
- **Algorithm Verification**: Uses `scipy.optimize.minimize` (Nelder-Mead) to validate the optimization strategy before deploying to hardware.
- **Physics Model**: Simulates laser reflection and intensity falloff based on mirror angles.

### Hardware Architecture
#### Pre-Existing Hardware Setup
The standard hardware setup for the automatic laser alignment system includes the following components.
- One photodetector (e.g., Thorlabs PDA100A) to measure laser intensity TODO: include an image
- Two 2-axis mirror mounts (K100)  TODO: include an image
- One 760nm laser source TODO: include an image & asking Thomas for more details

#### New Components
The hardware setup consists of the following components:
- Four 2048-step steppers (e.g., 28BYJ-48) for two 2-axis mirrors (controls) TODO: include an image
- Four motor drivers (TODO: list the specific models)
- One photodetector (e.g., Thorlabs PDA100A) to measure laser intensity TODO: include an image
- Microcontroller (e.g., Arduino Uno) to control the stepper motors and read photodetector values TODO: include an image

#### Power and Networking
- The photodetector transmits to the Arduino Mega via analog input pin A0, outputting an analog voltage proportional to the detected light intensity.
- The stepper motors communicate through the 
- The Arguino Mega coordinates the motor drivers
- All power comes from

## Testing and Validation

### Theoretical Resolution
Based on the hardware specifications defined in the firmware (`arduino_full.ino`), we can calculate the theoretical mechanical resolution of the system.

*   **Stepper Motor**: 2048 steps per revolution (28BYJ-48).
*   **Mirror Mount Gearing**: 2 motor revolutions per 1 degree of mirror tilt (`REV_PER_DEGREE`).
*   **System Steps per Degree**: $2048 \times 2 = 4096$ steps/degree.

**Angular Resolution**:
$$ \theta_{res} = \frac{1}{4096} \approx 0.000244^\circ \approx 4.26 \mu\text{rad} $$

**Lateral Resolution**:
At a target distance of 50 cm, the theoretical beam displacement per step is:
$$ \Delta x = d \times \tan(\theta_{res}) \approx 500 \text{ mm} \times 4.26 \times 10^{-6} \approx 2.13 \mu\text{m} $$

This resolution is sufficient for coupling into single-mode fibers (typical core diameter ~9 $\mu$m) and far exceeds the precision of manual adjustment knobs.

### Algorithm Validation
The Nelder-Mead optimization logic was validated using the firmware's `SIMULATED` mode, which models a Gaussian intensity profile.
*   **Convergence**: The algorithm consistently converges to the global maximum (>99.9% intensity) within 60-70 iterations.
*   **Accuracy**: In simulation, the final position error is typically < 1 step (< 0.0002°).
*   **Robustness**: The system successfully recovers optimal alignment from arbitrary starting positions within the capture range.

## Demonstration Videos
