# Context-aware Interference Analysis Tool - Case Study of 12 GHz Band

**#Project Overview:**

This project proposes a holistic, multi-disciplinary, context-aware spectrum-sharing approach to address the spectrum coexistence of broadband wireless systems in satellite bands. In this work, we developed a weather-dependent (sunny & rainy) simulation-based interference evaluation framework for the 12 GHz band considering a realistic deployment of the FSS receiver (1770 Forecast Drive, Blacksburg, Virginia), the exact position of all 33 MBS (Macro Base Station) from opencellID within a 5000m circular area, and Geolocation data of 8644 Buildings from OpenStreetMap. The DSA framework will transmit exclusion zone and weather details, base station information, and interference thresholds for both sunny and rainy weather to the simulator. Therefore, the simulator will start analyzing the interference for each MBS based on the industry-standardized beamforming at MBSs, directional signal reception at FSS receivers, channel scheduling at 5G MBSs, Set of Interfering Beams from 5G MBS, Path Loss between Interfering MBS and FSS, and Noise power. Moreover, the aggregate interference-to-noise ratio is calculated for each MBS. If the interference-to-noise ratio (I/N) for an MBS  meets the threshold requirements considering the weather and exclusion zones, the MBS will be marked as "Active"; otherwise, it will be marked as "Inactive."  

**#How to Run the code:**

  1. Load all the data from the folder **data or split_csv.py**
  2. Run the **rest_server.py** from the DSA folder that will send the necessary information to the **simulator.py**
  3. Run the **simulator.py** which will analyze the I/N for all the MBSs.
  4. Run the **app.py**



**#Acknowledgement:**

This project is generously supported by the National Science Foundation under Award # 2128584. Find more information about the project here: https://www.nextgwirelesslab.org/current-projects/spectrum-sharing-in-satellite-bands

**Relevant Publications**

 [1] T. R. Niloy, Z. Hassan, N. Stephenson, and V. K. Shah, Interference Analysis of Coexisting 5G Networks and NGSO FSS Receivers in the 12 GHz Band, IEEE Wireless Comms. Letters (WCL), 2023 Link: https://ieeexplore.ieee.org/document/10139318 
 
 [2] T. Niloy, S. Kumar, A. Hore, Z. Hassan, E. Burger, C. Dietrich, J. Reed, and V. K. Shah, ASCENT: A Context-Aware Spectrum Coexistence Design and 
Implementation Toolset for Policymakers in Satellite Bands, IEEE Conference of Dynamic Spectrum Access Networks (DySPAN), 2024. Link: https://arxiv.org/abs/2402.05273
