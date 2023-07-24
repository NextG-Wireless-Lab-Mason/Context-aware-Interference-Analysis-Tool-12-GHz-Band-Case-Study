# Context-aware Interference Analysis Tool - Case Study of 12 GHz Band
**#ProJect Overview:**
  This project proposes a holistic, multi-disciplinary, context-aware spectrum-sharing approach to address the spectrum coexistence of broadband wireless systems in satellite bands. In this work, we developed a weather-dependent (sunny & rainy) simulation-based interference evaluation framework for the 12 GHz band considering a realistic deployment of the FSS receiver (1770 Forecast Drive, Blacksburg, Virginia), the exact position of all 33 MBS (Macro Base Station)from opencellID within a 5000m circular area, and Geolocation data of 8644 Buildings from OpenStreetMap.The DSA framework will transmit exclusion zone and weather details, base station information, and interference thresholds for both sunny and rainy weather to the simulator. Therefore, the simulator will start analyzing the interference for each MBSs based on the industry-standardized beamforming at MBSs, directional signal reception at FSS receivers, channel scheduling at 5G MBSs, Set of Interfering Beams from 5G MBS, Path Loss between Interfering MBS and FSS, and Noise power. Moreover, the aggregate interference-to-noise ratio is calculated for each MBS. If the interference-to-noise ratio (I/N) for an MBS  meets the threshold requirements considering the weather and exclusion zones, the MBS will be marked as "Active"; otherwise, it will be marked as "Inactive."  

**#How to Run the code:**
  1. Load all the data from the folder data or 
  2. Run the res_server.py from the DSA folder that will send the necessary information to the simulator.py
  3. Run the Simulator.py 



**#Acknowledgement:**
  The National Science Foundation under Award # 2128584 generously supports this project. Link: https://www.nsf.gov/awardsearch/showAward?AWD_ID=2128584 

**Relevant Publications**
  Niloy T. R., Hassan Z., Stephenson N., and Shah, V. K., "Interference Analysis of Coexisting 5G Networks and NGSO FSS Receivers in the 12 GHz Band," in IEEE Wireless Communications Letters, doi: 10.1109/LWC.2023.3281769. Link: https://ieeexplore.ieee.org/document/10139318 
