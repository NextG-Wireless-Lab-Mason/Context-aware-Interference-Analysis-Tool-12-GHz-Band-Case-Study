-------------------------------------------------------------------------------------------------------

## Notes
1. Only radios of type GSM are being considered in the simulator.
2. UE information is randomly generated in the simulator by dividing region per BS into 3 sectors and then randomly 
placing UEs in them.
3. Simulator accepts "BS_Count" parameter from the DSA. This is used to randomly select BS_Count number of 
base stations from the CSV file. However, in actual implementation `(line 491)`, simulator just takes the first 
BS_Count lines from the CSV.
4. FSS co-ordinates are expected from the DSA framework.
5. I/N threshold:
   1. Rainy Weather: -12db
   2. Normal Weather: -8.5db
   3. Those violating this would need to be turned off. Hence, creating kind of an exclusion zone

-------------------------------------------------------------------------------------------------------

## Open Issues to Discuss
1. Where is the file `310.csv` and who is generating this file? 
2. Currently, simulator reads BS information from the csv file and then randomly generates the UEs.
   1. DSA framework currently does not have any idea about this distribution.
   2. Is DSA framework required to know about this distribution?
3. First we are reading entire CSV file and then selecting a few BS based on the `"base_station_count"` parameter.
   1. Should the DSA framework just read the `first "base_station_count" lines` from the CSV?
4. Should the FSS co-ordinates be hardcoded `for now`?
   1. What's the strategy to obtain this parameter in the future?
5. Simulator takes a parameter `radius` as input. What is this?
   1. Is it same as exclusion zone radius? 

-------------------------------------------------------------------------------------------------------