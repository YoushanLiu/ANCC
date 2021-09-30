# ANCC
Some scripts for ambient noise tomography, such as computing cross-correlation (with parallel computing), extracting dispersion curves.

1. convert the reftek file format (or other file formats) to the SAC file format (scripts)

    reftek2sac.py    -> file convertor with the folder structure, such as /root/period_and_station_name/day/
    
    reftek2sac_p.py  -> a parallel processing of reftek2sac.py
    
    reftek2sac2.py   -> file converter with the folder structure, such as /root/period/station_name/day/
    
    reftek2sac2_p.py -> a parallel processing of reftek2sac2.py

2. cut daily data into predefined-segment data (scripts)

   This program first merges daily data, then cut it into predefined-segment. Simultaneously, it resamples data at each natural sampling point.

    cutdata.py       -> cut daily data into predefined-segment data with the folder structure, such as /root/period_and_station_name/day/
    
    cutdata_p.py     -> a parallel processing of cutdata.py
    
    cutdata2.py      -> cut daily data into predefined-segment data with the folder structure, such as /root/period/station_name/day/
    
    cutdata2_p.py    -> a parallel processing of cutdata2.py

3. cross-correlation and frequency-time analysis (Ambient Noise Denoise)

    It utilizes the AFTAN.

4. generating instrument response automatically (makePZfiles)

    In makePZfiles folder, you can generate automatically instrument response polezero files using resp2pz.py or makePZs.py 
    once fill a excel sheet based on your project log. The former can take responses of sensor and DAS into consideration, 
    it will merge them into one PZ file. While, the latter only take the response of sensor into consideration, then genertes
    corresponding PZ files automatically.


# Prerequirements
    The SAC software is required by the AND.
    The obspy is required by the scripts.


ANCC is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ANCC is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
