![Lanuguage](https://img.shields.io/badge/-Fortran-734f96?logo=fortran)
![Python](https://img.shields.io/badge/-Python3-3776AB?style=flat&logo=python&logoColor=white)
![Version](https://img.shields.io/static/v1?label=version&message=6.6&color=blue)
![Downloads](https://img.shields.io/github/downloads/YoushanLiu/ANCC/total)
![Commits](https://img.shields.io/github/commit-activity/m/YoushanLiu/ANCC)
![Size](https://img.shields.io/github/languages/code-size/YoushanLiu/ANCC)
![RepoSize](https://img.shields.io/github/repo-size/YoushanLiu/ANCC)
![License](https://img.shields.io/github/license/YoushanLiu/ANCC)
![Stars](https://img.shields.io/github/stars/YoushanLiu/ANCC)
![Forks](https://img.shields.io/github/forks/YoushanLiu/ANCC)
# ANCC
Some scripts and programs for ambient noise tomography, such as computing cross-correlation, auto-correlation (with parallel computing), and extracting dispersion curves.
It can also output prestack cross-correlation functions.

1. Convert the REFTEK file format (or other file formats) to the SAC file format (scripts)

    reftek2sac.py     -> file convertor with the folder structure, such as /root/stage_and_station_name/day/
    
    reftek2sac_p.py   -> parallel processing of reftek2sac.py
    
    reftek2sac2.py    -> file convertor with the folder structure, such as /root/stage/station_name/day/
    
    reftek2sac2_p.py  -> a parallel processing of reftek2sac2.py

	Actually, it is easily modified to convert other file formats by using the flexibility of ObsPy's read function, i.e. change the FORMAT parameter in the read of the ObsPy.
	A file format should be given, it will improve the performance.

2. Cut daily data into predefined-segment data (scripts)

   This program first merges daily data then cut it into predefined segments. Simultaneously, it resamples data at each natural sampling point.

    cutdata.py        -> cut daily data into predefined-segment data with the folder structure, such as /root/stage_and_station_name/day/
    
    cutdata_p.py      -> parallel processing of cutdata.py
    
    cutdata2.py       -> cut daily data into predefined-segment data with the folder structure, such as /root/stage/station_name/day/
    
    cutdata2_p.py     -> parallel processing of cutdata2.py

3. Auto- or cross-correlation and frequency-time analysis (Ambient Noise Auto- or Cross-Correlation)

    It includes the AFTAN to extract the dispersion curve.

4. Generate instrument response automatically (makePZfiles)

    In the makePZfiles folder, you can generate instrument response polezero files automatically using resp2pz.py or makePZs.py, 
    once fill an Excel sheet based on your project log. The former can take responses of sensor and DAS into consideration, 
    it then merges them into one PZ file. While the latter only takes the response of the sensor into consideration, then genertes
    corresponding PZ files automatically.


# Prerequisites
    The MPICH/OpenMPI is required. There are many MPI release versions. Usually, MPICH is recommended.
	The number of processes that can be used by the former (MPICH) is twice the number of physical 
    cores, while the number of processes that can be used by the latter (OpenMPI) is equal to the 
	number of physical cores.
    The FFTW3, SAC, and MPI (MPICH2/OpenMPI) softwares are required by the ANCC.
    The ObsPy is required by the scripts.
	The xlrd is required by the makePZs.py and resp2pz.py.



# LICENSE
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
