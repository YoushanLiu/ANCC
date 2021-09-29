# ANT_Preprocessing
Some scripts for ambient noise tomography (support parallel computing)

1. convert reftek file format (or other file formats) to SAC file format

  reftek2sac.py    -> file convertor with the folder structure, such as /root/period_and_station_name/day/
  reftek2sac_p.py  -> a parallel processing of reftek2sac.py
  reftek2sac2.py   -> file converter with the folder structure, such as /root/period/station_name/day/
  reftek2sac2_p.py -> a parallel processing of reftek2sac2.py

2. cut daily data into predefined-segment data

  This program first merges daily data, then cut it into predefined-segment. Simultaneously, it resamples data at each sampling point.

  cutdata.py       -> cut daily data into predefined-segment data with the folder structure, such as /root/period_and_station_name/day/
  cutdata_p.py     -> a parallel processing of cutdata.py
  cutdata2.py      -> cut daily data into predefined-segment data with the folder structure, such as /root/period/station_name/day/
  cutdata2_p.py    -> a parallel processing of cutdata2.py

3. cross-correlation and frequency-time analysis

  It utilizes the AFTAN.
