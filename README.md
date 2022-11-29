# The development and validation of an Inhomogeneous Wind Scheme for Urban Street
\
\
IWSUS-v0.1: A model for wind profile calculation at typical positions of a long street.
Inhomogeneous Wind Scheme for Urban Street (IWSUS) is developed for genrating wind profiles at typical position of a long street. 
\
\
The files includes the equations of these wind profiles and am example for a daily variation of heat flux with the use of daily urban surface temperature data in data folder. They has been coupled in WRF/WRF-chem (Verson 4.2.1)by being drived in the land surface schemes. 
\
\
To apply the IWSUS in WRF, please unfold the package "Registry"、"phys" and "dyn_em" and replace the files with the same in original folders, then compile again the WRFV3 module. 
\
\
The flag "urban_ahs_on=1" should be added in the file namelist.input to call IWSUS scheme. 
\
\
The Gothernburg case is used to validate surface flux from a real street measurement in Gothenburg (Offerle, B., Eliasson, I., Grimmond, C. S. B., and Holmer, B.: Surface heating in relation to air temperature, wind and turbulence in an urban street canyon, Boundary-Layer Meteorology, 122, 273-292.https://doi.org/10.1007/s10546-006-9099-8, 2006.).
