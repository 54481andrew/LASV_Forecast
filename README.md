# LASV_Forecast

The RunForecast.r script generates a forecast for Lassa virus on a
particular month-year. Run

Rscript RunForecast.r -d '2022-04'

to get the forecast for April, 2022. The script will automatically look for and download any new CHIRPS data that is required for the specified forecast date. The code should gracefully exit if there is not enough precipitation data for the requested forecast (e.g., if you tried '2023-04'). 