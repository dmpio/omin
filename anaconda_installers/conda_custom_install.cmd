:: Installer for Anaconda
:: James Draper 10/4/2017

REM This directory is subject to change.
set "conda_inst_dir=%userprofile%\Downloads"
set "conda_inst_fname=Anaconda3-5.0.0-Windows-x86_64.exe"
set "conda_inst_path=%conda_inst_dir%\%conda_inst_fname%"
REM Run Anaconda silient install for all users, adding variables to the path and registering system python.
REM MUST BE RUN AS ADMIN!
start /wait "" %conda_inst_path% /InstallationType=AllUsers /RegisterPython=1 /AddToPath=1 /S /D=C:\Anaconda3
