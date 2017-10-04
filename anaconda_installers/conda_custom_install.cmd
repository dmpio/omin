:: Installer for Anaconda
:: James Draper 10/4/2017

REM This directory is subject to change.
set "conda_inst_dir=%userprofile%\Downloads"
set "conda_inst_fname=Anaconda3-5.0.0-Windows-x86_64.exe"
set "conda_inst_path=%conda_inst_dir%\%conda_inst_fname%"
REM Run Anaconda silient install for all users, adding variables to the path and registering system python.
REM MUST BE RUN AS ADMIN!
start /wait "" %conda_inst_path% /InstallationType=AllUsers /RegisterPython=1 /AddToPath=1 /S /D=C:\Anaconda3

@REM Upgrade conda if needed.
conda upgrade -y conda

@REM Uninstall any old verision of omin.
pip uninstall -y omin

@REM These dependencies need to be installed via conda.
conda install -c conda-forge ipywidgets -y
conda install -c conda-forge bokeh -y
conda install -c conda-forge jupyter_contrib_nbextensions -y
conda install -c conda-forge jupyter_nbextensions_configurator -y

@REM Turn on hinterland.
jupyter nbextension enable hinterland/hinterland

@REM Make sure ipywidgets is turned on.
jupyter nbextension enable --py --sys-prefix widgetsnbextension

@REM Create jupyter notebook config file.
jupyter notebook -y --generate-config

@REM Create a standard JupyterNotebooks dir.
cd %userprofile%\Documents
mkdir JupyterNotebooks

@REM Add lines to .jupyter/jupyter_notebook_config.py file so jupyter notebook always opens to the new default dir.
cd %userprofile%\.jupyter
(echo. && echo # Make jupyter notebook open in ~\Documents\JupyterNotebook && echo import os) > scratch.txt
type scratch.txt >> jupyter_notebook_config.py
echo c.NotebookApp.notebook_dir = os.path.expanduser("~/Documents/JupyterNotebooks") > scratch.txt
type scratch.txt >> jupyter_notebook_config.py
del scratch.txt

@echo I'M FINISHED!

timeout 5
