REM Upgrade conda if needed.
conda upgrade -y conda 

REM Uninstall any old verision of omin.
pip uninstall -y omin

REM These dependencies need to be installed via conda.
conda install -cy conda-forge ipywidgets
conda install -y bokeh
conda install -cy conda-forge jupyter_contrib_nbextensions
conda install -cy conda-forge jupyter_nbextensions_configurator

REM Turn on hinterland.
jupyter nbextension enable hinterland/hinterland

REM Make sure ipywidgets is turned on.
jupyter nbextension enable --py --sys-prefix widgetsnbextension

REM Create jupyter notebook config file.
jupyter notebook -y --generate-config

REM Install omin.
python setup.py install

REM omin install files left behind.
rmdir build /s /q
rmdir dist /s /q
rmdir omin.egg-info /s /q

REM Create a standard JupyterNotebooks dir.
cd %userprofile%\Documents
mkdir JupyterNotebooks

REM Add lines to .jupyter/jupyter_notebook_config.py file so jupyter notebook always opens to the new default dir.
cd %userprofile%\.jupyter
(echo. && echo # Make jupyter notebook open in ~\Documents\JupyterNotebook && echo import os) > scratch.txt

type scratch.txt >> jupyter_notebook_config.py

echo c.NotebookApp.notebook_dir = os.path.expanduser("~/Documents/JupyterNotebooks") > scratch.txt

type scratch.txt >> jupyter_notebook_config.py

del scratch.txt

@echo I'M FINISHED!

timeout 5