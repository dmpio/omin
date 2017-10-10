:: Upgrade conda if needed.
conda upgrade -y conda

:: These dependencies need to be installed via conda.
conda install -c conda-forge ipywidgets -y
conda install -c conda-forge bokeh -y
conda install -c conda-forge jupyter_contrib_nbextensions -y
conda install -c conda-forge jupyter_nbextensions_configurator -y

:: Turn on hinterland.
jupyter nbextension enable hinterland/hinterland

:: Make sure ipywidgets is turned on.
jupyter nbextension enable --py --sys-prefix widgetsnbextension

:: Create jupyter notebook config file.
jupyter notebook -y --generate-config

:: Create a standard JupyterNotebooks dir.
cd %userprofile%\Documents
mkdir JupyterNotebooks

:: Add lines to .jupyter/jupyter_notebook_config.py file so jupyter notebook always opens to the new default dir.
cd %userprofile%\.jupyter
(echo. && echo # Make jupyter notebook open in ~\Documents\JupyterNotebook && echo import os) > scratch.txt
type scratch.txt >> jupyter_notebook_config.py
echo c.NotebookApp.notebook_dir = os.path.expanduser("~/Documents/JupyterNotebooks") > scratch.txt
type scratch.txt >> jupyter_notebook_config.py
del scratch.txt

@echo I'M FINISHED!

timeout 5
