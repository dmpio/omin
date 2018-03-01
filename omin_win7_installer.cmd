REM Copyright James Draper 2018
REM Installing Omin via the Klugey method.
:: FIXME: Call python from its root.

:: Activating miniconda

:: set miniconda_activate="%programfiles%\Miniconda3\Scripts\activate.bat"
:: set anaconda_activate="%programfiles%\Anaconda3\Scripts\activate.bat"

:: if exist "%programfiles%\Miniconda3\Scripts\activate.bat" echo activating miniconda3 & "%programfiles%\Miniconda3\Scripts\activate.bat"

:: IF exist miniconda_activate( miniconda_activate & echo activating Miniconda3... ) ELSE ( anaconda_activate & echo activating Anaconda3... )

:: REM Upgrading conda if needed...
:: conda upgrade -y conda

:: REM Uninstalling any old versions of omin.
:: pip uninstall -y omin

:: REM Installing omin...
:: python setup.py install

:: omin install files left behind.
:: rmdir build /s /q
:: rmdir dist /s /q
:: rmdir omin.egg-info /s /q

:: @echo I'M FINISHED!

timeout 35
