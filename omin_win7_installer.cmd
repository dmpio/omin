:: Upgrade conda if needed.
conda upgrade -y conda

:: Uninstall any old verision of omin.
pip uninstall -y omin

:: Install omin.
python setup.py install

:: omin install files left behind.
rmdir build /s /q
rmdir dist /s /q
rmdir omin.egg-info /s /q

@echo I'M FINISHED!

timeout 5
