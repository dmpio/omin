"Creates file in .ipython\profile_default\startup that adds package to path."

import os
from pathlib import PureWindowsPath

dir_path = os.path.dirname(os.path.realpath(__file__))
dir_path = os.path.join(dir_path)

dir_path = dir_path.encode('unicode_escape').decode('latin1')
usr_path = os.path.expanduser('~')
ipython_startup = os.path.join(usr_path, '.ipython', 'profile_default', 'startup')

code ="""
import sys
sys.path.append("{}")
""".format(dir_path)

if os.path.exists(ipython_startup):
    with open(os.path.join(ipython_startup, 'append_omin_to_path.py'), 'w') as f:
        f.write(code)
else:
    print('Could not find:', ipython_startup)
