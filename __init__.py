# -*- coding: utf-8 -*-
import os
import pickle

from . import normalize


__docformat__ = 'restructuredtext'

# from omin.info import __doc__
# from omin.ceive import *
# from omin.norm import *
# from omin.trunch import *
# from omin.mitoCartaCall import *
# from omin.replicate import *
# from omin.excellor import *
# from omin.fdr import *
# from omin.experiment import *
# from omin.norm_to_input import *
# from omin.norm_to_pool import *
# import omin.vis
# from omin.comparison import *

# === LOADING mod_dict ===

this_dir, this_filename = os.path.split(__file__)
mpd = "mod_dict.pickle"
DATA_PATH = os.path.join(this_dir, mpd)
f = open(DATA_PATH, "rb")
mod_dict = pickle.load(f)

# # === import all submodules ===
# import importlib
# import pkgutil
#
#
# def import_submodules(package, recursive=True):
#     """ Import all submodules of a module, recursively, including subpackages
#
#     :param package: package (name or actual module)
#     :type package: str | module
#     :rtype: dict[str, types.ModuleType]
#     """
#     if isinstance(package, str):
#         package = importlib.import_module(package)
#     results = {}
#
#     for loader, name, is_pkg in pkgutil.walk_packages(package.__path__):
#         full_name = package.__name__ + '.' + name
#         results[full_name] = importlib.import_module(full_name)
#         if recursive and is_pkg:
#             results.update(import_submodules(full_name))
#     return results
#
# import_submodules(__name__)
