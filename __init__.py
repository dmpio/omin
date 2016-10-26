__all__ = ["norm","up","excellor","vis","trunch","replicate","ceive","mitoCartaCall","fdr","experiment","norm_to_pool",
"norm_to_input","transform"]

from omin.up import *
from omin.norm import *
from omin.trunch import *
#import omin.mitoCartaCall
from omin.mitoCartaCall import *

from omin.replicate import *
from omin.ceive import *
from omin.excellor import *
from omin.fdr import *
from omin.experiment import *
from omin.norm_to_input import *
from omin.norm_to_pool import *
import omin.vis
import omin.transform

import os
import pickle
###LOADING mod_dict###
this_dir, this_filename = os.path.split(__file__)
mpd = "mod_dict.pickle"
DATA_PATH = os.path.join(this_dir, mpd)
f = open(DATA_PATH,"rb")
mod_dict = pickle.load(f)
