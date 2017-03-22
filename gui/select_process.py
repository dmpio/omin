# -*- coding: utf-8 -*-
import easygui as eg

def select_process():
    selected = eg.choicebox("Please pick a process:",
                            "Select Process",
                            ["RawData","PreProcess","Process"])
    return selected

