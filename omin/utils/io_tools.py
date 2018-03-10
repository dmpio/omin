# -*- coding: utf-8 -*-
"""omin.io_tools

Provides
--------
Tools for file input output handling.
"""
# Copyright 2018 James Draper, Paul Grimsrud, Deborah Muoio, Colette Blach, Blair Chesnut, and Elizabeth Hauser.

import os
import re
import shutil


class IOTools(object):
    """Tools for file handling and dir building.

    """

    @staticmethod
    def sanitize_file_path(path):
        """Return an os safe file string.
        """
        result = re.sub('[^0-9a-zA-Z]+', '_', path)
        return result

    @classmethod
    def file_path_here(cls, path, ext=None):
        """Return a sanitized file path for the current working directory.
        """
        path = cls.sanitize_file_path(path)

        if ext is not None:
            path = '.'.join([path, ext])

        here = os.path.abspath('.')
        result = os.path.join(here, path)
        return result

    @staticmethod
    def mkdir(directory):
        """Create a directory.

        Parameters
        ----------
        directory : str
        """
        # FIXME: This class is weak :/
        # FIXME: Should just use shutil
        assert type(directory) == str
        directory = directory.replace(" ", "_")

        if os.path.isdir(directory):
            pass
        else:
            print("Created dir:", directory)
            os.mkdir(directory)


class UserProfile(object):
    """Creates user profile files for configuration and databases.
    """
    user_profile = os.path.expanduser("~")
    omin_profile_dir = os.path.join(user_profile, ".omin")
    omin_database_profile_dir = os.path.join(omin_profile_dir, "databases")
    # Get this dir as a string.
    this_dir, _ = os.path.split(__file__)
    this_dir, _ = os.path.split(this_dir)
    # Create the path string.
    carta_src = os.path.join(this_dir, 'databases', 'mitocarta', 'complete_mitocarta_2.p')
    carta_dst = os.path.join(omin_database_profile_dir, 'complete_mitocarta_2.p')


    @classmethod
    def _create_omin_profile_dir(cls):
        IOTools.mkdir(cls.omin_profile_dir)


    @classmethod
    def _create_omin_database_profile_dir(cls):
        IOTools.mkdir(cls.omin_database_profile_dir)


    @classmethod
    def _copy_mitocarta(cls):
        if os.path.exists(cls.carta_dst):
            pass
        else:
            print("Copying the MitoCarta2.0 Database to:", cls.carta_dst)
            shutil.copyfile(cls.carta_src, cls.carta_dst)
            # shutil copy


    @classmethod
    def _create_profile_dirs(cls):
        cls._create_omin_profile_dir()
        cls._create_omin_database_profile_dir()
        cls._copy_mitocarta()


UserProfile._create_profile_dirs()
