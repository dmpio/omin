# -*- coding: utf-8 -*-
"""
DEPRECATED: This will be eventually integrated into th econtainer class.
"""

# import guipyter
from ..core import guipyter
from xml.etree import cElementTree as et


class PDStudyTools(object):

    def __init__(self, file_name=None):
        self.file_name = file_name
        self.file_reader()
        self.file_to_tree()


    def file_reader(self):
        if self.file_name == None:
            self.file_name = guipyter.filedialog.askopenfilename()

        with open(self.file_name, 'rb') as f:
            # NOTE: ENCODED AS BYTES!
            data = f.read()
        self.decoded_data = data.decode("utf-8-sig")

    def file_to_tree(self):
        self.xml_tree = et.fromstring(self.decoded_data)

    @property
    def study_factors(self):
        """Return a dict with study factors as keys and factor options as values."""
        result = dict()
        for i in list(self.xml_tree.iter('Factor')):
            # print(i.attrib['Name'])
            sf_name = i.attrib['Name']
            if sf_name not in result:
                result[sf_name] = []
            for j in list(i.iter('FactorOption')):
                result[sf_name].append(j.attrib['Value'])

        return result
