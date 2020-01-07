""" SPL document.

"""
from lxml import etree


class SplDocument(object):
    NAMESPACES = {"x": "urn:hl7-org:v3"}
    SPLDescriptor = {
        "document": {
            "xpath": "/x:document[x:code[@code = '64124-1']]",
            "uri": "document/code[code=64124-1]",
            "content": ["section"],
            "mandatory": True,
            "cardinality": 1
        },
        "section": {
            "xpath": "./x:component/x:structuredBody/x:component/x:section[x:code[@code='48779-3']]",
            "uri": "section/code[code=48779-3]",
            "content": ["substance-main", "substance-other"],
            "mandatory": True,
            "cardinality": 1
        },
        "substance-main": {
            "xpath": "./x:subject/x:identifiedSubstance/x:identifiedSubstance[x:code[@codeSystem = '2.16.840.1.113883.4.9']]",
            "mandatory": True,
            "cardinality": 1
        },
        "substance-other": {
            "xpath": "./x:subject/x:identifiedSubstance/x:identifiedSubstance[x:code[@codeSystem != '2.16.840.1.113883.4.9']]",
            "mandatory": False
        },
    }

    def __init__(self, file_path):
        with open(file_path, 'r') as src:
            doc = src.read()
        self.dom = etree.fromstring(doc)
        self.document_ = None
        self.section_ = None
        self.substance_ = None
        self.substance_other_ = None

    def document(self):
        """ Return document element.
        """
        if self.document_ is None:
            desc = self.SPLDescriptor["document"]
            nodes = self.dom.xpath(desc["xpath"], namespaces=self.NAMESPACES)
            if len(nodes) != 1:
                raise SPLDocumentError("Document element must be present and unique")
            self.document_ = nodes[0]
        return self.document_

    def section(self):
        """ Return section element
        """
        if self.section_ is None:
            desc = self.SPLDescriptor["section"]
            nodes = self.document().xpath(desc["xpath"], namespaces=self.NAMESPACES)
            if len(nodes) != 1:
                raise SPLDocumentError("Section element must be present and unique")
            self.section_ = nodes[0]
        return self.section_

    def substance(self):
        """ Return main substance element.
        """
        if self.substance_ is None:
            desc = self.SPLDescriptor["substance-main"]
            nodes = self.section().xpath(desc["xpath"], namespaces=self.NAMESPACES)
            if len(nodes) != 1:
                raise SPLDocumentError("Main substance element must be present and unique")
            self.substance_ = nodes[0]
        return self.substance_

    def substance_other(self):
        if self.substance_other_ is None:
            desc = self.SPLDescriptor["substance-other"]
            nodes = self.section().xpath(desc["xpath"], namespaces=self.NAMESPACES)
            self.substance_other_ = nodes
        return self.substance_other_

    def select(self, base, query):
        """
        """
        return base.xpath(query, namespaces=self.NAMESPACES)


class SPLDocumentError(Exception):
    def __init__(self, message):
        super().__init__(message)