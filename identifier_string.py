""" layers.py: proof-of-a-concept implementation of a layered identifier string
               for SPL XML documents.
               Implementation consists of two components: model and view.

               1. View.
               
               View is responsible for presentation.
               Vuew is a collection of named string templates.
               Each named string template element
               - has collection of attributes and provides ability to set/get value of an attribute by name.
               - has a template string that defines element's representation.
               A renderer will substitute attributes' names in a template string with values of corresponding attributes.
               View may contain named group elements that serve as containers.

               Example:

               {
                    "protein_identifier" : {
                                            "attributes" : ["chains", "polymers", "modifications", "substitutions", "attachments"],
                                            "str": "/chains={{ chains|sort|join:";" }}/poly={{ polymers|sort|join:";" }}/mods={{ modifications|sort|join:";" }}"
                    },
                    "chain" : {
                        "attributes" : ["name", "value"],
                        "str" : "{{ name }}:{{ value }}"
                    },
                    "polymer" : {
                        "attributes" : ["name", "value"],
                        "str" : "{{ name }}:{{ value }}"
                    },
                    "modification" : {
                        "attributes" : ["name", "chain", "position", "type"],
                        "str" : "{{ name }}:{{ chain }}:{{ position }}:{{ type }}"
                    },
                    "substitution" : {
                        "attributes" : ["name", "modification", "polymer", "connection_point", "index"],
                        "str" : "{{ name }}:{{ modification }}:{{ polymer }}:{{ connection_point }}:{{ index }}"
                    }

                }

                2. Model.
                A model component represents SPL XML document instance.
                A visitor pattern will be used to to navigate SP document model's structure and collect
                data necessary to make an identifier string:
                
                ```
                SplDocumentModel* spl_document = new ProteinDocumentModel(file_path);
                StringTemplate identifier = templates.makeInstanceOf("protein_identifier");
                IdentifierStringVisitor visitor(identifier);
                spl_document->accept(visitor);
                string s = identifier.toString();

                ...

                void ProteinModel::accept(Vistor& visitor)
                {
                    visitor.visit(this->chains_);
                    visitor.visit(this->polymers_);
                    visitor.visit(this->modifications_);
                    visitor.visit(this->substitutions_);
                    visitor.visit(this->attachments_);
                    visitor.visit(this->quantities_);
                }

                ...

                class IdentifierStringVisitor : public Visitor
                {
                public:
                    IdentifierStringVisitor(StringTemplate& identifier_template)
                    : identifier_template_(identifier_template)
                    {
                    }
                
                    void visit(Chains& );
                    void visit(Polymers&);
                    ...
                
                private:
                    StringTemplate& identifier_template_;
                };

                void IdentifierStringVisitor::visit(Chains& chains) 
                {
                    for ( auto&& chain: chains ) {
                        StringTemplate t = templates.makeInstanceOf("chain");
                        t.setAttribute("name", chain.getName());
                        t.setAttribute("value", chanin.getValue();

                        identifier_template_.setAttribute("chains", t);
                    }
                }

                ```
"""

import os
import sys
import json
from hashlib import md5
from collections import defaultdict

from lxml import etree



class SplDocument(object):
    NAMESPACES = {"x": "urn:hl7-org:v3"}
    SPLDescriptor = {
        "document" : {
            "xpath"       : "/x:document[x:code[@code = '64124-1']]",
            "uri"         : "document/code[code=64124-1]",
            "content"     : ["section"],
            "mandatory"   : True,
            "cardinality" : 1
        },
        "section" : {
            "xpath"    : "./x:component/x:structuredBody/x:component/x:section[x:code[@code='48779-3']]",
            "uri"      : "section/code[code=48779-3]",
            "content"  : ["substance-main", "substance-other"],
            "mandatory": True,
            "cardinality" : 1
        },
        "substance-main" : {
            "xpath"      : "./x:subject/x:identifiedSubstance/x:identifiedSubstance[x:code[@codeSystem = '2.16.840.1.113883.4.9']]",
            "mandatory": True,
            "cardinality": 1
        },
        "substance-other" : {
            "xpath"      : "./x:subject/x:identifiedSubstance/x:identifiedSubstance[x:code[@codeSystem != '2.16.840.1.113883.4.9']]",
            "mandatory" : False
        },
    }
    
    def __init__(self, filepath):
        with open(docpath, 'r') as src:
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
        """ Return secrion element
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


class SplModelProtein(object):
    def __init__(self, xmldoc):
        self.chains = Chains(xmldoc)
        self.polymers = Polymers(xmldoc)
        self.modifications = []

    def accept(self, visitor):
        """
        :visitor: Vistor object
        """
        visitor.visit("chains", self.chains)
        visitor.visit("polymers", self.polymers)
        visitor.visit("modifications", self.modifications)


class Chains(object):
    """ SPL document protein chains.
    """
    xpath_moiety = "./x:moiety[x:code[@code=\"C118424\"]]"
    xpath_localid = "./x:partMoiety/x:id/@extension"
    xpath_aa = "./x:subjectOf/x:characteristic[x:code[@code=\"C103240\"]]/x:value[@mediaType=\"application/x-aa-seq\"]/text()"

    def __init__(self, doc):
        """ 
        """
        self._counter = 0
        self._pos = 0
        self.chains = []
        self._load(doc)

    def _load(self, doc):
        """ Load chains defined in the SPL XML document.
        :doc: SPL XML DOM object
        """
        substance = doc.substance()
        nodes = substance.xpath(self.xpath_moiety, namespaces=doc.NAMESPACES)
        for node in nodes:
            local_id = node.xpath(self.xpath_localid, namespaces=doc.NAMESPACES)
            if not local_id:
                raise SPLDocumentError("local id not found")

            value = node.xpath(self.xpath_aa, namespaces=doc.NAMESPACES)
            if not value:
                raise SPLDocumentError("Polypeptide chain AA sequence not found")
            self.chains.append(Chain(local_id[0], value[0]))
        self.chains = sorted(self.chains, key=lambda x: (x.value, x.local_id))
        for index, chain in enumerate(self.chains):
            chain.name = "chain{}".format(self._counter)
            self._counter += 1
        self._lookup = {x.local_id: x for x in self.chains}

    def __getitem__(self, key):
        """ Return chain descriptor by local id
        """
        return self._lookup[key]

    def __iter__(self):
        self._pos = 0
        return self

    def __next__(self):
        if self._pos < len(self.chains):
            i = self._pos
            self._pos += 1
            return self.chains[i]
        else:
            raise StopIteration()

class Chain(object):
    def __init__(self, local_id, value):
        self.local_id = local_id
        self.value = value
        self.name = None


class Polymers(object):
    """ SPL document polymers / irregular AA.
    """
    xpath_localid = "./x:code/@code"

    def __init__(self, doc):
        """ 
        """
        self._counter = 0
        self._pos = 0
        self.polymers = []
        self._load(doc)


    def _load(self, doc):
        """ Load polymers / irrgegula AA moleculs defined in the SPL XML document (other substance(s)).
        :doc: SPL XML DOM object
        """
       
        def read_code(subject):
            """ Reurn aux substance code"""
            code = subject.xpath(self.xpath_localid, namespaces=doc.NAMESPACES)
            if len(code) != 1:
                raise SPLDocumentError("Aux substance code not found")
            return code[0]
        
        def get_moiety(subject):
            """ Return moiety that representis subject's chemichal structure.
            """
            moiety_code = sub.xpath("./x:asSpecializedKind/x:generalizedMaterialKind/x:code/@code", namespaces=doc.NAMESPACES)
            if len(moiety_code) != 1:
                raise SPLDocumentError("Moiety code not found")
            moiety = subject.xpath("./x:moiety[x:partMoiety/x:code[@code=\"{}\"]]".format(moiety_code[0]),
                               namespaces=doc.NAMESPACES)
            if len(moiety) != 1:
                raise SPLDocumentError("Moiety \"{}\" not found".format(moiety_code[0]))
            return moiety[0]
    
        def get_connection_points(subject):
            """
            """
            class ConnectionPoint(object):
                def __init__(self, amino_group, carboxyl_group):
                    self.amino_group = amino_group
                    self.carboxyl_group = carboxyl_group

                def to_string(self):
                    return "N{}C{}".format(self.amino_group, self.carboxyl_group)
            
            points = []
            nodes = subject.xpath("./x:moiety[x:code[@code=\"C118427\"]]", namespaces=doc.NAMESPACES)
            for node in nodes:
                positions = node.xpath("./x:positionNumber/@value", namespaces=doc.NAMESPACES)
                points.append(ConnectionPoint(positions[0], positions[1]))
            return points


        subjects = doc.substance_other()
        if subjects is not None:
            for sub in subjects:
                code = read_code(sub)
                moiety = get_moiety(sub)
                conn_points = get_connection_points(sub)
                value = get_chem_structure(moiety, None)
                quantity = get_quantity(moiety)
                self.polymers.append(Polymer(code, value, conn_points, quantity))
        self.polymers = sorted(self.polymers, key=lambda x: (x.value, x.code))
        self._lookup = {x.code: x for x in self.polymers}
        for index, x in enumerate(self.polymers):
            x.name = "poly{}".format(self._counter)
            self._counter += 1


    def __getitem__(self, key):
        """ Return chain descriptor by local id
        """
        return self._lookup[key]

    def __iter__(self):
        self._pos = 0
        return self

    def __next__(self):
        if self._pos < len(self.polymers):
            i = self._pos
            self._pos += 1
            return self.polymers[i]
        else:
            raise StopIteration()


class Polymer(object):
    def __init__(self, code, value, conn_points, quantity):
        self.code = code
        self.value = value
        self.conn_points = conn_points
        self.quantity = quantity
        self.name = None
    
    @property
    def connection_points(self):
        """ Return connection points
        """
        return ",".join([x.to_string() for x in self.conn_points])


class Modifcations(object):
    
    xpath_moiety = "./x:moiety[x:code[@code=\"C118425\"]]/partMoiety" 

    def __init__(self, doc):
        """ 
        """
        self._counter = 0
        self._pos = 0
        self.polymers = []
        self._load(doc)

    def _load(self, doc):
        substance = doc.substance()
        nodes = substance.xpath(self.xpath_moiety, namespaces=doc.NAMESPACES)
        for node in nodes:
            code = node.xpath("./x:code/@code", namespaces=doc.NAMESPACES)
            if len(code) != 1:
                raise SPLDocumentError("Moiety substance code not found")
            bonds = node.xpath("./x:bond[x:code[@code=\"C118426\"]]") # AA substitutions
            if bonds:
                self.mods.append(Substitution(code, bonds))
            bonds = node.xpath("./x:bond[x:code[@code=\"C14050\"]]") # Attachments
            if bonds:
                self.mods.append(Attachment(code, bonds))

    def __iter__(self):
        self._pos = 0
        return self

    def __next__(self):
        if self._pos < len(self.mods):
            i = self.pos_
            self.pos_ += 1
            return self.mods[i]
        else:
            raise StopIteration() 


class AminoAcidSubstitution(object):
    def __init__(self, bonds, irreg_aa, chains_lookup):
        self.bonds = []
        self._load(bonds, irreg_aa, chains_lookup)

    def _load(self, bonds, irreg_aa, chains_lookup):
        """
        """
        for bond in bonds:
            local_id = bond.xpath("./x:distalMoiety/x:id/@extension")[0]
            self.chain = chain_lookup(local_id)
            positions = bond.xpath("./x:positionNumber/@value")
            if len(positions) != 2:
                raise SPLDocumentError("Expecting two position per bond")
            positions = map(int, positions)
            self.bonds.append(Bond(irreg_aa, positions[0], chain, positions[1]))
        self.bonds = sorted(self.bonds, key=lambda x: (x.irreg_aa.name, x.irreg_aa_pos, x.chain.name, x.chain_pos))

    def to_string(self):
        return "sub:" + ",".join([x.to_string() for x in self.bonds])

    @property
    def value(self):
        return self.to_string()


class Bond(object):
    def __init__(self, irreg_aa, irreg_aa_pos, chain, chain_pos):
        self.irreg_aa = irreg_aa
        self.irreg_aa_pos = irreg_aa_pos
        self.chain = chain
        self.chain_pos = chain_pos

    def to_string(self):
        return "{}:{}:{}:{}".format(self.irreg_aa.name,
                                    self.irreg_aa_pos,
                                    self.chain.name,
                                    self.chain_pos)

CHEMICAL_STRUCT = [
    ("x-inchi-key", "./x:subjectOf/x:characteristic[x:code[@code=\"C103240\"]]/x:value[@mediaType=\"application/x-inchi-key\"]/text()"),
    ("x-inchi", "./x:subjectOf/x:characteristic[x:code[@code=\"C103240\"]]/x:value[@mediaType=\"application/x-inchi\"]/text()"),
    ("x-mdl-molfile", "./x:subjectOf/x:characteristic[x:code[@code=\"C103240\"]]/x:value[@mediaType=\"application/x-mdl-molfile\"]/text()"),
    ("x-aa=seq", "./x:subjectOf/x:characteristic[x:code[@code=\"C103240\"]]/x:value[@mediaType=\"application/x-aa-seq\"]/text()"),
    ("x-na-seq", "./x:subjectOf/x:characteristic[x:code[@code=\"C103240\"]]/x:value[@mediaType=\"application/x-na-seq\"]/text()"),
]

def get_chem_structure(moiety, mediaType, namespaces=SplDocument.NAMESPACES):
    """ Return chemical struct balue.
        :moiety: xml dom element 
        :mediaType: concrete mediaType or None if any
    """
    def get_value(query):
        nodes = moiety.xpath(query, namespaces=namespaces)
        if nodes:
            return nodes[0]
        else:
            return None

    for media, query in CHEMICAL_STRUCT:
        if mediaType is None or mediaType == media:
            value = get_value(query)
            if value is not None or mediaType is not None:
                return value
    return None

def get_quantity(moiety, namespaces=SplDocument.NAMESPACES):
    """
    """
    class Quantity(object):
        def __init__(self, num, denom, unit):
            self.num = num
            self.denom = denom
            self.unit = unit

        def to_string(self):
            return "{}:{}:{}".format(self.num, self.denom. self.unit)

    class QuanityRange(object):
        def __init__(self, low, low_inclusive, high, high_inclusive,  denom, unit):
            self.low = low
            self.low_incl = low_inclusive
            self.high = high
            self.high_incl = high_inclusive
            self.denom = denom
            self.unit = unit

        def to_string(self):
            return "{}{},{}{}:{}:{}".format("[" if self.low_incl else "(",
                                            self.low,
                                            self.high,
                                            "]" if self.high_incl else ")",
                                            self.denom. 
                                            self.unit)


    def is_inclusive(el):
        rc = True
        if "inclusive" in el:
            if "false" == el["inclusive"]:
                rc = False
        return rc

    numerator = moiety.xpath("./x:quantity/x:numerator", namespaces=namespaces)[0]
    denominator = moiety.xpath("./x:quantity/x:denominator", namespaces=namespaces)[0]
    
    attributes = denominator.attrib
    unit = attributes["unit"]
    denom_value = attributes["value"]
    
    attributes = numerator.attrib
    if "value" in attributes:
        num_value = attributes["value"]
        return Quantity(num_value, denom_value, unit)
    else: # range
        low = numerator.xapth("./x:low")[0]
        high = numerator.xapth("./x:high")[0]
        return QuantityRange(low.attrib["value"],
                             is_inclusive(low),
                             high.attrib["value"],
                             is_inclusive(high),
                             denom_value,
                             unit)


class IdentifierStringTemplate(object):
    """
    """
    def __init__(self, templates):
        self.templates = templates
        self.context = {}

    def visit(self, name, target):
        """
        :name: object name
        :target: 
        """
        if name == "chains":
            self.visit_chains(target)
        elif name == "polymers":
            self.visit_polymers(target)
        elif name == "modifications":
            self.visit_modifications(target)
        else:
            pass

    def visit_chains(self, target):
        """
        """
        self.context["chains"] = []
        for chain in target:
            t = self.templates.make_instance_of("chain")
            t.load(chain)
            self.context["chains"].append(t)

    def visit_polymers(self, target):
        """
        """
        self.context["polymers"] = []
        for polymer in target:
            t = self.templates.make_instance_of("polymer")
            t.load(polymer)
            self.context["polymers"].append(t)

    def visit_modifications(self, target):
        """
        """
        self.context["modifications"] = []
        for modification in target:
            t = self.templates.make_instance_of("modification")
            t.load(modification)
            self.context["modifications"].append(t)

    def to_string(self, template_str):
        """
        """
        gen = TextGenerator(template_str)
        return gen.to_string(self.context)


class TextGenerator(object):
    """ Generates text according to specified template by
        substituting template's variables with values
        of correspondng variables stored in the context.
    """
    def __init__(self, template_str):
        """
        :template_str: text to emit. Variables in the text defined as {{ varname }}.
        """
        self.nodes = parse(template_str)
        
    def to_string(self, context):
        """
        :context: dictionary with named variables
        """
        return ''.join([x.to_string(context) for x in self.nodes])


class TextGeneratorElementString(object):
    def __init__(self, value):
        self.value = value

    def to_string(self, context):
        return self.value


class TextGeneratorElementVariable(object):
    def __init__(self, variable):
        self.variable = variable

    def to_string(self, context):
        """
        """
        try:
            item = context[self.variable.name]
        except KeyError:
            return "{{ {} }}".format(self.variable.name)
        
        if isinstance(item, str):
            return item
        elif isinstance(item, list):
            coll = [x.to_string() for x in item]
            return ";".join(coll)
        else:
            return item.to_string()


def parse(input_str):
    """
    """
    class Variable(object):
        def __init__(self, name, transforms):
            self.name = name
            self.transforms = transforms

    class Transform(object):
        def __init__(self, func_name, params):
            self.func_name = func_name
            self.params = params

    nodes = []
    pos = 0
    N = len(input_str)

    def parse_variable(offset):
        """ Parse variable: "{{ " <variable> " }}"
        """
        offset = match_pattern("{{ ", offset)
        offset, var = match_variable(offset)
        offset = match_pattern(" }}", offset)
        return offset, var

    def match_pattern(pattern, offset):
        """ Match specified pattern and advance
            position within input buffer.
        """
        n = len(pattern)
        if offset + n <= N:
            if pattern == input_str[offset: offset + n]:
                return offset + n
        raise ParserError()

    def match_variable(offset):
        """ Match variable pattern and 
            retrun Variable object.
            Variable pattern: name(|transform)*
        """
        offset, var_name = match_name(offset)
        transforms = []
        while input_str[offset] == '|':
            offset, transform = match_transformation(offset + 1)
            transforms.append(transform)
        return offset, Variable(var_name, transforms)

    def match_transformation(offset):
        """ Match transformation pattern.
            Return transformation: function name, parameters.

            Transformation pattern: "|func_name(:parameter(,parameter)*)?
        """
        offset, func_name = match_name(offset)
        params = []
        if input_str[offset] == ':':
            offset, param = match_parameter(offset)
            params.append(param)
            while input_str[offset] == ',':
                offset, param = match_parameter(offset+1)
                params.append(param)
        return offset, Transform(func_name, params)

    def match_name(offset):
        """ Return name token.
            Name token is a sequence of characters up to 
            a first occurrence of characters '|' or ' '.
        """
        text = []
        while offset < N:
            c = input_str[offset]
            if c == '|' or c == ' ':
                break
            text.append(c)
            offset += 1
        return offset, ''.join(text)

    def match_parameter(offset):
        """ Return parameter value
        """
        text = []
        quoted = False

        while offset < N:
            c = input_str[offset]
            if c == '"':
                if quoted:
                    offset += 1
                    break
                elif not text:
                    quoted = True
                    offset += 1
            elif c == '|' or c == ' ' or c == ',':
                if not quoted:
                    break
                text.append(c)
                offset += 1
            else:
                text.append(c)
                offset += 1
        return offset, ''.join(text)


    nodes = []
    while pos < N:
        offset = input_str[pos: ].find("{{ ")
        if offset != -1:
            nodes.append(TextGeneratorElementString(input_str[pos: pos+offset]))
            pos, variable = parse_variable(pos+offset)
            nodes.append(TextGeneratorElementVariable(variable))
        else:
            nodes.append(TextGeneratorElementString(input_str[pos:]))
            pos = N
    return nodes


class StringTemplate(object):
    """
    """
    def __init__(self, template):
        self.attributes = template["attributes"]
        self.string_template = template["str"]
        self._context = {}

    def load(self, target):
        """ Load context with value of target's attributes.
        """
        for name in self.attributes:
            try:
                value = getattr(target, name)
            except AttributeError:
                pass
            finally:
                self._context[name] = value

    def __getattr__(self, name):
        """ Return value of specified attribute.
        """
        return self._context[name]

    def to_string(self):
        """ Generate text according to a template.
            Use context to replace template's variables with
            concete values.
        """
        gen = TextGenerator(self.string_template)
        return gen.to_string(self._context)


def make_dom(docpath):
    """
    REturn XML DOM object.
    :docpath: path tp SPL XML document
    """
    with open(docpath, 'r') as src:
        doc = src.read()
    return etree.fromstring(doc)


class Templates(object):
    def __init__(self, rules):
        self.rules = rules

    def make_instance_of(self, name):
        if name in self.rules:
            return StringTemplate(self.rules[name].copy())
        else:
            raise KeyError(name)

Rules = {
    "protein_identifier" : {
        "attributes" : ["chains", "polymers", "modifications", "substitutions", "attachments"],
        "str": '/chains={{ chains|sort|join:";" }}/poly={{ polymers|sort|join:";" }}/mods={{ modifications|sort|join:";" }}'
    },
    "chain" : {
        "attributes" : ["name", "value"],
        "str" : "{{ name }}:{{ value }}"
    },
    "polymer" : {
        "attributes" : ["name", "value", "connection_points"],
        "str" : "{{ name }}:{{ value }}:{{ connection_points }}"
    },
    "modification" : {
        "attributes" : ["name", "chain", "position", "type"],
        "str" : "{{ name }}:{{ chain }}:{{ position }}:{{ type }}"
    },
    "substitution" : {
        "attributes" : ["name", "modification", "polymer", "connection_point", "index"],
        "str" : "{{ name }}:{{ modification }}:{{ polymer }}:{{ connection_point }}:{{ index }}"
    }
}


class MockModel(object):
    def __init__(self):
        pass

    def accept(self, visitor):
        """
        """
        chains = [MockChain("c1", "A"), MockChain("c2", "AA"), MockChain("c3", "ABC")]
        visitor.visit("chains", chains)

        polymers = [MockPoly("poly1", "A"), MockPoly("poly2", "AA")]
        visitor.visit("polymers", polymers)


class MockChain(object):
    def __init__(self, name, value):
        self.name = name
        self.value = value

class MockPoly(object):
    def __init__(self, name, value):
        self.name = name
        self.value = value


if __name__ == '__main__':
    docpath = sys.argv[1]
    identifier = IdentifierStringTemplate(Templates(Rules))
    model = SplModelProtein(SplDocument(docpath))
    model.accept(identifier)
    #mock_model = MockModel()
    #mock_model.accept(identifier)
    print(identifier.to_string(Rules["protein_identifier"]["str"]))
    


