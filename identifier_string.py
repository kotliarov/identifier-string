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

NAMESPACES = {"x": "urn:hl7-org:v3"}

class Chains(object):
    """ SPL document protein chains.
    """
    xpath_moiety = "/x:moiety[x:code[@code=\"C118424\"]"
    xpath_local_id = "./x:partMoiety/x:id/@extension"
    xpath_aa = "./x:subjectOf/x:characteristic[x:code[@code=\"C103240\"]]/x:value[@media-type=\"application/x-aa-seq\"]"

    def __init__(self, doc):
        """ 
        """
        self._counter = 0
        self._chains = []
        self._load(doc)


    def _load(self, doc):
        """ Load chains defined in the SPL XML document.
        :doc: SPL XML DOM object
        """
        nodes = doc.xpath(self.xpath_moiety, namesoaces=NAMESPACES)
        for node in nodes:
            local_id = node.xpath(self.xpath_localid, namespaces=NAMESPACES)
            if not local_id:
                raise SPLDocumentError("local id not found")

            sequence = node.xpath(self.xpath_aa, namespaces=NAMESPACES)
            if not sequence:
                raise SPLDocumentError("Polypeptide chain AA sequence not found")


class ProteinSplModel(object):
    def __init__(self, xmldoc):
        pass

    def accept(self, visitor):
        """
        :visitor: Vistor object
        """
        visitor.visit("chains", self.chains)
        visitor.visit("polymers", self.polymers)
        visitor.visit("modifications", self.modifications)


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
        offset, var_name = name(offset)
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
        offset, func_name = name(offset)
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



class SPLDocumentError(Exception):
    def __init__(self, message):
        super().__init__(message)


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
    spl_doc = make_dom(docpath)
    identifier = IdentifierStringTemplate(Templates(Rules))
    ##model = make_document_model(spl_doc)
    ##model.accept(identifer)
    mock_model = MockModel()
    mock_model.accept(identifier)
    print(identifier.to_string(Rules["protein_identifier"]["str"]))
    

