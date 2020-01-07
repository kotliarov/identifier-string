from spl import SplDocument, SPLDocumentError


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
            """ Return aux substance code"""
            code = subject.xpath(self.xpath_localid, namespaces=doc.NAMESPACES)
            if len(code) != 1:
                raise SPLDocumentError("Aux substance code not found")
            return code[0]

        def get_moiety(subject):
            """ Return moiety that represents subject's chemical structure.
            """
            moiety_code = sub.xpath("./x:asSpecializedKind/x:generalizedMaterialKind/x:code/@code",
                                    namespaces=doc.NAMESPACES)
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


class Modifications(object):
    xpath_moiety = "./x:moiety[x:code[@code=\"C118425\"]]/partMoiety"

    def __init__(self, doc):
        """
        """
        self._counter = 0
        self._pos = 0
        self.mods = []
        self._load(doc)

    def _load(self, doc):
        substance = doc.substance()
        nodes = substance.xpath(self.xpath_moiety, namespaces=doc.NAMESPACES)
        for node in nodes:
            code = node.xpath("./x:code/@code", namespaces=doc.NAMESPACES)
            if len(code) != 1:
                raise SPLDocumentError("Moiety substance code not found")
            bonds = node.xpath("./x:bond[x:code[@code=\"C118426\"]]")  # AA substitutions
            if bonds:
                self.mods.append(Substitution(code, bonds))
            bonds = node.xpath("./x:bond[x:code[@code=\"C14050\"]]")  # Attachments
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
    ("x-inchi-key",
     "./x:subjectOf/x:characteristic[x:code[@code=\"C103240\"]]/x:value[@mediaType=\"application/x-inchi-key\"]/text()"),
    ("x-inchi",
     "./x:subjectOf/x:characteristic[x:code[@code=\"C103240\"]]/x:value[@mediaType=\"application/x-inchi\"]/text()"),
    ("x-mdl-molfile",
     "./x:subjectOf/x:characteristic[x:code[@code=\"C103240\"]]/x:value[@mediaType=\"application/x-mdl-molfile\"]/text()"),
    ("x-aa=seq",
     "./x:subjectOf/x:characteristic[x:code[@code=\"C103240\"]]/x:value[@mediaType=\"application/x-aa-seq\"]/text()"),
    ("x-na-seq",
     "./x:subjectOf/x:characteristic[x:code[@code=\"C103240\"]]/x:value[@mediaType=\"application/x-na-seq\"]/text()"),
]


def get_chem_structure(moiety, mediaType, namespaces=SplDocument.NAMESPACES):
    """ Return chemical struct value.
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
            return "{}:{}:{}".format(self.num, self.denom.self.unit)

    class QuanityRange(object):
        def __init__(self, low, low_inclusive, high, high_inclusive, denom, unit):
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
    else:  # range
        low = numerator.xapth("./x:low")[0]
        high = numerator.xapth("./x:high")[0]
        return QuantityRange(low.attrib["value"],
                             is_inclusive(low),
                             high.attrib["value"],
                             is_inclusive(high),
                             denom_value,
                             unit)
