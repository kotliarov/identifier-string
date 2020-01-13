import unittest
from idstring.identifier_string import IdentifierStringTemplate, Templates, TextGenerator


class TestTextGenerator(unittest.TestCase):

    def test_generator(self):
        gen = TextGenerator("/chains={{ chains }}")
        self.assertEqual(2, len(gen.nodes))
        self.assertEqual("/chains=", gen.nodes[0].to_string({}))
        self.assertEqual(["chains"], gen.variable_names())

    def test_idstring(self):
        identifier = IdentifierStringTemplate(Templates(Rules))
        model = MockModel()
        model.accept(identifier)
        self.assertEqual(identifier.to_string("/chains={{ chains }}/poly={{ polymers }}"), 
                        "/chains=c1:A;c2:AA;c3:ABC/poly=poly1:A;poly2:AA")


Rules = {
    "chain" : "{{ name }}:{{ value }}",
    "polymer" : "{{ name }}:{{ value }}",
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
    unittest.main()
