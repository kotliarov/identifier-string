import sys

from idstring.identifier_string import IdentifierStringTemplate, Templates, Rules
from idstring.model import SplModelProtein
from idstring.spl import SplDocument


if __name__ == "__main__":
    docpath = sys.argv[1]
    identifier = IdentifierStringTemplate(Templates(Rules))
    model = SplModelProtein(SplDocument(docpath))
    model.accept(identifier)
    print(identifier.to_string(Rules["protein_identifier"]["str"]))

