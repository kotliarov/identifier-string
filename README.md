## Identifier String for SPL XML documents.

In order to be able to specify structure of an identifier string in a configuration file
and produce an instance of identifier string for a concrete SPL XML document we will use 
a string template engine to enforce model-view separation.
This means that all logic and SPL document data model is kept in code and all the output
text in the templates.

### String template.

String template is a chunks of text and expressions enclosed in double-curly braces, e.g. {{ expression }}.
String template treats everything outside the expressions as text to emit.
The idea behind constructing output is that we create a template and inject it with variables.
Nested template hierarchy can be created by injecting templates as variables into other templates. 

1. A template is a text string.
2. A template string expresses presentation, not a program logic.
3. A template contains variables, which get replaced with values when the template is rendered/evaluated.

Examples of string templates::
```
    "protein_identifier" : "/chains={{ chains }}/poly={{ polymers }}/subs={{ substitutions }}",
    "chain" : "{{ name }}:{{ value }}:{{ quantity }}",
    "polymer" : "{{ name }}:{{ value }}:{{ connection_points }}:{{ quantity }}",
    "substitution" : "{{ name }}:{{ chain }}:{{ position }}:{{ polymer }}:{{ connection_point }}",
    "attachment" :  "{{ chain }}:{{ position }}:{{ glycan }}",
```


### Variables

Template variables are defined by the context dictionary passed to the template. 
They are passed in by the application. Variables may have attributes. 
What attributes a variable has depends on the application providing that variable.

There is a delimiter to print to the template output: {{ variable }}. 
Outer double-curly braces are not part of the variable, but the print statement.


### Model.

A model component represents SPL XML document instance.
A visitor pattern will be used to to navigate SPL document model's structure, collect
data and pass it to corresponding templates.

```python
# Make a document model.
spl_document = ProteinSplModel(file_path)

# Make an instance of identifier string template.
template = templates.makeInstanceOf("protein_identifier")

# Make an instance of identifier string visitor that "knows" document structure. 
visitor = IdentifierStringVisitor()

# Visit document' elements.
spl_document.accept(visitor);

# Generate output according to the template.
s = identifier.toString(template)
print(s)

...

# Data model.

class ProteinSplModel(object):
  ...
  def accept(self, visitor)
    visitor.visit("chains", self.chains_)
    visitor.visit("polymers", self.polymers_)
    visitor.visit("substitutions", self.substitutions_)
    visitor.visit("attachments", self.attachments_)

...
# Visitor

class IdentifierStringVisitor(object):
  def __iint__(self):
    self.context = {}

  def visit(self, name, target):
    if name == "chains"
      self.visit_chains(target)
    elif ...

  def visit_chains(self, chains):
    for chain in chains:
      t = templates.makeInstanceOf("chain")
      t.setAttribute("name", chain.getName())
      t.setAttribute("value", chanin.getValue()
      self.context["chains"].append(t)

  def toString(self, template):
    # Ask string template to evaluate the template
    return template.toString(self.context_)
```

### Install and try out.
```sh
git clone https://github.com/kotliarov/identifier-string.git
virtualenv --python=python3.6 env3
source env3/bin/activate
cd identifier-string
pip install .
python -m idstring ./protein.xml
```

Windows:
```
git clone https://github.com/kotliarov/identifier-string.git
virtualenv  env3
env3\Scripts\activate.bat
cd identifier-string
pip install .
python -m idstring ./protein.xml
```
