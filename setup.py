from setuptools import setup, find_packages

setup(name='idstring',
      version='1.0',
      description='SPL Document Id String',
      author='Alex Kotliarov',
      packages=find_packages('src'),
      package_dir={'': 'src'},
      install_requires=['lxml']
      )
