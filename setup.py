from setuptools import setup

setup(
   name='MiSDEED',
   version='1.0.0',
   description='Microbiome Synthetic Data Engine for Experimental Design ',
   author='Philippe Chlenski',
   author_email='pac@cs.columbia.edu',
   url='http://www.github.com/pchlenski/misdeed',
   packages=['misdeed'],
   package_dir={"misdeed": "src"},
   install_requires=['numpy', 'pandas', 'sklearn', 'scipy', 'matplotlib'],
)
