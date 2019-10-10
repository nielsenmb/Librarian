import setuptools

exec(open('asterosearch/version.py').read())

setuptools.setup(
    name="asterosearch",
    version=__version__,
    author="Martin Nielsen",
    author_email="m.b.nielsen.1@bham.ac.uk",
    description="A package for querying multiple stellar catalogs",
    long_description=open("README.md").read(),
    packages=setuptools.find_packages(),
    install_requires=open('requirements.txt').read().splitlines(),
    include_package_data=True,   
)
