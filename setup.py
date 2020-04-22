import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="atlas-metadata-validator",
    version="1.3.0",
    author="Anja Füllgrabe",
    author_email="anjaf@ebi.ac.uk",
    description="A MAGE-TAB validator for Expression Atlas and Single Cell Expression Atlas",
    license="Apache Software License",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ebi-gene-expression-group/atlas-metadata-validator",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
    scripts=[
        'atlas_validation.py'
    ],
    include_package_data=True,
    install_requires=[
        "requests>=2.20.1"
    ],
    python_requires=">=3.6",
)
