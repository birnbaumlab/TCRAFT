from setuptools import setup, find_packages

def parse_requirements(filename):
    with open(filename, "r") as f:
        return [line.strip() for line in f if line.strip() and not line.startswith("#")]

setup(
    name="TCRAFT",
    version="1.0.0",
    packages=find_packages(),
    include_package_data=True,
    package_data={
        "TCRAFT.references": ['*'],
    },
    install_requires=parse_requirements("requirements.txt"),
    entry_points={
        "console_scripts": [
            "TCRAFT-generate=TCRAFT.generate_cdr3_oligos:main",
            "TCRAFT-validate=TCRAFT.validate_cdr3_oligos:main"
        ],
    },
    author="Rachit Mukkamala",
    description="Software package to generate and validate CDR3 oligonucleotide pools for TCRAFT (TCR Rapid Assembly for Functional Testing)",
    url="https://github.com/birnbaumlab/TCRAFT",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
)
