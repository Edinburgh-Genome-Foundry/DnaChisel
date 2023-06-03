# This will try to import setuptools. If not here, it will reach for the embedded
# ez_setup (or the ez_setup package). If none, it fails with a message
try:
    from setuptools import setup
except ImportError:
    try:
        import ez_setup

        ez_setup.use_setuptools()
    except ImportError:
        raise ImportError(
            "DnaChisel could not be installed, probably because"
            " neither setuptools nor ez_setup are installed on"
            "this computer. \nInstall ez_setup "
            "([sudo] pip install ez_setup) and try again."
        )

from setuptools import setup, find_packages

exec(open("dnachisel/version.py").read())  # loads __version__

setup(
    name="dnachisel",
    version=__version__,
    author="Zulko",
    description="Optimize DNA sequences under constraints.",
    long_description=open("pypi-readme.rst").read(),
    license="MIT",
    url="https://github.com/Edinburgh-Genome-Foundry/DnaChisel",
    keywords=("DNA gene design codon optimization constraints " "synthetic biology"),
    packages=find_packages(exclude="docs"),
    include_package_data=True,
    scripts=["scripts/dnachisel"],
    install_requires=[
        "numpy",
        "Biopython",
        "proglog",
        "docopt",
        "flametree",
        "python_codon_tables",
    ],
    extras_require={
        "reports": [
            "pdf_reports",
            "sequenticon",
            "matplotlib",
            "dna_features_viewer",
            "pandas",
        ],
        "tests": [
            "pytest",
            "pytest-cov",
            "coveralls",
            "geneblocks",
            "genome_collector",
            "matplotlib",
            "primer3-py",
        ],
    },
)
