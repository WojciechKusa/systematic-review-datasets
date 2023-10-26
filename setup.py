from setuptools import setup, find_packages

from csmed import __version__

if __name__ == "__main__":
    setup(
        name="csmed",
        version=__version__,
        packages=find_packages(),
        author="Wojciech Kusa",
        author_email="wojciech.kusa@tuwien.ac.at",
        description="CSMeD - Citation Screening Meta-Dataset",
        license="Apache 2.0",
    )
