from setuptools import setup, find_packages

from ec2s import __version__

if __name__ == "__main__":
    setup(
        name="ec2s",
        version=__version__,
        packages=find_packages(),
        author="Wojciech Kusa",
        author_email="wojciech.kusa@tuwien.ac.at",
        description="Systematic reviews datasets",
        license="Apache 2.0",
    )
