from setuptools import setup, find_packages
import codecs
import os.path


def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), "r") as fp:
        return fp.read()


def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith("__version__"):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")


setup(
    name="protinter",
    author="Maxime Borry",
    version=get_version("protinter/__init__.py"),
    description="ProtInter : Protein interaction calculator",
    long_description=open("README.md").read(),
    url="https://github.com/maxibor/protinter",
    long_description_content_type="text/markdown",
    license="GNU-GPLv3",
    python_requires=">=3.6",
    install_requires=["click", "biopython"],
    packages=find_packages(include=["protinter"]),
    entry_points={"console_scripts": ["protinter = protinter.cli:cli"]},
)
