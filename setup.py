from setuptools import setup, find_packages
from os import path

setup(
    name="janas",
    description="description",
    url="https://github.com/fuyans/janas",
    author="Ian Fu",
    license="MIT",
    packages=find_packages(),
    install_requires=['pandas', 'scipy', 'matplotlib', 'numpy'],
    data_files=[ ('dat', ['project/dat/*'])]
)