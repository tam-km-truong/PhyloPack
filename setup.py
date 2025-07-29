from setuptools import setup, find_packages

setup(
    name="phylopack",
    version="0.1",
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'phylopack=phylopack.phylopack:main',
        ]
    },
    install_requires=[],
)