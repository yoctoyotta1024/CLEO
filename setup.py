from setuptools import setup, find_packages

setup(
    name='CLEO 4.0',
    version='4.0',
    packages=find_packages(),
    install_requires=[
        'pytest',
        'sphinx',
        'matplotlib',
    ],
)
