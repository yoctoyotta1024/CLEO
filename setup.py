from setuptools import setup, find_packages

setup(
    name='version3.0',
    version='3.0',
    packages=find_packages(),
    install_requires=[
        'pytest',
        'sphinx',
        'matplotlib',
    ],
)
