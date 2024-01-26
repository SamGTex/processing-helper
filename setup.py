from setuptools import setup, find_packages

setup(
    name='processing_helper',
    version='0.1',  # update the version as needed
    packages=find_packages(),
    install_requires=[
        'ic3_labels',
        'click',
        'pyyaml',
        'numpy',
        # Specify any dependencies your package may have
    ],
)