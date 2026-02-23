from setuptools import setup, find_packages

setup(
    name='ddi_prediction',
    version='0.1.0',
    packages=find_packages(),       # Finds all folders with __init__.py and includes them as packages

    install_requires=[              # Dependencies required to run the code
        "pandas",
        "numpy",
        "scikit-learn",
        "matplotlib",
        "seaborn"
    ],
)