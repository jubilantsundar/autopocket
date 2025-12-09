"""
AutoPocket: Automated Binding Pocket Discovery
"""
from setuptools import setup, find_packages
from pathlib import Path

# Read the contents of README file
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name="autopocket",
    version="1.0.0",
    author="Sundar Jubilant",
    author_email="jubilantsundar@gmail.com",
    description="Automated binding pocket discovery using AI-powered cofolding",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Jubilantsundar/autopocket",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.8",
    install_requires=[
        "pandas>=1.5.0",
        "numpy>=1.23.0",
        "pyyaml>=6.0",
        "gemmi>=0.6.0",
        "scipy>=1.9.0",
        "scikit-learn>=1.2.0",
    ],
    extras_require={
        "viz": [
            "matplotlib>=3.6.0",
            "seaborn>=0.12.0",
        ],
        "dev": [
            "pytest>=7.0.0",
            "pytest-cov>=4.0.0",
            "black>=22.0.0",
            "flake8>=5.0.0",
            "mypy>=0.990",
        ],
    },
    entry_points={
        "console_scripts": [
            "autopocket=autopocket:main",
        ],
    },
    include_package_data=True,
    keywords="drug-discovery protein-ligand binding-pocket structure-prediction AI SAR",
    project_urls={
        "Bug Reports": "https://github.com/Jubilantsundar/autopocket/issues",
        "Source": "https://github.com/Jubilantsundar/autopocket",
        "Documentation": "https://github.com/Jubilantsundar/autopocket#readme",
    },
)
