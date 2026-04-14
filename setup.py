from setuptools import setup, find_packages
import os
import sys

# Read README
readme_path = "README.md"
if os.path.exists(readme_path):
    with open(readme_path, "r", encoding="utf-8") as f:
        long_description = f.read()
else:
    long_description = "Tools for coarse-grained molecular dynamics and Bayesian potentials"

# Read requirements
requirements_path = "requirements.txt"
if os.path.exists(requirements_path):
    with open(requirements_path, "r", encoding="utf-8") as f:
        requirements = [line.strip() for line in f if line.strip() and not line.startswith("#")]
else:
    # Default requirements if requirements.txt doesn't exist
    requirements = [
        "numpy>=1.20.0",
        "scipy>=1.7.0",
        "matplotlib>=3.4.0",
        "pandas>=1.3.0",
        "mdtraj>=1.9.0",
        "MDAnalysis>=2.0.0",
    ]

# Function to find all data files
def find_package_data():
    """Find all data files to include in the package"""
    package_data = {
        "bayesian_potentials": [
            "bin/*.sh",
            "scripts/*.py",
            "data/**/*",
            "examples/**/*",
        ],
    }
    
    # Include specific data files
    data_dirs = [
        "data/ff_files",
        "data/mdp",
        "data/definitions",
        "examples/example_1_map_aa_to_cg",
        "examples/example_1_map_aa_to_cg/ff_files",
        "examples/example_1_map_aa_to_cg/mdp",
        "examples/example_1_map_aa_to_cg/ndx",
        "examples/example_1_map_aa_to_cg/setup",
    ]
    
    for data_dir in data_dirs:
        if os.path.exists(data_dir):
            package_data["bayesian_potentials"].append(f"{data_dir}/**/*")
    
    return package_data

setup(
    name="bayesian_potentials",
    version="0.1.0",
    author="Anacleto",
    author_email="anacleto@example.com",
    description="Tools for coarse-grained molecular dynamics and Bayesian potentials",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/anacleto/bayesian_potentials",
    license="MIT",
    packages=find_packages(where=".", exclude=["tests", "tests.*", "docs", "examples"]),
    include_package_data=True,
    package_data=find_package_data(),
    install_requires=requirements,
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov>=2.0",
            "black>=21.0",
            "flake8>=3.9",
            "mypy>=0.910",
        ],
        "docs": [
            "sphinx>=4.0",
            "sphinx-rtd-theme>=1.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "bayesian-potentials=bayesian_potentials.cli:main",
            "bp-map=bayesian_potentials.scripts.map_aa_to_cg:main",
            "bp-gen-top=bayesian_potentials.scripts.generate_cg_top:main",
            "bp-analyze=bayesian_potentials.scripts.generate_bonds_angles_dihedrals:main",
            "bp-distributions=bayesian_potentials.scripts.bp_distributions:main",
        ],
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS :: MacOS X",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.8",
    zip_safe=False,
    # Platform specific
    platforms=["Linux", "MacOS-X"],
    # Project URLs
    project_urls={
        "Bug Reports": "https://github.com/anacleto/bayesian_potentials/issues",
        "Source": "https://github.com/anacleto/bayesian_potentials",
        "Documentation": "https://github.com/anacleto/bayesian_potentials/README.md",
    },
)
