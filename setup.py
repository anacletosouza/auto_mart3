from setuptools import setup, find_packages
import os

# Read README
with open("README.md", "r", encoding="utf-8") as f:
    long_description = f.read()

# Read requirements
with open("requirements.txt", "r", encoding="utf-8") as f:
    requirements = [line.strip() for line in f if line.strip() and not line.startswith("#")]

setup(
    name="bayesian_potentials",
    version="0.1.0",
    author="Anacleto",
    description="Tools for coarse-grained molecular dynamics and Bayesian potentials",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/anacleto/bayesian_potentials",
    packages=find_packages(),
    include_package_data=True,
    package_data={
        "bayesian_potentials": [
            "bin/*.sh",
            "scripts/*.py",
        ],
    },
    install_requires=requirements,
    entry_points={
        "console_scripts": [
            "bayesian-potentials=bayesian_potentials.cli:main",
            "bp-map=bayesian_potentials.scripts.map_aa_to_cg:main",
            "bp-gen-top=bayesian_potentials.scripts.generate_cg_top:main",
        ],
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    python_requires=">=3.8",
    zip_safe=False,
)
