from setuptools import setup, find_packages
import os

# Read README
readme_path = "README.md"
if os.path.exists(readme_path):
    with open(readme_path, "r", encoding="utf-8") as f:
        long_description = f.read()
else:
    long_description = "Auto_Mart3: Automated CG Mapping and Parameterization"

# Read requirements
requirements_path = "requirements.txt"
if os.path.exists(requirements_path):
    with open(requirements_path, "r", encoding="utf-8") as f:
        requirements = [line.strip() for line in f if line.strip() and not line.startswith("#")]
else:
    requirements = [
        "numpy>=1.20.0",
        "scipy>=1.7.0",
        "matplotlib>=3.4.0",
        "pandas>=1.3.0",
        "mdtraj>=1.9.0",
        "MDAnalysis>=2.0.0",
    ]

setup(
    name="auto_mart3",
    version="1.0.0",
    author="Anacleto Silva de Souza",
    author_email="anacletosilvadesouza@usp.br",
    description="Automated CG Mapping and Parameterization Pipeline for Martini 3",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/anacletosouza/auto_mart3.git",
    license="MIT",

    packages=find_packages(exclude=["tests", "docs"]),

    include_package_data=True,

    package_data={
        "auto_mart3": [
            "data/**/*",   
            "bin/*.sh",
        ],
    },

    install_requires=requirements,

    entry_points={
        "console_scripts": [
            "auto_mart3=auto_mart3.cli:main",
            "auto-map=auto_mart3.scripts.map_aa_to_cg:main",
            "auto-gen-top=auto_mart3.scripts.generate_cg_top:main",
            "auto-analyze=auto_mart3.scripts.generate_bonds_angles_dihedrals:main",
            "auto-distributions=auto_mart3.scripts.bp_distributions:main",
            "auto-adapt-itp=auto_mart3.scripts.adaptation_gro_itp:main",
            "auto-prep=auto_mart3.scripts.bp_prep:main",
            "auto-plot-distributions=auto_mart3.scripts.plot_distributions:main",
            "bayes-potential-adjust=auto_mart3.scripts.potential_adjustment:main",
        ],
    },

    python_requires=">=3.8",
    zip_safe=False,
)
