from setuptools import setup, find_packages

setup(
    name="scTDRP",
    version="0.1.0",
    author="",
    author_email="",
    description="single-cell Terminal Differentiation Repair Pathway based on Optimal Transport",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.20.0",
        "scipy>=1.7.0",
        "scanpy>=1.9.0",
        "anndata>=0.8.0",
        "POT>=0.9.0",
        "scikit-learn>=1.0.0",
        "matplotlib>=3.5.0",
        "seaborn>=0.11.0",
        "pandas>=1.3.0",
    ],
)
