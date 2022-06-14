import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="bed2genes",
    version="0.1.2",
    author="Alexander Lachmann",
    author_email="alexander.lachmann@mssm.edu",
    description="Package for turning bed files into gene lists.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/maayanlab/bed2genes",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ],
    package_data={
        "bed2genes": ["data/*"]
    },
    include_package_data=True,
    install_requires=[
        'pandas>=1.1.5',
        'tqdm',
        'numpy'
    ],
    python_requires='>=3.6',
)