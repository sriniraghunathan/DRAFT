from setuptools import find_packages, setup

setup(
    name="FlatSkyLab",
    version="1.0",
    description="Tools for analysing flatsky maps",
    author="Srini Raghunathan",
    author_email="sriniraghuna@gmail.com",
    zip_safe=True,
    packages=find_packages(),
    python_requires=">=3.6",
    install_requires=[
        "pyparsing>=2.0.2",
    ],
    package_data={
        f"{lkl}": ["*"]
        for lkl in ["FlatSkyLab"]
    },
)

