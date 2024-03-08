from setuptools import setup

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

pkg = 'reliontomotools'
exec(open(f'{pkg}/_version.py').read())

setup(
      name=pkg,
      version=__version__,
      author="Joaquín Otón",
      description="Additional tools for Relion tomo",
      long_description=long_description,
      long_description_content_type="text/markdown",
      url="https://github.com/joton/reliontomotools",
      project_urls={
        "Bug Tracker": "https://github.com/joton/reliontomotools/issues",
                    },
      classifiers=[
          "Programming Language :: Python :: 3",
          "License :: OSI Approved :: MIT License",
          "Operating System :: OS Independent",
                  ],
      packages=[pkg],
      license='MIT',
      entry_points={
          "console_scripts": [
              f"warptomo2relion = {pkg}.warptomo2relion:warpTomo2RelionProgram",
              f"motl2relion = {pkg}.motl2relion:motl2RelionProgram",
              f"relion2motl = {pkg}.relion2motl:relion2motlProgram",
          ]
                   },
      python_requires=">=3.7",
      install_requires=[
          'numpy',
          'pandas',
          'xmltodict',
          'docopt',
          'scipy',
          'tqdm',
          'transforms3d',
          'mrcfile',
          'emfile',
          ]

     )
