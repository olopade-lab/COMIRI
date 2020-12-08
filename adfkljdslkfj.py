from setuptools import setup, find_packages

setup(
    name='comiri',
    version=0.0.1,
    description='A meta-caller for fusion transcript detection',
    url='https://github.com/olopade-lab/COMIRI',
    author='The Olopade Lab',
    author_email='jbreynier@uchicago.edu',
    license='MIT',
    include_package_data=True,
    packages=find_packages(),
    install_requires=['parsl', 'seaborn', 'joblib', 'sklearn', 'tables', 'pyfaidx', 'upsetplot', 'twine'],
    keywords=['Workflows', 'Scientific computing', 'fusion transcript detection', 'bioinformatics'],
    entry_points = {
        'console_scripts': ['polyfuse=polyfuse.polyfuse:run'],
    }
)