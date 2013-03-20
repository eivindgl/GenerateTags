from setuptools import setup

setup(
    name='GenerateTags',
    version='0.1.0',
    author='Eivind Gard Lund',
    author_email='gardlund@gmail.com',
    url='http://github.com/eivindgl/GenerateTags',
    packages=['generatetags'],
    scripts=['bin/generate_tags.py'],
    description='Generate 4C tag sequences from ',
    long_description=open('README.txt').read(),
    install_requires=[
        'Logbook',
        'bx-python',
        'pybedtools',
        'biopython'],
    test_suite='nose.collector',
    tests_require=['nose']
)
