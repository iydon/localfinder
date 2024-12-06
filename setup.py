# File: setup.py
from setuptools import setup, find_packages

setup(
    name='localfinder',
    version='0.1.14',
    author='Pengfei Yin',
    author_email='12133074@mail.sustech.edu.cn',
    description='A tool calculating local correlation and enrichment significance of two tracks and finding significantly different genomic regions',
    long_description=open('README.md', encoding='utf-8').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/astudentfromsustech/localfinder',
    packages=find_packages(),
    include_package_data=True,
    classifiers=[
        'Development Status :: 4 - Beta',
        'Programming Language :: Python :: 3',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.8',
    install_requires=[
        'numpy',
        'pandas',
        'scipy',
        'matplotlib',
        'plotly',
        'pyGenomeTracks',
        'scikit-learn',
        'argcomplete'
    ],
    entry_points={
        'console_scripts': [
            'localfinder=localfinder.__main__:main',
        ],
    },
)
