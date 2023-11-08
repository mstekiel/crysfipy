from setuptools import setup

setup(name='crysfipy',
    version='1.0',
    description='Tools for analysis of physical properties arising from crystal-field Hamiltonian.',
    url='https://github.com/mstekiel/crysfipy',
    author=['Michał Stękiel', 'Petr Čermák'],
    author_email='michal.stekiel@gmail.com',
    packages=['crysfipy'],
    python_requires='==3.9.*',
    install_requires=[
        'numpy',
        'scipy'
    ],
    extras_require={
        "dev": ['sphinx'],
    },
    include_package_data=True,
    zip_safe=False)