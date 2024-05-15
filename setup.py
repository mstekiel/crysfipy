from setuptools import setup

def get_version():
    import crysfipy
    return crysfipy.__version__

setup(name='crysfipy',
    version=get_version(),
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
        "doc": ['sphinx', 'furo'],
    },
    include_package_data=True,
    zip_safe=False)