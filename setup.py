from setuptools import setup

setup(name='diva_seq_opt',
      version='0.0.2',
      description='Predict the library loading concentration to use for mi-seq runs based on bioanalyzer data.',
      url='https://github.com/JBEI/library_loading_concentration_prediction',
      author='Zak Costello',
      author_email='zak.costello@gmail.com',
      license='BSD',
      packages=['diva_seq_opt'],      
      entry_points = {
        'console_scripts': ['predict_loading_concentration=diva_seq_opt.predict_loading_concentration:main'],
      },
      install_requires=[
          'pandas','numpy','scipy','matplotlib',
          'peakutils','sklearn','tpot'
      ],
      include_package_data=True,
      zip_safe=False)