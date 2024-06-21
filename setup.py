from cx_Freeze import setup, Executable

executables = [Executable('main.py',
                              target_name='Protein.exe',
                              base='Win32GUI')]

options = {"build_exe":
           {
               "packages": ["flask",  "flask_cors", "PyQt5",],
               'includes': ['PyQt5'],
               "include_files": ["server/"],
           }
           }


setup(name='Protein3D',
      version='0.0.1',
      description='3D Molecular Visualizer',
      executables=executables,
      options=options)
