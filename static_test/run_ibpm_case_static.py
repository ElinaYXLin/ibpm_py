import sys, types, pathlib
repo_root = pathlib.Path('/Users/elina/Desktop/SURF2026/ibpm_py-main')
sys.path.insert(0, str(repo_root))
pkg = types.ModuleType('py_static')
pkg.__path__ = [str(repo_root / 'py_static')]
sys.modules['py_static'] = pkg
from py_static.ibpm import main
sys.exit(main(['py_static.ibpm'] + sys.argv[1:]))
