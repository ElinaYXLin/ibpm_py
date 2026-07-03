import sys, types, pathlib
repo_root = pathlib.Path('/Users/elina/Desktop/SURF2026/ibpm_py-main')
sys.path.insert(0, str(repo_root))
pkg = types.ModuleType('py')
pkg.__path__ = [str(repo_root / 'py')]
sys.modules['py'] = pkg
from py.ibpm import main
sys.exit(main(['py.ibpm'] + sys.argv[1:]))
