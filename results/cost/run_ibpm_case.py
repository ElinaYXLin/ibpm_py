import sys, types, pathlib

repo_root = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(repo_root))

# Some Python environments have an unrelated PyPI package literally named
# 'py' (an old pytest dependency) installed in site-packages, which shadows
# this repo's py/ namespace package when using `python -m py.ibpm`. Forcibly
# register this repo's py/ directory as sys.modules['py'] before importing,
# bypassing that shadow -- see results/vortall/README.md for the original
# diagnosis of this environment quirk.
pkg = types.ModuleType("py")
pkg.__path__ = [str(repo_root / "py")]
sys.modules["py"] = pkg

from py.ibpm import main

if __name__ == "__main__":
    sys.exit(main(["py.ibpm"] + sys.argv[1:]))
