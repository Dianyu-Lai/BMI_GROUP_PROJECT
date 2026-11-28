"""Run the single-family generator directly from a repository checkout."""

from pathlib import Path
import sys


CODE = Path(__file__).resolve().parent.parent
sys.path[:0] = [str(CODE / "rvne" / "src"), str(CODE)]

try:
    from data_generation.rvne_simulation import main
except ModuleNotFoundError as error:
    if error.name != "numpy":
        raise
    requirements = CODE / "rvne" / "requirements.txt"
    raise SystemExit(f'Install NumPy with `python -m pip install -r "{requirements}"`.') from None


if __name__ == "__main__":
    raise SystemExit(main())
