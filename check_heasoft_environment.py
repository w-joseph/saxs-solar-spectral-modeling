'''
Checks whether PyXspec is importable in the current environment — i.e.
whether HEASOFT has been initialized in this shell (e.g. via
`source $HEADAS/headas-init.sh`) before Python was launched.

Run this first if you're unsure whether your environment is set up
correctly, before running any script that imports from xspec (coronal_regions.py
and everything that depends on it — see the pipeline diagram in README.md).

Usage: python3 check_heasoft_environment.py
'''
import sys

try:
    from xspec import Xset
    # a trivial call, just to confirm the object is actually functional,
    # not just importable
    Xset.chatter
    print("PyXspec imported successfully — HEASOFT appears to be active in this environment.")
    print("You should be able to run the fitting/plotting scripts from here.")
except ImportError:
    print("Could not import PyXspec.")
    print("This usually means HEASOFT hasn't been initialized in this shell.")
    print("Try running 'source $HEADAS/headas-init.sh' (or the .csh equivalent")
    print("for your shell) in this terminal before running Python again.")
    sys.exit(1)
