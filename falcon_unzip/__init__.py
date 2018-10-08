__version__ = '1.1.4'

try:
    import sys, pkg_resources
    sys.stderr.write('{}\n'.format(pkg_resources.get_distribution('falcon-unzip')))
except Exception:
    pass
