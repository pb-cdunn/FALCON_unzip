__version__ = '1.1.8'

try:
    import sys, pkg_resources
    sys.stderr.write('falcon-unzip {} (pip thinks "{}")\n'.format(__version__, pkg_resources.get_distribution('falcon-unzip')))
except Exception:
    pass
