#!/usr/bin/env python

from setuptools import setup, find_packages

from distutils.core import Extension

import glob

install_requires = ["networkx >= 1.7", "pysam >= 0.8.4", "msgpack-python"]

#scripts = glob.glob("src/py_scripts/*.py")

setup(name='falcon_unzip',
      version='0.4.0',
      description='Falcon unzip',
      author='Jason Chin',
      author_email='jchin@pacificbiosciences.com',
      maintainer='Christopher Dunn',
      maintainer_email='pb.cdunn@gmail.com',
      packages=find_packages(),
      package_dir={'falcon_unzip': 'falcon_unzip/'},
      entry_points={'console_scripts': [
          # We do not really need most of these, but mobs might depend on them.
          'fc_dedup_h_tigs.py=falcon_unzip.mains.dedup_h_tigs:main',
          'fc_get_read_hctg_map.py=falcon_unzip.mains.get_read_hctg_map:main',
          'fc_graphs_to_h_tigs.py=falcon_unzip.mains.graphs_to_h_tigs:main',
          'fc_ovlp_filter_with_phase.py=falcon_unzip.mains.ovlp_filter_with_phase:main',
          'fc_phased_ovlp_to_graph.py=falcon_unzip.mains.phased_ovlp_to_graph:main',
          'fc_phasing.py=falcon_unzip.mains.phasing:main',
          'fc_phasing_readmap.py=falcon_unzip.mains.phasing_readmap:main',
          'fc_quiver.py=falcon_unzip.mains.start_quiver:main',
          'fc_rr_hctg_track.py=falcon_unzip.mains.rr_hctg_track:main',
          'fc_rr_hctg_track2.py=falcon_unzip.mains.rr_hctg_track:main2',
          'fc_select_reads_from_bam.py=falcon_unzip.mains.select_reads_from_bam:main',  # not used?
          'fc_unzip.py=falcon_unzip.mains.start_unzip:main',
          'fc_unzip_gen_gfa_v1.py=falcon_unzip.mains.unzip_gen_gfa_v1:main',
      ]},
      #scripts = scripts,
      zip_safe=False,
      install_requires=install_requires
      )
