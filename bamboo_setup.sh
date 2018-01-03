type module >& /dev/null || . /mnt/software/Modules/current/init/bash

module load gcc
module load ccache
module load python/2-UCS4
which python
which pip # assumed

export PYTHONUSERBASE=$(pwd)/LOCAL
mkdir -p LOCAL/

export PATH=$PYTHONUSERBASE/bin:$PATH
