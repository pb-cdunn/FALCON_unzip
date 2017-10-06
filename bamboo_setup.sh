type module >& /dev/null || . /mnt/software/Modules/current/init/bash

module load gcc/6.4.0
module load python/2.7.13-UCS2
which pip

export PYTHONUSERBASE=$(pwd)/LOCAL
mkdir -p LOCAL/

export PATH=$PYTHONUSERBASE/bin:$PATH
