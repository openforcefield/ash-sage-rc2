#!/bin/bash

# rsync -arv -P hpc3:/dfs9/dmobley-lab/lilyw7/ash-sage-rc2/04_benchmark/qm/rmsd-tfd .
# rsync -arv -P hpc3:/dfs9/dmobley-lab/lilyw7/ash-sage-rc2/04_benchmark/qm/all-to-all-rmsd .

# rsync -arv -P hpc3:/dfs9/dmobley-lab/lilyw7/ash-sage-rc2/04_benchmark/qm/images .
# rsync -arv -P hpc3:/dfs9/dmobley-lab/lilyw7/ash-sage-rc2/04_benchmark/qm/*.py .
# rsync -arv -P hpc3:/dfs9/dmobley-lab/lilyw7/ash-sage-rc2/04_benchmark/qm/*.sh .

rsync -rl *.py *.sh hpc3:/dfs9/dmobley-lab/lilyw7/ash-sage-rc2/04_benchmark/qm/