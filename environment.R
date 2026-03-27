# 1. In your project root, initialise renv (once)
renv::init()

# 2. After installing/using packages normally, snapshot
renv::snapshot()   # writes renv.lock

# 3. On a new machine / HPC node, restore exactly
renv::restore()    # installs exact versions from renv.lock