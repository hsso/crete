Installing the code
-------

From the top-level directory, run the `configure.sh` shell script to install
`rat4com`.

```
./configure.sh
```

This script is a simplified version of the original `ratran`'s `configure`
script that works for `bash` and `zsh` shells on Linux, and probably also
MacOS, using the `gfortran` compiler.

If you are using the C-shell or the improved version tcsh for interactive use,
you can instead run the csh script to install the program:

```
./configure
```

After running either of these scripts you can add the definition of the RATRAN
and RATRANRUN environment variables to initialize these variables for login
shells in the shell-specific files `~/.bashrc`, `~/.zshrc` or `~/.cshrc`.


The IR radiation is specified as pumping coefficients in the file
`amc/getmatrix.f`.  One needs to use the `g_ir` array for ortho- or
para-water by uncommenting the corresponding definition and include
the heliocentric distance to scale the coefficients in lines 103 and 104.
