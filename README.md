This folder is a loose aggregation of optimization model code and plotting
scripts related to the joint mapping and routing project from 2010--11.

The programs require the mosek optimization library (https://www.mosek.com/) to
successfully compile and execute. A Linux-like environment (with gcc and python)
will also likely be necessary.

The best place to start understanding the code is from the main function in
optimality.c. The helper files compile-commands.txt and plotting-commands.txt
contain useful information to compile and execute the software.

This software is offered as is, and is not guaranteed to be fit for any given
purpose. It is not supported or maintained. It is made available simply to allow
people to use it as a starting point to develop their own optimization models.
