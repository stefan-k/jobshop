@stefan-k Let's write all documentation in [Markdown](http://warpedvisions.org/projects/markdown-cheat-sheet/)?


# Contents #

* The files in the **open_jobshop** directory consider the Open Shop Scheduling Problem (OSSP) and use the evolib.


# TODO #

[ ] Maybe use alias Int/Uint instead of Int64?
[ ] Test OSSP results against benchmark db
[ ] OSSP efficiency

# Results #

First attempt to plug OSSP into genetic algorithm:

5 jobs, 9 machines, initial makespan = total time = 482 (time if all jobs are exeuted one after another)

Makespan of permutation [1,2,3,...,50]: 169 (reduced to 35%)

Popsize= 10, Max generations=   1, optimal makespan: 153, reduced to: 31.74%
Popsize= 10, Max generations=  10, optimal makespan: 151, reduced to: 31.33%
Popsize= 10, Max generations= 100, optimal makespan: 164, reduced to: 34.02%

Popsize=100, Max generations=   1, optimal makespan: 150, reduced to: 31.12%
Popsize=100, Max generations=  10, optimal makespan: 148, reduced to: 30.71%
Popsize=100, Max generations= 100, optimal makespan: 146, reduced to: 30.29%

Popsize=100, Max generations=1000, optimal makespan: 148, reduced to: 30.71% (elapsed time: 83.28931593894958 seconds)