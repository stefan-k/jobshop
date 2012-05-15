# Documentation #

@stefan-k Let's write all documentation in [Markdown](http://warpedvisions.org/projects/markdown-cheat-sheet/)?


## Contents ##

* The files in the **open_jobshop** directory consider the Open Shop Scheduling Problem (OSSP) and use the evolib.

 
## TODO ##

* Maybe use alias Int/Uint instead of Int64?
* Test OSSP results against benchmark db
* OSSP efficiency

## Results ##

### First attempt to plug OSSP into genetic algorithm ###

5 jobs, 9 machines, initial makespan = total time = 482 (time if all jobs are exeuted one after another)

Makespan of permutation [1,2,3,...,50]: 169 (reduced to 35%)

| Popsize | Max generations| best found makespan | reduced to |   @time |
| ------: | -------------: | ------------------: | ---------: | ------- |
|      10 |              1 |                 153 |     31.74% |         |
|      10 |             10 |                 151 |     31.33% |         |
|      10 |            100 |                 164 |     34.02% |         |
|     100 |              1 |                 150 |     31.12% |         |
|     100 |             10 |                 148 |     30.71% |         |
|     100 |            100 |                 146 |     30.29% |         |
|     100 |           1000 |                 148 |     30.71% |   83.29s|


Evaluation of 10 random problems with different genetic probabilities:

