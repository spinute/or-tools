************************************************************************
file with basedata            : md272_.bas
initial value random generator: 924074084
************************************************************************
projects                      :  1
jobs (incl. supersource/sink ):  20
horizon                       :  146
RESOURCES
  - renewable                 :  2   R
  - nonrenewable              :  2   N
  - doubly constrained        :  0   D
************************************************************************
PROJECT INFORMATION:
pronr.  #jobs rel.date duedate tardcost  MPM-Time
    1     18      0       27       17       27
************************************************************************
PRECEDENCE RELATIONS:
jobnr.    #modes  #successors   successors
   1        1          3           2   3   4
   2        3          2           5  13
   3        3          3           5  10  11
   4        3          3           5   6   7
   5        3          2          16  18
   6        3          3           8  14  19
   7        3          2          10  12
   8        3          3           9  11  16
   9        3          1          15
  10        3          3          13  15  16
  11        3          2          12  17
  12        3          2          13  15
  13        3          1          18
  14        3          1          17
  15        3          1          18
  16        3          1          17
  17        3          1          20
  18        3          1          20
  19        3          1          20
  20        1          0        
************************************************************************
REQUESTS/DURATIONS:
jobnr. mode duration  R 1  R 2  N 1  N 2
------------------------------------------------------------------------
  1      1     0       0    0    0    0
  2      1     5       6    9    8    0
         2     7       3    8    4    0
         3    10       2    8    2    0
  3      1     1       4    6    0   10
         2     3       3    5    4    0
         3     4       1    2    0    8
  4      1     4       7    7    0    5
         2     4       7    6    4    0
         3    10       2    3    0    5
  5      1     1       3    6    0    9
         2     5       1    6    0    8
         3     5       2    6    0    6
  6      1     3      10    7    3    0
         2     5       9    6    3    0
         3     6       9    3    0    2
  7      1     1       5    8    6    0
         2     3       4    6    5    0
         3     7       4    4    0    5
  8      1     1       8    3    6    0
         2     4       8    3    0    6
         3     8       8    3    3    0
  9      1     2       9    9    9    0
         2     9       7    8    9    0
         3    10       6    8    9    0
 10      1     6       7    9    6    0
         2     8       7    7    0    3
         3    10       7    7    2    0
 11      1     2       4    7    5    0
         2     2       6    5    6    0
         3     8       4    4    0    2
 12      1     5       7    7    0    6
         2     8       6    7    8    0
         3    10       6    7    0    3
 13      1     1       9    7    0    8
         2     2       9    6    4    0
         3     7       8    6    2    0
 14      1     1      10    8   10    0
         2     8       6    6    6    0
         3     8       7    6    0    6
 15      1     6       9    7    0    7
         2     7       5    4    0    4
         3     9       2    2    6    0
 16      1     7       9    3    0    2
         2     9       7    2    1    0
         3     9       6    3    0    1
 17      1     2       9    6    0    2
         2     2      10    6    4    0
         3     8       8    6    0    2
 18      1     6       7   10    6    0
         2     7       6    8    4    0
         3     8       5    8    2    0
 19      1     3       6    8    0    3
         2     4       6    7    0    3
         3     9       6    6    0    2
 20      1     0       0    0    0    0
************************************************************************
RESOURCEAVAILABILITIES:
  R 1  R 2  N 1  N 2
   31   33   52   42
************************************************************************
