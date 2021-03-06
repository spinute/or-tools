************************************************************************
file with basedata            : md298_.bas
initial value random generator: 382594544
************************************************************************
projects                      :  1
jobs (incl. supersource/sink ):  20
horizon                       :  152
RESOURCES
  - renewable                 :  2   R
  - nonrenewable              :  2   N
  - doubly constrained        :  0   D
************************************************************************
PROJECT INFORMATION:
pronr.  #jobs rel.date duedate tardcost  MPM-Time
    1     18      0       18       12       18
************************************************************************
PRECEDENCE RELATIONS:
jobnr.    #modes  #successors   successors
   1        1          3           2   3   4
   2        3          2           8  18
   3        3          2           5   8
   4        3          3           6  10  15
   5        3          3           6   9  11
   6        3          3           7  16  19
   7        3          1          13
   8        3          2          11  13
   9        3          3          10  12  13
  10        3          3          17  18  19
  11        3          1          15
  12        3          2          14  17
  13        3          1          17
  14        3          2          15  16
  15        3          1          19
  16        3          1          18
  17        3          1          20
  18        3          1          20
  19        3          1          20
  20        1          0        
************************************************************************
REQUESTS/DURATIONS:
jobnr. mode duration  R 1  R 2  N 1  N 2
------------------------------------------------------------------------
  1      1     0       0    0    0    0
  2      1     4       7    0    2    9
         2     7       5    0    2    6
         3     8       3    0    1    5
  3      1     1       5    0    9    6
         2     7       0    8    7    3
         3     9       0    6    6    1
  4      1     1       9    0    7    8
         2     6       0    4    5    6
         3    10       8    0    4    5
  5      1     1       0    6    8    8
         2     4       5    0    4    8
         3     8       0    3    2    7
  6      1     4       0    9    7   10
         2     5       1    0    6    9
         3     7       0    6    4    9
  7      1     2       0    2    7    8
         2     3       9    0    4    6
         3     8       0    1    1    5
  8      1     5       0    4    7    2
         2     5       0    2    8    2
         3    10       6    0    5    2
  9      1     1       5    0    9    9
         2     4       4    0    9    9
         3    10       0    1    9    8
 10      1     2       8    0    8    5
         2     8       8    0    5    2
         3     8       0    2    2    3
 11      1     1       0    8    6    9
         2     2       0    5    6    8
         3     3       6    0    6    8
 12      1     4       0    9    6    1
         2     5       0    6    5    1
         3     8       0    4    3    1
 13      1     1       7    0    6    7
         2     3       6    0    4    5
         3     8       0    8    4    5
 14      1     4       8    0    6    4
         2     5       6    0    5    3
         3     7       0    8    5    2
 15      1     3       0    4    5    8
         2     4       9    0    5    7
         3     9       9    0    4    5
 16      1     6       8    0    4    8
         2     7       8    0    4    7
         3    10       7    0    2    6
 17      1     1       0    5    4    5
         2     6       8    0    3    5
         3     9       0    5    2    4
 18      1     1       8    0    6    4
         2     6       0    3    6    4
         3    10       7    0    4    4
 19      1     4       0    2    9    6
         2     7       0    2    9    4
         3    10       0    2    9    1
 20      1     0       0    0    0    0
************************************************************************
RESOURCEAVAILABILITIES:
  R 1  R 2  N 1  N 2
   16   12   95   99
************************************************************************
