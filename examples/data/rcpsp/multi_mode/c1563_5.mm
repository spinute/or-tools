************************************************************************
file with basedata            : c1563_.bas
initial value random generator: 273480164
************************************************************************
projects                      :  1
jobs (incl. supersource/sink ):  18
horizon                       :  135
RESOURCES
  - renewable                 :  2   R
  - nonrenewable              :  2   N
  - doubly constrained        :  0   D
************************************************************************
PROJECT INFORMATION:
pronr.  #jobs rel.date duedate tardcost  MPM-Time
    1     16      0       18        6       18
************************************************************************
PRECEDENCE RELATIONS:
jobnr.    #modes  #successors   successors
   1        1          3           2   3   4
   2        3          3           5   7   9
   3        3          2          10  13
   4        3          1           8
   5        3          3           6  11  12
   6        3          1          10
   7        3          1          17
   8        3          1          14
   9        3          3          12  14  15
  10        3          1          14
  11        3          1          16
  12        3          2          16  17
  13        3          1          16
  14        3          1          17
  15        3          1          18
  16        3          1          18
  17        3          1          18
  18        1          0        
************************************************************************
REQUESTS/DURATIONS:
jobnr. mode duration  R 1  R 2  N 1  N 2
------------------------------------------------------------------------
  1      1     0       0    0    0    0
  2      1     2       3   10    5    6
         2     5       3    5    5    5
         3     5       2    6    4    4
  3      1     2       9    6    7    8
         2     5       8    4    6    7
         3     7       7    4    3    7
  4      1     3       3    7    2    8
         2     7       2    5    1    5
         3     8       2    3    1    4
  5      1     1       4    7    3    7
         2     1       5    7    3    5
         3     8       4    5    2    3
  6      1     2       9   10    4    7
         2     6       8   10    4    5
         3    10       8   10    2    5
  7      1     1       4    3   10    7
         2     4       3    2    6    5
         3     8       2    1    4    2
  8      1     8       6    9    6   10
         2     9       5    9    6    5
         3    10       5    8    2    2
  9      1     3       6    8    8    8
         2     8       3    7    8    4
         3    10       3    5    8    4
 10      1     5       5    6    4    9
         2     6       4    4    2    8
         3    10       4    3    2    6
 11      1     2       9    9    9    6
         2     4       8    9    7    3
         3     7       6    9    7    3
 12      1     5       7    7    6    9
         2     7       6    6    5    8
         3     8       5    5    3    8
 13      1     2      10    9    5    7
         2     5       9    5    4    7
         3     8       9    4    3    5
 14      1     5       5    6    8    6
         2     5       6    5    7    6
         3     7       3    3    6    6
 15      1     1      10    8    6    2
         2     3       5    7    5    2
         3    10       4    7    1    2
 16      1     1       6    8    8    6
         2     7       6    7    3    4
         3     9       4    4    1    2
 17      1     2       9    7    9    8
         2     6       8    6    8    4
         3    10       7    4    5    3
 18      1     0       0    0    0    0
************************************************************************
RESOURCEAVAILABILITIES:
  R 1  R 2  N 1  N 2
   32   36  100  114
************************************************************************