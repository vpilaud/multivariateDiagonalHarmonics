Precalculated_F=Family({
Partition( [0] ):(
+tensor([s([0]),e( [] )])
),
Partition( [1] ):(
+tensor([ s[1] ,e( [] )])
+tensor([s([0]),e( [1] )])
),
Partition( [1, 1] ):(
+tensor([ s[2] ,e( [] )])
+tensor([ s[1] ,e( [1] )])
+tensor([s([0]),e( [2] )])
),
Partition( [2] ):(
+tensor([ s[2] ,e( [] )])
+tensor([ s[0] + s[1] ,e( [1] )])
),
Partition( [1, 1, 1] ):(
+tensor([ s[3] ,e( [] )])
+tensor([ s[2] ,e( [1] )])
+tensor([ s[1] ,e( [2] )])
+tensor([s([0]),e( [3] )])
),
Partition( [2, 1] ):(
+tensor([ s[1, 1] + s[3] ,e( [] )])
+tensor([ s[1] + s[2] ,e( [1] )])
+tensor([s([0]),e( [1, 1] )])
+tensor([ s[1] ,e( [2] )])
),
Partition( [3] ):(
+tensor([ s[3] ,e( [] )])
+tensor([ s[0] + s[1] + s[2] ,e( [1] )])
),
Partition( [1, 1, 1, 1] ):(
+tensor([ s[4] ,e( [] )])
+tensor([ s[3] ,e( [1] )])
+tensor([ s[2] ,e( [2] )])
+tensor([ s[1] ,e( [3] )])
+tensor([s([0]),e( [4] )])
),
Partition( [2, 1, 1] ):(
+tensor([ s[2, 1] + s[4] ,e( [] )])
+tensor([ s[1, 1] + s[2] + s[3] ,e( [1] )])
+tensor([ s[1] ,e( [1, 1] )])
+tensor([ s[2] ,e( [2] )])
+tensor([s([0]),e( [2, 1] )])
+tensor([ s[1] ,e( [3] )])
),
Partition( [3, 1] ):(
+tensor([ s[2, 1] + s[4] ,e( [] )])
+tensor([ s[1] + s[1, 1] + s[2] + s[3] ,e( [1] )])
+tensor([ s[0] + s[1] ,e( [1, 1] )])
+tensor([ s[2] ,e( [2] )])
),
Partition( [4] ):(
+tensor([ s[4] ,e( [] )])
+tensor([ s[0] + s[1] + s[2] + s[3] ,e( [1] )])
),
Partition( [1, 1, 1, 1, 1] ):(
+tensor([ s[5] ,e( [] )])
+tensor([ s[4] ,e( [1] )])
+tensor([ s[3] ,e( [2] )])
+tensor([ s[2] ,e( [3] )])
+tensor([ s[1] ,e( [4] )])
+tensor([s([0]),e( [5] )])
),
Partition( [2, 1, 1, 1] ):(
+tensor([ s[3, 1] + s[5] ,e( [] )])
+tensor([ s[2, 1] + s[3] + s[4] ,e( [1] )])
+tensor([ s[2] ,e( [1, 1] )])
+tensor([ s[1, 1] + s[3] ,e( [2] )])
+tensor([ s[1] ,e( [2, 1] )])
+tensor([ s[2] ,e( [3] )])
+tensor([s([0]),e( [3, 1] )])
+tensor([ s[1] ,e( [4] )])
),
Partition( [2, 2, 1] ):(
+tensor([ s[3, 1] + s[5] ,e( [] )])
+tensor([ s[2, 1] + s[3] + s[4] ,e( [1] )])
+tensor([ s[2] ,e( [1, 1] )])
+tensor([ s[1] + s[1, 1] + s[3] ,e( [2] )])
+tensor([ s[0] + s[1] ,e( [2, 1] )])
+tensor([ s[2] ,e( [3] )])
),
Partition( [3, 2] ):(
+tensor([ s[3, 1] + s[5] ,e( [] )])
+tensor([ s[2] + s[2, 1] + s[3] + s[4] ,e( [1] )])
+tensor([ s[0] + s[1] + s[2] ,e( [1, 1] )])
+tensor([ s[1] + s[1, 1] + s[3] ,e( [2] )])
),
Partition( [3, 1, 1] ):(
+tensor([ s[5] + s[3, 1] + s[2, 2] ,e( [] )])
+tensor([ s[4] + s[3] + s[2] + 2*s[2, 1] ,e( [1] )])
+tensor([ s[1] + s[2] + s[1, 1] ,e( [1, 1] )])
+tensor([ s[3] ,e( [2] )])
+tensor([ s[2] ,e( [3] )])
+tensor([ s[0]+s[1] ,e( [2,1] )])
),
Partition( [4, 1] ):(
+tensor([ s[3, 1] + s[5] ,e( [] )])
+tensor([ s[1] + s[1, 1] + s[2] + s[2, 1] + s[3] + s[4] ,e( [1] )])
+tensor([ s[0] + s[1] + s[2] ,e( [1, 1] )])
+tensor([ s[3] ,e( [2] )])
),
Partition( [5] ):(
+tensor([ s[5] ,e( [] )])
+tensor([ s[0] + s[1] + s[2] + s[3] + s[4] ,e( [1] )])
),
Partition( [1, 1, 1, 1, 1, 1] ):(
+tensor([ s[6] ,e( [] )])
+tensor([ s[5] ,e( [1] )])
+tensor([ s[4] ,e( [2] )])
+tensor([ s[3] ,e( [3] )])
+tensor([ s[2] ,e( [4] )])
+tensor([ s[1] ,e( [5] )])
+tensor([s([0]),e( [6] )])
),
Partition( [2, 1, 1, 1, 1] ):(
+tensor([ s[4, 1] + s[6] ,e( [] )])
+tensor([ s[3, 1] + s[4] + s[5] ,e( [1] )])
+tensor([ s[3] ,e( [1, 1] )])
+tensor([ s[2, 1] + s[4] ,e( [2] )])
+tensor([ s[2] ,e( [2, 1] )])
+tensor([ s[1, 1] + s[3] ,e( [3] )])
+tensor([ s[1] ,e( [3, 1] )])
+tensor([ s[2] ,e( [4] )])
+tensor([s([0]),e( [4, 1] )])
+tensor([ s[1] ,e( [5] )])
),
Partition( [2, 2, 1, 1] ):(
+tensor([ s[2, 2] + s[4, 1] + s[6] ,e( [] )])
+tensor([ s[2, 1] + s[3, 1] + s[4] + s[5] ,e( [1] )])
+tensor([ s[1, 1] + s[3] ,e( [1, 1] )])
+tensor([ -s[1, 1] + s[2] + s[2, 1] + s[4] ,e( [2] )])
+tensor([ s[1] + s[2] ,e( [2, 1] )])
+tensor([s([0]),e( [2, 2] )])
+tensor([ s[1, 1] + s[3] ,e( [3] )])
+tensor([ s[1] ,e( [3, 1] )])
+tensor([ s[2] ,e( [4] )])
),
Partition( [3, 2, 1] ):(
+tensor([ s[1, 1, 1] + s[3, 1] + s[4, 1] + s[6] ,e( [] )])
+tensor([ s[1, 1] + s[2, 1] + s[3] + s[3, 1] + s[4] + s[5] ,e( [1] )])
+tensor([ s[1] + s[2] + s[3] ,e( [1, 1] )])
+tensor([s([0]),e( [1, 1, 1] )])
+tensor([ s[1, 1] + s[2] + s[2, 1] + s[4] ,e( [2] )])
+tensor([ 2*s[1] + s[2] ,e( [2, 1] )])
+tensor([ s[1, 1] + s[3] ,e( [3] )])
),
Partition( [4, 2] ):(
+tensor([ s[2, 2] + s[4, 1] + s[6] ,e( [] )])
+tensor([ s[2] + s[2, 1] + s[3] + s[3, 1] + s[4] + s[5] ,e( [1] )])
+tensor([ s[0] + 2*s[1] + s[1, 1] + s[2] + s[3] ,e( [1, 1] )])
+tensor([ s[2] + s[2, 1] + s[4] ,e( [2] )])
),
Partition( [5, 1] ):(
+tensor([ s[4, 1] + s[6] ,e( [] )])
+tensor([ s[1] + s[1, 1] + s[2] + s[2, 1] + s[3] + s[3, 1] + s[4] + s[5] ,e( [1] )])
+tensor([ s[0] + s[1] + s[2] + s[3] ,e( [1, 1] )])
+tensor([ s[4] ,e( [2] )])
),
Partition( [6] ):(
+tensor([ s[6] ,e( [] )])
+tensor([ s[0] + s[1] + s[2] + s[3] + s[4] + s[5] ,e( [1] )])
),
Partition( [1, 1, 1, 1, 1, 1, 1] ):(
+tensor([ s[7] ,e( [] )])
+tensor([ s[6] ,e( [1] )])
+tensor([ s[5] ,e( [2] )])
+tensor([ s[4] ,e( [3] )])
+tensor([ s[3] ,e( [4] )])
+tensor([ s[2] ,e( [5] )])
+tensor([ s[1] ,e( [6] )])
+tensor([s([0]),e( [7] )])
),
Partition( [2, 1, 1, 1, 1, 1] ):(
+tensor([ s[5, 1] + s[7] ,e( [] )])
+tensor([ s[4, 1] + s[5] + s[6] ,e( [1] )])
+tensor([ s[4] ,e( [1, 1] )])
+tensor([ s[3, 1] + s[5] ,e( [2] )])
+tensor([ s[3] ,e( [2, 1] )])
+tensor([ s[2, 1] + s[4] ,e( [3] )])
+tensor([ s[2] ,e( [3, 1] )])
+tensor([ s[1, 1] + s[3] ,e( [4] )])
+tensor([ s[1] ,e( [4, 1] )])
+tensor([ s[2] ,e( [5] )])
+tensor([s([0]),e( [5, 1] )])
+tensor([ s[1] ,e( [6] )])
),
Partition( [2, 2, 1, 1, 1] ):(
+tensor([ s[3, 2] + s[5, 1] + s[7] ,e( [] )])
+tensor([ s[2, 2] + s[3, 1] + s[4, 1] + s[5] + s[6] ,e( [1] )])
+tensor([ s[2, 1] + s[4] ,e( [1, 1] )])
+tensor([ s[3] + s[3, 1] + s[5] ,e( [2] )])
+tensor([ s[1, 1] + s[2] + s[3] ,e( [2, 1] )])
+tensor([ s[1] ,e( [2, 2] )])
+tensor([ -s[1, 1] + s[2, 1] + s[4] ,e( [3] )])
+tensor([ s[2] ,e( [3, 1] )])
+tensor([s([0]),e( [3, 2] )])
+tensor([ s[1, 1] + s[3] ,e( [4] )])
+tensor([ s[1] ,e( [4, 1] )])
+tensor([ s[2] ,e( [5] )])
),
Partition( [3, 2, 1, 1] ):(
+tensor([ s[2, 1, 1] + s[3, 2] + s[4, 1] + s[5, 1] + s[7] ,e( [] )])
+tensor([ s[1, 1, 1] + s[2, 1] + s[2, 2] + 2*s[3, 1] + s[4] + s[4, 1] + s[5] + s[6] ,e( [1] )])
+tensor([ s[1, 1] + s[2] + s[2, 1] + s[3] + s[4] ,e( [1, 1] )])
+tensor([ s[1] ,e( [1, 1, 1] )])
+tensor([ s[2, 1] + s[3] + s[3, 1] + s[5] ,e( [2] )])
+tensor([ s[1, 1] + 2*s[2] + s[3] ,e( [2, 1] )])
+tensor([s([0]),e( [2, 1, 1] )])
+tensor([ s[1] ,e( [2, 2] )])
+tensor([ s[2, 1] + s[4] ,e( [3] )])
+tensor([ s[1] + s[2] ,e( [3, 1] )])
+tensor([ s[1, 1] + s[3] ,e( [4] )])
),
Partition( [4, 2, 1] ):(
+tensor([ s[2, 1, 1] + s[3, 2] + s[4, 1] + s[5, 1] + s[7] ,e( [] )])
+tensor([ s[1, 1] + s[1, 1, 1] + s[2, 1] + s[2, 2] + s[3] + 2*s[3, 1] + s[4] + s[4, 1] + s[5] + s[6] ,e( [1] )])
+tensor([ s[1] + s[1, 1] + 2*s[2] + s[2, 1] + s[3] + s[4] ,e( [1, 1] )])
+tensor([ s[0] + s[1] ,e( [1, 1, 1] )])
+tensor([ s[2, 1] + s[3] + s[3, 1] + s[5] ,e( [2] )])
+tensor([ s[1] + s[1, 1] + 2*s[2] + s[3] ,e( [2, 1] )])
+tensor([ s[2, 1] + s[4] ,e( [3] )])
),
Partition( [5, 2] ):(
+tensor([ s[3, 2] + s[5, 1] + s[7] ,e( [] )])
+tensor([ s[2] + s[2, 1] + s[2, 2] + s[3] + s[3, 1] + s[4] + s[4, 1] + s[5] + s[6] ,e( [1] )])
+tensor([ s[0] + 2*s[1] + s[1, 1] + 2*s[2] + s[2, 1] + s[3] + s[4] ,e( [1, 1] )])
+tensor([ s[3] + s[3, 1] + s[5] ,e( [2] )])
),
Partition( [6, 1] ):(
+tensor([ s[5, 1] + s[7] ,e( [] )])
+tensor([ s[1] + s[1, 1] + s[2] + s[2, 1] + s[3] + s[3, 1] + s[4] + s[4, 1] + s[5] + s[6] ,e( [1] )])
+tensor([ s[0] + s[1] + s[2] + s[3] + s[4] ,e( [1, 1] )])
+tensor([ s[5] ,e( [2] )])
),
Partition( [7] ):(
+tensor([ s[7] ,e( [] )])
+tensor([ s[0] + s[1] + s[2] + s[3] + s[4] + s[5] + s[6] ,e( [1] )])
),
Partition( [1, 1, 1, 1, 1, 1, 1, 1] ):(
+tensor([ s[8] ,e( [] )])
+tensor([ s[7] ,e( [1] )])
+tensor([ s[6] ,e( [2] )])
+tensor([ s[5] ,e( [3] )])
+tensor([ s[4] ,e( [4] )])
+tensor([ s[3] ,e( [5] )])
+tensor([ s[2] ,e( [6] )])
+tensor([ s[1] ,e( [7] )])
+tensor([s([0]),e( [8] )])
),
Partition( [2, 1, 1, 1, 1, 1, 1] ):(
+tensor([ s[6, 1] + s[8] ,e( [] )])
+tensor([ s[5, 1] + s[6] + s[7] ,e( [1] )])
+tensor([ s[5] ,e( [1, 1] )])
+tensor([ s[4, 1] + s[6] ,e( [2] )])
+tensor([ s[4] ,e( [2, 1] )])
+tensor([ s[3, 1] + s[5] ,e( [3] )])
+tensor([ s[3] ,e( [3, 1] )])
+tensor([ s[2, 1] + s[4] ,e( [4] )])
+tensor([ s[2] ,e( [4, 1] )])
+tensor([ s[1, 1] + s[3] ,e( [5] )])
+tensor([ s[1] ,e( [5, 1] )])
+tensor([ s[2] ,e( [6] )])
+tensor([s([0]),e( [6, 1] )])
+tensor([ s[1] ,e( [7] )])
),
Partition( [2, 2, 1, 1, 1, 1] ):(
+tensor([ s[4, 2] + s[6, 1] + s[8] ,e( [] )])
+tensor([ s[3, 2] + s[4, 1] + s[5, 1] + s[6] + s[7] ,e( [1] )])
+tensor([ s[3, 1] + s[5] ,e( [1, 1] )])
+tensor([ s[2, 2] + s[4] + s[4, 1] + s[6] ,e( [2] )])
+tensor([ s[2, 1] + s[3] + s[4] ,e( [2, 1] )])
+tensor([ s[2] ,e( [2, 2] )])
+tensor([ s[3, 1] + s[5] ,e( [3] )])
+tensor([ s[1, 1] + s[3] ,e( [3, 1] )])
+tensor([ s[1] ,e( [3, 2] )])
+tensor([ -s[1, 1] + s[2, 1] + s[4] ,e( [4] )])
+tensor([ s[2] ,e( [4, 1] )])
+tensor([s([0]),e( [4, 2] )])
+tensor([ s[1, 1] + s[3] ,e( [5] )])
+tensor([ s[1] ,e( [5, 1] )])
+tensor([ s[2] ,e( [6] )])
),
Partition( [2, 2, 2, 1, 1] ):(
+tensor([ s[4, 2] + s[6, 1] + s[8] ,e( [] )])
+tensor([ s[3, 2] + s[4, 1] + s[5, 1] + s[6] + s[7] ,e( [1] )])
+tensor([ s[3, 1] + s[5] ,e( [1, 1] )])
+tensor([ s[2, 2] + s[4] + s[4, 1] + s[6] ,e( [2] )])
+tensor([ s[2, 1] + s[3] + s[4] ,e( [2, 1] )])
+tensor([ s[2] ,e( [2, 2] )])
+tensor([ s[2] + s[3, 1] + s[5] ,e( [3] )])
+tensor([ s[1] + s[1, 1] + s[3] ,e( [3, 1] )])
+tensor([ s[0] + s[1] ,e( [3, 2] )])
+tensor([ s[2, 1] + s[4] ,e( [4] )])
+tensor([ s[2] ,e( [4, 1] )])
+tensor([ s[3] ,e( [5] )])
),
Partition( [3, 2, 2, 1] ):(
+tensor([ s[3, 1, 1] + s[4, 2] + s[5, 1] + s[6, 1] + s[8] ,e( [] )])
+tensor([ s[2, 1, 1] + s[3, 1] + s[3, 2] + 2*s[4, 1] + s[5] + s[5, 1] + s[6] + s[7] ,e( [1] )])
+tensor([ s[2, 1] + s[3] + s[3, 1] + s[4] + s[5] ,e( [1, 1] )])
+tensor([ s[2] ,e( [1, 1, 1] )])
+tensor([ s[1, 1] + s[1, 1, 1] + s[2, 2] + s[3, 1] + s[4] + s[4, 1] + s[6] ,e( [2] )])
+tensor([ s[1] + s[1, 1] + s[2, 1] + 2*s[3] + s[4] ,e( [2, 1] )])
+tensor([ s[0] + s[1] ,e( [2, 1, 1] )])
+tensor([ s[2] ,e( [2, 2] )])
+tensor([ s[2] + s[2, 1] + s[3, 1] + s[5] ,e( [3] )])
+tensor([ s[1] + s[1, 1] + s[2] + s[3] ,e( [3, 1] )])
+tensor([ s[2, 1] + s[4] ,e( [4] )])
),
Partition( [4, 3, 1] ):(
+tensor([ s[3, 1, 1] + s[4, 2] + s[5, 1] + s[6, 1] + s[8] ,e( [] )])
+tensor([ s[2, 1] + s[2, 1, 1] + s[3, 1] + s[3, 2] + s[4] + 2*s[4, 1] + s[5] + s[5, 1] + s[6] + s[7] ,e( [1] )])
+tensor([ s[1] + s[1, 1] + s[2] + s[2, 1] + 2*s[3] + s[3, 1] + s[4] + s[5] ,e( [1, 1] )])
+tensor([ s[0] + s[1] + s[2] ,e( [1, 1, 1] )])
+tensor([ s[1, 1] + s[1, 1, 1] + s[2] + s[2, 1] + s[2, 2] + s[3, 1] + s[4] + s[4, 1] + s[6] ,e( [2] )])
+tensor([ s[1] + s[1, 1] + s[2] + s[2, 1] + 2*s[3] + s[4] ,e( [2, 1] )])
+tensor([ s[3, 1] + s[5] ,e( [3] )])
),
Partition( [5, 3] ):(
+tensor([ s[4, 2] + s[6, 1] + s[8] ,e( [] )])
+tensor([ s[3] + s[3, 1] + s[3, 2] + s[4] + s[4, 1] + s[5] + s[5, 1] + s[6] + s[7] ,e( [1] )])
+tensor([ s[0] + 2*s[1] + s[1, 1] + 2*s[2] + s[2, 1] + 2*s[3] + s[3, 1] + s[4] + s[5] ,e( [1, 1] )])
+tensor([ s[2] + s[2, 1] + s[2, 2] + s[4] + s[4, 1] + s[6] ,e( [2] )])
),
Partition( [6, 2] ):(
+tensor([ s[4, 2] + s[6, 1] + s[8] ,e( [] )])
+tensor([ s[2] + s[2, 1] + s[2, 2] + s[3] + s[3, 1] + s[3, 2] + s[4] + s[4, 1] + s[5] + s[5, 1] + s[6] + s[7] ,e( [1] )])
+tensor([ s[0] + 2*s[1] + s[1, 1] + 2*s[2] + s[2, 1] + 2*s[3] + s[3, 1] + s[4] + s[5] ,e( [1, 1] )])
+tensor([ s[4] + s[4, 1] + s[6] ,e( [2] )])
),
Partition( [7, 1] ):(
+tensor([ s[6, 1] + s[8] ,e( [] )])
+tensor([ s[1] + s[1, 1] + s[2] + s[2, 1] + s[3] + s[3, 1] + s[4] + s[4, 1] + s[5] + s[5, 1] + s[6] + s[7] ,e( [1] )])
+tensor([ s[0] + s[1] + s[2] + s[3] + s[4] + s[5] ,e( [1, 1] )])
+tensor([ s[6] ,e( [2] )])
),
Partition( [8] ):(
+tensor([ s[8] ,e( [] )])
+tensor([ s[0] + s[1] + s[2] + s[3] + s[4] + s[5] + s[6] + s[7] ,e( [1] )])
),
Partition( [1, 1, 1, 1, 1, 1, 1, 1, 1] ):(
+tensor([ s[9] ,e( [] )])
+tensor([ s[8] ,e( [1] )])
+tensor([ s[7] ,e( [2] )])
+tensor([ s[6] ,e( [3] )])
+tensor([ s[5] ,e( [4] )])
+tensor([ s[4] ,e( [5] )])
+tensor([ s[3] ,e( [6] )])
+tensor([ s[2] ,e( [7] )])
+tensor([ s[1] ,e( [8] )])
+tensor([s([0]),e( [9] )])
),
Partition( [2, 1, 1, 1, 1, 1, 1, 1] ):(
+tensor([ s[7, 1] + s[9] ,e( [] )])
+tensor([ s[6, 1] + s[7] + s[8] ,e( [1] )])
+tensor([ s[6] ,e( [1, 1] )])
+tensor([ s[5, 1] + s[7] ,e( [2] )])
+tensor([ s[5] ,e( [2, 1] )])
+tensor([ s[4, 1] + s[6] ,e( [3] )])
+tensor([ s[4] ,e( [3, 1] )])
+tensor([ s[3, 1] + s[5] ,e( [4] )])
+tensor([ s[3] ,e( [4, 1] )])
+tensor([ s[2, 1] + s[4] ,e( [5] )])
+tensor([ s[2] ,e( [5, 1] )])
+tensor([ s[1, 1] + s[3] ,e( [6] )])
+tensor([ s[1] ,e( [6, 1] )])
+tensor([ s[2] ,e( [7] )])
+tensor([s([0]),e( [7, 1] )])
+tensor([ s[1] ,e( [8] )])
),
Partition( [2, 2, 1, 1, 1, 1, 1] ):(
+tensor([ s[5, 2] + s[7, 1] + s[9] ,e( [] )])
+tensor([ s[4, 2] + s[5, 1] + s[6, 1] + s[7] + s[8] ,e( [1] )])
+tensor([ s[4, 1] + s[6] ,e( [1, 1] )])
+tensor([ s[3, 2] + s[5] + s[5, 1] + s[7] ,e( [2] )])
+tensor([ s[3, 1] + s[4] + s[5] ,e( [2, 1] )])
+tensor([ s[3] ,e( [2, 2] )])
+tensor([ s[2, 2] + s[4, 1] + s[6] ,e( [3] )])
+tensor([ s[2, 1] + s[4] ,e( [3, 1] )])
+tensor([ s[2] ,e( [3, 2] )])
+tensor([ s[3, 1] + s[5] ,e( [4] )])
+tensor([ s[1, 1] + s[3] ,e( [4, 1] )])
+tensor([ s[1] ,e( [4, 2] )])
+tensor([ -s[1, 1] + s[2, 1] + s[4] ,e( [5] )])
+tensor([ s[2] ,e( [5, 1] )])
+tensor([s([0]),e( [5, 2] )])
+tensor([ s[1, 1] + s[3] ,e( [6] )])
+tensor([ s[1] ,e( [6, 1] )])
+tensor([ s[2] ,e( [7] )])
),
Partition( [2, 2, 2, 1, 1, 1] ):(
+tensor([ s[3, 3] + s[5, 2] + s[7, 1] + s[9] ,e( [] )])
+tensor([ s[3, 2] + s[4, 2] + s[5, 1] + s[6, 1] + s[7] + s[8] ,e( [1] )])
+tensor([ s[2, 2] + s[4, 1] + s[6] ,e( [1, 1] )])
+tensor([ -s[2, 2] + s[3, 1] + s[3, 2] + s[5] + s[5, 1] + s[7] ,e( [2] )])
+tensor([ s[2, 1] + s[3, 1] + s[4] + s[5] ,e( [2, 1] )])
+tensor([ s[1, 1] + s[3] ,e( [2, 2] )])
+tensor([ -s[2, 1] + s[2, 2] + s[3] + s[4, 1] + s[6] ,e( [3] )])
+tensor([ -s[1, 1] + s[2] + s[2, 1] + s[4] ,e( [3, 1] )])
+tensor([ s[1] + s[2] ,e( [3, 2] )])
+tensor([s([0]),e( [3, 3] )])
+tensor([ s[3, 1] + s[5] ,e( [4] )])
+tensor([ s[1, 1] + s[3] ,e( [4, 1] )])
+tensor([ s[1] ,e( [4, 2] )])
+tensor([ s[2, 1] + s[4] ,e( [5] )])
+tensor([ s[2] ,e( [5, 1] )])
+tensor([ s[3] ,e( [6] )])
),
Partition( [3, 2, 2, 1, 1] ):(
+tensor([ s[2, 2, 1] + s[4, 1, 1] + s[4, 2] + s[5, 2] + s[6, 1] + s[7, 1] + s[9] ,e( [] )])
+tensor([ s[2, 1, 1] + s[2, 2] + s[3, 1, 1] + s[3, 2] + 2*s[4, 1] + s[4, 2] + 2*s[5, 1] + s[6] + s[6, 1] + s[7] + s[8] ,e( [1] )])
+tensor([ s[1, 1, 1] + s[2, 1] + 2*s[3, 1] + s[4] + s[4, 1] + s[5] + s[6] ,e( [1, 1] )])
+tensor([ s[1, 1] + s[3] ,e( [1, 1, 1] )])
+tensor([ -s[1, 1, 1] + s[2, 1, 1] + s[2, 2] + s[3, 2] + s[4, 1] + s[5] + s[5, 1] + s[7] ,e( [2] )])
+tensor([ -s[1, 1] + s[2] + 2*s[2, 1] + s[3, 1] + 2*s[4] + s[5] ,e( [2, 1] )])
+tensor([ s[1] + s[2] ,e( [2, 1, 1] )])
+tensor([ s[3] ,e( [2, 2] )])
+tensor([s([0]),e( [2, 2, 1] )])
+tensor([ s[1, 1] + s[1, 1, 1] + s[2, 2] + s[3] + s[3, 1] + s[4, 1] + s[6] ,e( [3] )])
+tensor([ 2*s[1, 1] + s[2] + s[2, 1] + s[3] + s[4] ,e( [3, 1] )])
+tensor([ s[1] ,e( [3, 1, 1] )])
+tensor([ s[1] + s[2] ,e( [3, 2] )])
+tensor([ -s[1, 1] + s[2, 1] + s[3, 1] + s[5] ,e( [4] )])
+tensor([ s[1, 1] + s[2] + s[3] ,e( [4, 1] )])
+tensor([ s[2, 1] + s[4] ,e( [5] )])
),
Partition( [3, 3, 2, 1] ):(
+tensor([ s[3, 3] + s[4, 1, 1] + s[5, 2] + s[6, 1] + s[7, 1] + s[9] ,e( [] )])
+tensor([ s[3, 1, 1] + s[3, 2] + s[4, 1] + s[4, 2] + 2*s[5, 1] + s[6] + s[6, 1] + s[7] + s[8] ,e( [1] )])
+tensor([ s[2, 2] + s[3, 1] + s[4] + s[4, 1] + s[5] + s[6] ,e( [1, 1] )])
+tensor([ s[3] ,e( [1, 1, 1] )])
+tensor([ s[2, 1] + s[2, 1, 1] - s[2, 2] + s[3] + s[3, 1] + s[3, 2] + s[4, 1] + s[5] + s[5, 1] + s[7] ,e( [2] )])
+tensor([ s[1] + s[1, 1] + 2*s[2] + 2*s[2, 1] + s[3, 1] + 2*s[4] + s[5] ,e( [2, 1] )])
+tensor([ s[0] + s[1] + s[2] ,e( [2, 1, 1] )])
+tensor([ s[1] + s[1, 1] + s[3] ,e( [2, 2] )])
+tensor([ s[1, 1] + s[1, 1, 1] + s[2, 2] + s[3] + s[3, 1] + s[4, 1] + s[6] ,e( [3] )])
+tensor([ s[2] + s[2, 1] + s[3] + s[4] ,e( [3, 1] )])
+tensor([ s[3, 1] + s[5] ,e( [4] )])
),
Partition( [4, 3, 2] ):(
+tensor([ s[3, 3] + s[4, 1, 1] + s[5, 2] + s[6, 1] + s[7, 1] + s[9] ,e( [] )])
+tensor([ s[3, 1] + s[3, 1, 1] + s[3, 2] + s[4, 1] + s[4, 2] + s[5] + 2*s[5, 1] + s[6] + s[6, 1] + s[7] + s[8] ,e( [1] )])
+tensor([ s[2] + s[2, 1] + s[2, 2] + s[3] + s[3, 1] + 2*s[4] + s[4, 1] + s[5] + s[6] ,e( [1, 1] )])
+tensor([ s[0] + s[1] + s[2] + s[3] ,e( [1, 1, 1] )])
+tensor([ s[2, 1] + s[2, 1, 1] + s[3] + s[3, 1] + s[3, 2] + s[4, 1] + s[5] + s[5, 1] + s[7] ,e( [2] )])
+tensor([ 2*s[1] + 2*s[1, 1] + 2*s[2] + 2*s[2, 1] + s[3] + s[3, 1] + 2*s[4] + s[5] ,e( [2, 1] )])
+tensor([ s[1, 1] + s[1, 1, 1] + s[3] + s[3, 1] + s[4, 1] + s[6] ,e( [3] )])
),
Partition( [5, 3, 1] ):(
+tensor([ s[2, 2, 1] + s[4, 1, 1] + s[4, 2] + s[5, 2] + s[6, 1] + s[7, 1] + s[9] ,e( [] )])
+tensor([ s[2, 1] + s[2, 1, 1] + s[2, 2] + s[3, 1] + s[3, 1, 1] + s[3, 2] + s[4] + 2*s[4, 1] + s[4, 2] + s[5] + 2*s[5, 1] + s[6] + s[6, 1] + s[7] + s[8] ,e( [1] )])
+tensor([ s[1] + 2*s[1, 1] + s[1, 1, 1] + 2*s[2] + 2*s[2, 1] + 2*s[3] + 2*s[3, 1] + 2*s[4] + s[4, 1] + s[5] + s[6] ,e( [1, 1] )])
+tensor([ s[0] + 2*s[1] + s[1, 1] + s[2] + s[3] ,e( [1, 1, 1] )])
+tensor([ s[2, 1] + s[2, 1, 1] + s[2, 2] + s[3] + s[3, 1] + s[3, 2] + s[4, 1] + s[5] + s[5, 1] + s[7] ,e( [2] )])
+tensor([ 2*s[2] + 2*s[2, 1] + s[3] + s[3, 1] + 2*s[4] + s[5] ,e( [2, 1] )])
+tensor([ s[2, 2] + s[4, 1] + s[6] ,e( [3] )])
),
Partition( [6, 3] ):(
+tensor([ s[3, 3] + s[5, 2] + s[7, 1] + s[9] ,e( [] )])
+tensor([ s[3] + s[3, 1] + s[3, 2] + s[4] + s[4, 1] + s[4, 2] + s[5] + s[5, 1] + s[6] + s[6, 1] + s[7] + s[8] ,e( [1] )])
+tensor([ s[0] + 2*s[1] + s[1, 1] + 3*s[2] + 2*s[2, 1] + s[2, 2] + 2*s[3] + s[3, 1] + 2*s[4] + s[4, 1] + s[5] + s[6] ,e( [1, 1] )])
+tensor([ s[3] + s[3, 1] + s[3, 2] + s[5] + s[5, 1] + s[7] ,e( [2] )])
),
Partition( [7, 2] ):(
+tensor([ s[5, 2] + s[7, 1] + s[9] ,e( [] )])
+tensor([ s[2] + s[2, 1] + s[2, 2] + s[3] + s[3, 1] + s[3, 2] + s[4] + s[4, 1] + s[4, 2] + s[5] + s[5, 1] + s[6] + s[6, 1] + s[7] + s[8] ,e( [1] )])
+tensor([ s[0] + 2*s[1] + s[1, 1] + 2*s[2] + s[2, 1] + 2*s[3] + s[3, 1] + 2*s[4] + s[4, 1] + s[5] + s[6] ,e( [1, 1] )])
+tensor([ s[5] + s[5, 1] + s[7] ,e( [2] )])
),
Partition( [8, 1] ):(
+tensor([ s[7, 1] + s[9] ,e( [] )])
+tensor([ s[1] + s[1, 1] + s[2] + s[2, 1] + s[3] + s[3, 1] + s[4] + s[4, 1] + s[5] + s[5, 1] + s[6] + s[6, 1] + s[7] + s[8] ,e( [1] )])
+tensor([ s[0] + s[1] + s[2] + s[3] + s[4] + s[5] + s[6] ,e( [1, 1] )])
+tensor([ s[7] ,e( [2] )])
),
Partition( [9] ):(
+tensor([ s[9] ,e( [] )])
+tensor([ s[0] + s[1] + s[2] + s[3] + s[4] + s[5] + s[6] + s[7] + s[8] ,e( [1] )])
),
Partition( [1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ):(
+tensor([ s[10] ,e( [] )])
+tensor([ s[9] ,e( [1] )])
+tensor([ s[8] ,e( [2] )])
+tensor([ s[7] ,e( [3] )])
+tensor([ s[6] ,e( [4] )])
+tensor([ s[5] ,e( [5] )])
+tensor([ s[4] ,e( [6] )])
+tensor([ s[3] ,e( [7] )])
+tensor([ s[2] ,e( [8] )])
+tensor([ s[1] ,e( [9] )])
+tensor([s([0]),e( [10] )])
),
Partition( [2, 1, 1, 1, 1, 1, 1, 1, 1] ):(
+tensor([ s[8, 1] + s[10] ,e( [] )])
+tensor([ s[7, 1] + s[8] + s[9] ,e( [1] )])
+tensor([ s[7] ,e( [1, 1] )])
+tensor([ s[6, 1] + s[8] ,e( [2] )])
+tensor([ s[6] ,e( [2, 1] )])
+tensor([ s[5, 1] + s[7] ,e( [3] )])
+tensor([ s[5] ,e( [3, 1] )])
+tensor([ s[4, 1] + s[6] ,e( [4] )])
+tensor([ s[4] ,e( [4, 1] )])
+tensor([ s[3, 1] + s[5] ,e( [5] )])
+tensor([ s[3] ,e( [5, 1] )])
+tensor([ s[2, 1] + s[4] ,e( [6] )])
+tensor([ s[2] ,e( [6, 1] )])
+tensor([ s[1, 1] + s[3] ,e( [7] )])
+tensor([ s[1] ,e( [7, 1] )])
+tensor([ s[2] ,e( [8] )])
+tensor([s([0]),e( [8, 1] )])
+tensor([ s[1] ,e( [9] )])
),
Partition( [2, 2, 1, 1, 1, 1, 1, 1] ):(
+tensor([ s[6, 2] + s[8, 1] + s[10] ,e( [] )])
+tensor([ s[5, 2] + s[6, 1] + s[7, 1] + s[8] + s[9] ,e( [1] )])
+tensor([ s[5, 1] + s[7] ,e( [1, 1] )])
+tensor([ s[4, 2] + s[6] + s[6, 1] + s[8] ,e( [2] )])
+tensor([ s[4, 1] + s[5] + s[6] ,e( [2, 1] )])
+tensor([ s[4] ,e( [2, 2] )])
+tensor([ s[3, 2] + s[5, 1] + s[7] ,e( [3] )])
+tensor([ s[3, 1] + s[5] ,e( [3, 1] )])
+tensor([ s[3] ,e( [3, 2] )])
+tensor([ s[2, 2] + s[4, 1] + s[6] ,e( [4] )])
+tensor([ s[2, 1] + s[4] ,e( [4, 1] )])
+tensor([ s[2] ,e( [4, 2] )])
+tensor([ s[3, 1] + s[5] ,e( [5] )])
+tensor([ s[1, 1] + s[3] ,e( [5, 1] )])
+tensor([ s[1] ,e( [5, 2] )])
+tensor([ -s[1, 1] + s[2, 1] + s[4] ,e( [6] )])
+tensor([ s[2] ,e( [6, 1] )])
+tensor([s([0]),e( [6, 2] )])
+tensor([ s[1, 1] + s[3] ,e( [7] )])
+tensor([ s[1] ,e( [7, 1] )])
+tensor([ s[2] ,e( [8] )])
),
Partition( [2, 2, 2, 1, 1, 1, 1] ):(
+tensor([ s[4, 3] + s[6, 2] + s[8, 1] + s[10] ,e( [] )])
+tensor([ s[3, 3] + s[4, 2] + s[5, 2] + s[6, 1] + s[7, 1] + s[8] + s[9] ,e( [1] )])
+tensor([ s[3, 2] + s[5, 1] + s[7] ,e( [1, 1] )])
+tensor([ s[4, 1] + s[4, 2] + s[6] + s[6, 1] + s[8] ,e( [2] )])
+tensor([ s[2, 2] + s[3, 1] + s[4, 1] + s[5] + s[6] ,e( [2, 1] )])
+tensor([ s[2, 1] + s[4] ,e( [2, 2] )])
+tensor([ -s[2, 2] + s[3, 2] + s[4] + s[5, 1] + s[7] ,e( [3] )])
+tensor([ s[3] + s[3, 1] + s[5] ,e( [3, 1] )])
+tensor([ s[1, 1] + s[2] + s[3] ,e( [3, 2] )])
+tensor([ s[1] ,e( [3, 3] )])
+tensor([ -s[2, 1] + s[2, 2] + s[4, 1] + s[6] ,e( [4] )])
+tensor([ -s[1, 1] + s[2, 1] + s[4] ,e( [4, 1] )])
+tensor([ s[2] ,e( [4, 2] )])
+tensor([s([0]),e( [4, 3] )])
+tensor([ s[3, 1] + s[5] ,e( [5] )])
+tensor([ s[1, 1] + s[3] ,e( [5, 1] )])
+tensor([ s[1] ,e( [5, 2] )])
+tensor([ s[2, 1] + s[4] ,e( [6] )])
+tensor([ s[2] ,e( [6, 1] )])
+tensor([ s[3] ,e( [7] )])
),
Partition( [3, 2, 2, 1, 1, 1] ):(
+tensor([ s[3, 2, 1] + s[4, 3] + s[5, 1, 1] + s[5, 2] + s[6, 2] + s[7, 1] + s[8, 1] + s[10] ,e( [] )])
+tensor([ s[2, 2, 1] + s[3, 1, 1] + s[3, 2] + s[3, 3] + s[4, 1, 1] + 2*s[4, 2] + 2*s[5, 1] + s[5, 2] + 2*s[6, 1] + s[7] + s[7, 1] + s[8] + s[9] ,e( [1] )])
+tensor([ s[2, 1, 1] + s[2, 2] + s[3, 1] + s[3, 2] + 2*s[4, 1] + s[5] + s[5, 1] + s[6] + s[7] ,e( [1, 1] )])
+tensor([ s[2, 1] + s[4] ,e( [1, 1, 1] )])
+tensor([ s[3, 1, 1] + s[3, 2] + s[4, 1] + s[4, 2] + s[5, 1] + s[6] + s[6, 1] + s[8] ,e( [2] )])
+tensor([ s[1, 1, 1] + s[2, 2] + s[3] + 3*s[3, 1] + s[4, 1] + 2*s[5] + s[6] ,e( [2, 1] )])
+tensor([ s[1, 1] + s[2] + s[3] ,e( [2, 1, 1] )])
+tensor([ s[2, 1] + s[4] ,e( [2, 2] )])
+tensor([ s[1] ,e( [2, 2, 1] )])
+tensor([ -s[1, 1, 1] + s[2, 1, 1] + s[3, 2] + s[4] + s[4, 1] + s[5, 1] + s[7] ,e( [3] )])
+tensor([ -s[1, 1] + 2*s[2, 1] + s[3] + s[3, 1] + s[4] + s[5] ,e( [3, 1] )])
+tensor([ s[2] ,e( [3, 1, 1] )])
+tensor([ s[1, 1] + s[2] + s[3] ,e( [3, 2] )])
+tensor([s([0]),e( [3, 2, 1] )])
+tensor([ s[1] ,e( [3, 3] )])
+tensor([ s[1, 1] + s[1, 1, 1] - s[2, 1] + s[2, 2] + s[3, 1] + s[4, 1] + s[6] ,e( [4] )])
+tensor([ s[1, 1] + s[2, 1] + s[3] + s[4] ,e( [4, 1] )])
+tensor([ s[1] ,e( [4, 1, 1] )])
+tensor([ s[2] ,e( [4, 2] )])
+tensor([ -s[1, 1] + s[2, 1] + s[3, 1] + s[5] ,e( [5] )])
+tensor([ s[1, 1] + s[2] + s[3] ,e( [5, 1] )])
+tensor([ s[2, 1] + s[4] ,e( [6] )])
),
Partition( [3, 3, 2, 1, 1] ):(
+tensor([ s[3, 2, 1] + s[4, 3] + s[5, 1, 1] + s[5, 2] + s[6, 2] + s[7, 1] + s[8, 1] + s[10] ,e( [] )])
+tensor([ s[2, 2, 1] + s[3, 1, 1] + s[3, 2] + s[3, 3] + s[4, 1, 1] + 2*s[4, 2] + 2*s[5, 1] + s[5, 2] + 2*s[6, 1] + s[7] + s[7, 1] + s[8] + s[9] ,e( [1] )])
+tensor([ s[2, 1, 1] + s[2, 2] + s[3, 1] + s[3, 2] + 2*s[4, 1] + s[5] + s[5, 1] + s[6] + s[7] ,e( [1, 1] )])
+tensor([ s[2, 1] + s[4] ,e( [1, 1, 1] )])
+tensor([ s[2, 1] + s[3, 1, 1] + s[3, 2] + s[4] + s[4, 1] + s[4, 2] + s[5, 1] + s[6] + s[6, 1] + s[8] ,e( [2] )])
+tensor([ s[1, 1] + s[1, 1, 1] + s[2] + s[2, 2] + 2*s[3] + 3*s[3, 1] + s[4, 1] + 2*s[5] + s[6] ,e( [2, 1] )])
+tensor([ s[1] + s[1, 1] + s[2] + s[3] ,e( [2, 1, 1] )])
+tensor([ s[2] + s[2, 1] + s[4] ,e( [2, 2] )])
+tensor([ s[0] + s[1] ,e( [2, 2, 1] )])
+tensor([ -s[1, 1] - s[1, 1, 1] + s[2, 1] + s[2, 1, 1] + s[3, 2] + s[4] + s[4, 1] + s[5, 1] + s[7] ,e( [3] )])
+tensor([ 2*s[2, 1] + s[3] + s[3, 1] + s[4] + s[5] ,e( [3, 1] )])
+tensor([ s[2] ,e( [3, 1, 1] )])
+tensor([ s[1] + s[1, 1] + s[2] + s[3] ,e( [3, 2] )])
+tensor([ s[1, 1] + s[1, 1, 1] - s[2, 1] + s[2, 2] + s[3, 1] + s[4, 1] + s[6] ,e( [4] )])
+tensor([ s[2, 1] + s[3] + s[4] ,e( [4, 1] )])
+tensor([ s[3, 1] + s[5] ,e( [5] )])
),
Partition( [4, 3, 2, 1] ):(
+tensor([ s[1, 1, 1, 1] + s[3, 1, 1] + s[4, 1, 1] + s[4, 2] + s[4, 3] + s[5, 1, 1] + s[6, 1] + s[6, 2] + s[7, 1] + s[8, 1] + s[10] ,e( [] )])
+tensor([ s[1, 1, 1] + s[2, 1, 1] + s[3, 1] + s[3, 1, 1] + s[3, 2] + s[3, 3] + 2*s[4, 1] + s[4, 1, 1] + s[4, 2] + 2*s[5, 1] + s[5, 2] + s[6] + 2*s[6, 1] + s[7] + s[7, 1] + s[8] + s[9] ,e( [1] )])
+tensor([ s[1, 1] + s[2, 1] + s[3] + 2*s[3, 1] + s[3, 2] + s[4] + s[4, 1] + 2*s[5] + s[5, 1] + s[6] + s[7] ,e( [1, 1] )])
+tensor([ s[1] + s[2] + s[3] + s[4] ,e( [1, 1, 1] )])
+tensor([s([0]),e( [1, 1, 1, 1] )])
+tensor([ s[1, 1, 1] + s[2, 1] + s[2, 1, 1] + s[2, 2] + s[3, 1] + s[3, 1, 1] + s[4] + 2*s[4, 1] + s[4, 2] + s[5, 1] + s[6] + s[6, 1] + s[8] ,e( [2] )])
+tensor([ 2*s[1, 1] + 2*s[2] + 3*s[2, 1] + s[2, 2] + 2*s[3] + 2*s[3, 1] + s[4] + s[4, 1] + 2*s[5] + s[6] ,e( [2, 1] )])
+tensor([ 3*s[1] + 2*s[2] + s[3] ,e( [2, 1, 1] )])
+tensor([ s[1, 1] + s[2] + s[2, 1] + s[4] ,e( [2, 2] )])
+tensor([ s[1, 1, 1] + s[2, 1] + s[2, 1, 1] + s[3, 1] + s[3, 2] + s[4] + s[4, 1] + s[5, 1] + s[7] ,e( [3] )])
+tensor([ 2*s[1, 1] + s[2, 1] + 2*s[3] + s[3, 1] + s[4] + s[5] ,e( [3, 1] )])
+tensor([ s[1, 1, 1] + s[3, 1] + s[4, 1] + s[6] ,e( [4] )])
),
Partition( [5, 3, 2] ):(
+tensor([ s[3, 2, 1] + s[4, 3] + s[5, 1, 1] + s[5, 2] + s[6, 2] + s[7, 1] + s[8, 1] + s[10] ,e( [] )])
+tensor([ s[2, 1] + s[2, 1, 1] + s[2, 2] + s[2, 2, 1] + s[3, 1] + s[3, 1, 1] + s[3, 2] + s[3, 3] + s[4, 1] + s[4, 1, 1] + 2*s[4, 2] + s[5] + 2*s[5, 1] + s[5, 2] + s[6] + 2*s[6, 1] + s[7] + s[7, 1] + s[8] + s[9] ,e( [1] )])
+tensor([ s[1, 1] + s[1, 1, 1] + s[2] + 2*s[2, 1] + s[2, 1, 1] + s[2, 2] + 2*s[3] + 2*s[3, 1] + s[3, 2] + 2*s[4] + 2*s[4, 1] + 2*s[5] + s[5, 1] + s[6] + s[7] ,e( [1, 1] )])
+tensor([ s[0] + 2*s[1] + s[1, 1] + 2*s[2] + s[2, 1] + s[3] + s[4] ,e( [1, 1, 1] )])
+tensor([ s[3, 1] + s[3, 1, 1] +s[3, 2] + s[4] + s[4, 1] + s[4, 2] + s[5, 1] + s[6] + s[6, 1] + s[8] ,e( [2] )])
+tensor([ s[1] + s[1, 1] + 2*s[2] + 2*s[2, 1] + s[2, 2] + 3*s[3] + 3*s[3, 1] + s[4] + s[4, 1] + 2*s[5] + s[6] ,e( [2, 1] )])
+tensor([ s[3, 2] + s[4] + s[4, 1] + s[5, 1] + s[7] ,e( [3] )])
),
Partition( [6, 3, 1] ):(
+tensor([ s[3, 2, 1] + s[4, 3] + s[5, 1, 1] + s[5, 2] + s[6, 2] + s[7, 1] + s[8, 1] + s[10] ,e( [] )])
+tensor([ s[2, 1] + s[2, 1, 1] + s[2, 2] + s[2, 2, 1] + s[3, 1] + s[3, 1, 1] + s[3, 2] + s[3, 3] + s[4] + 2*s[4, 1] + s[4, 1, 1] + 2*s[4, 2] + s[5] + 2*s[5, 1] + s[5, 2] + s[6] + 2*s[6, 1] + s[7] + s[7, 1] + s[8] + s[9] ,e( [1] )])
+tensor([ s[1] + 2*s[1, 1] + s[1, 1, 1] + 2*s[2] + 3*s[2, 1] + s[2, 1, 1] + s[2, 2] + 3*s[3] + 3*s[3, 1] + s[3, 2] + 2*s[4] + 2*s[4, 1] + 2*s[5] + s[5, 1] + s[6] + s[7] ,e( [1, 1] )])
+tensor([ s[0] + 2*s[1] + s[1, 1] + 2*s[2] + s[2, 1] + s[3] + s[4] ,e( [1, 1, 1] )])
+tensor([ s[3, 1] + s[3, 1, 1] + s[3, 2] + s[4] + s[4, 1] + s[4, 2] + s[5, 1] + s[6] + s[6, 1] + s[8] ,e( [2] )])
+tensor([ s[2] + s[2, 1] + s[2, 2] + 2*s[3] + 2*s[3, 1] + s[4] + s[4, 1] + 2*s[5] + s[6] ,e( [2, 1] )])
+tensor([ s[3, 2] + s[5, 1] + s[7] ,e( [3] )])
),
Partition( [7, 3] ):(
+tensor([ s[4, 3] + s[6, 2] + s[8, 1] + s[10] ,e( [] )])
+tensor([ s[3] + s[3, 1] + s[3, 2] + s[3, 3] + s[4] + s[4, 1] + s[4, 2] + s[5] + s[5, 1] + s[5, 2] + s[6] + s[6, 1] + s[7] + s[7, 1] + s[8] + s[9] ,e( [1] )])
+tensor([ s[0] + 2*s[1] + s[1, 1] + 3*s[2] + 2*s[2, 1] + s[2, 2] + 3*s[3] + 2*s[3, 1] + s[3, 2] + 2*s[4] + s[4, 1] + 2*s[5] + s[5, 1] + s[6] + s[7] ,e( [1, 1] )])
+tensor([ s[4] + s[4, 1] + s[4, 2] + s[6] + s[6, 1] + s[8] ,e( [2] )])
),
Partition( [8, 2] ):(
+tensor([ s[6, 2] + s[8, 1] + s[10] ,e( [] )])
+tensor([ s[2] + s[2, 1] + s[2, 2] + s[3] + s[3, 1] + s[3, 2] + s[4] + s[4, 1] + s[4, 2] + s[5] + s[5, 1] + s[5, 2] + s[6] + s[6, 1] + s[7] + s[7, 1] + s[8] + s[9] ,e( [1] )])
+tensor([ s[0] + 2*s[1] + s[1, 1] + 2*s[2] + s[2, 1] + 2*s[3] + s[3, 1] + 2*s[4] + s[4, 1] + 2*s[5] + s[5, 1] + s[6] + s[7] ,e( [1, 1] )])
+tensor([ s[6] + s[6, 1] + s[8] ,e( [2] )])
),
Partition( [9, 1] ):(
+tensor([ s[8, 1] + s[10] ,e( [] )])
+tensor([ s[1] + s[1, 1] + s[2] + s[2, 1] + s[3] + s[3, 1] + s[4] + s[4, 1] + s[5] + s[5, 1] + s[6] + s[6, 1] + s[7] + s[7, 1] + s[8] + s[9] ,e( [1] )])
+tensor([ s[0] + s[1] + s[2] + s[3] + s[4] + s[5] + s[6] + s[7] ,e( [1, 1] )])
+tensor([ s[8] ,e( [2] )])
),
Partition( [10] ):(
+tensor([ s[10] ,e( [] )])
+tensor([ s[0] + s[1] + s[2] + s[3] + s[4] + s[5] + s[6] + s[7] + s[8] + s[9] ,e( [1] )])
),
})