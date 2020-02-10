# BrokenLineShellings

To get started using broken line shelling orders to produce new posets recording the h-vector of a matroid, we recommend using Cocalc.com, as this can used from any internet browser without any installation. To get started with Cocalc [go here](https://doc.cocalc.com/getting-started.html).

After logging into Cocalc, simply download the file `shellings.sage` and then upload it to your current Cocalc project folder.

Open a new sage worksheet in Cocalc and type `load("shellings.sage")` as your first command. After that you will have access to all the functions defined in the `shellings.sage` file. You may then run commands as described in the paper, for example:

```
load("shellings.sage")
V = [(1,1,1,0,0,0),
    (1,1,0,1,0,0),
    (1,1,0,0,1,0),
    (1,1,0,0,0,1),
    (1,0,1,1,0,0),
    (1,0,1,0,1,0),
    (1,0,1,0,0,1),
    (1,0,0,1,1,0),
    (1,0,0,1,0,1),
    (0,1,1,1,0,0),
    (0,1,1,0,1,0),
    (0,1,1,0,0,1),
    (0,1,0,1,1,0),
    (0,1,0,1,0,1)]
vFav = (1,1,1,0,0,0)
pivots = [3,6,7]
misses = 5
limit = 3
w = 5.0
result = search(V, vFav, pivots, limit, misses, w)
display_results(result)
```

One may also use the built-in matroids in `SAGE`, as in:

```
U24 = matroids.Uniform(2,4)
P = U24.matroid_polytope()
V = P.vertices_list()
```
This is equivalent to directly typing
```
V = [[0, 0, 1, 1], [0, 1, 0, 1], [0, 1, 1, 0],
          [1, 0, 0, 1], [1, 0, 1, 0], [1, 1, 0, 0]]
```

Finally, we note that as you update your search parameters, you should try to keep the number of pivots to 2, at least initially, as having 3 or more pivots will quickly cause longer computations. As you adjust the other parameters, or adjust where the pivots are, you will want to use the function `update_search` which takes the same arguments as `search` except the first positional argument should be the `result` of a previous computation. This way you do not throw away the posets obtained from previous searches, and the code will automatically compare each of the new posets to your old posets and each other, giving you an updated list of non-isomorphic posets obtained by the broken line shellings you have explored so far. An example of `update_search` is as follows, where we are assuming you have used `search(args...)` to obtain a `result`, which you use as input to `update_search` as in

```
result = update_search(result, V, vFav, pivots, limit, misses, w)
display_results(result)
```
