# BrokenLineShellings

To get started using broken line shelling orders to produce new posets recording the h-vector of a matroid, we recommend using Cocalc.com, as this can used from any internet browser without any installation. To get started with Cocalc [go here](https://doc.cocalc.com/getting-started.html).

After logging into Cocalc, simply download the file `shellings.sage` and then upload it to your current Cocalc project folder.

Open a new sage worksheet in Cocalc and type `load("shellings.sage")` as your first command. After that you will have access to all the functions defined in the `shellings.sage` file. You may then run commands as described in the paper, for example:

```python
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
