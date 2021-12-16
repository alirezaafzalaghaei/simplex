

# A bare bone simplex algorithm
In this repo, a basic simplex method visualizer/solver has been implemented. Unlike common Linear Programming solvers, this repo focuses on describing **tableau format of the simplex algorithm**. The implementation and examples are based on the well-known Bazaraa's book:
> Bazaraa, M.S., Jarvis, J.J. and Sherali, H.D., 2008. Linear Programming and Network Flows. John Wiley & Sons.

## Features
 - The solver prints the tableau for each simplex iteration.
 - The solver uses `Rational` numbers instead of python's `float` data type.
 - It supports both Two-Phase and the Big M methods.
 - The Big M method solves the problem symbolically (using sympy)
 - The solver can handle unbounded problems
 - The solver can detect the cycling of the basis.
## Requirements

- Python
- [Modified gilp package](https://github.com/alirezaafzalaghaei/gilp), clone and then install it using `pip install .`
- matplotlib
- seaborn
- pylab
- numpy
- sympy
- scipy
- tabulate
- gilp requiremetns


## Visualization
Plotting Basic Feasible Solutions (BFS) and Basic Infeasible Solutions (BIS) are done in `BFS.ipynb`.

## Simplex

- `Simplex.ipynb`: Solving `Example 3.9` of the book.
- `TwoPhase.ipynb`: Solving `Example 4.3` of the book.
- `BigM.ipynb`: Solving `Example 4.6` of the book.
- `Unbounded.ipynb`: Solving `Example 3.7` of the book.
- `CyclingDetection.ipynb`: Solving `Example 4.11` of the book.
- `simplex.py`: The module!

## Example

Solve the Example 3.9 of Bazaraa's book:

<img src="https://render.githubusercontent.com/render/math?math=%5Cbegin%7Barray%7D%7Bc%7D%0A%5Cdisplaystyle%5Cmin%20%26x_1%2Bx_2-4x_3%26%5C%5C%0As.t.%26%20x_1%2Bx_2%2B2x_3%20%26%5Cle%209%5C%5C%0A%26x_1%2Bx_2-x_3%20%26%5Cle%202%5C%5C%0A%26-x_1%2Bx_2%2Bx_3%20%26%5Cle%204%5C%5C%0A%26x_1%20%2C%20x_2%2Cx_3%20%26%5Cge%200%0A%5Cend%7Barray%7D">

- program.py:
```
from simplex import solve
LP = """
min 1 1 -4
1 1 2 <= 9
1 1 -1 <= 2
-1 1 1 <= 4
"""
z, A, b, c, x, B, Btrace = solve(LP)
```
- Output:
```
max iterations: 20
1th iteration:
╒═════╤═════╤═══════╤═══════╤═══════╤═══════╤═══════╤═══════╤═══════╕
│     │  z  │  x_1  │  x_2  │  x_3  │  x_4  │  x_5  │  x_6  │  RHS  │
╞═════╪═════╪═══════╪═══════╪═══════╪═══════╪═══════╪═══════╪═══════╡
│  z  │  1  │  -1   │  -1   │   4   │   0   │   0   │   0   │   0   │
├─────┼─────┼───────┼───────┼───────┼───────┼───────┼───────┼───────┤
│ x_4 │  0  │   1   │   1   │   2   │   1   │   0   │   0   │   9   │
├─────┼─────┼───────┼───────┼───────┼───────┼───────┼───────┼───────┤
│ x_5 │  0  │   1   │   1   │  -1   │   0   │   1   │   0   │   2   │
├─────┼─────┼───────┼───────┼───────┼───────┼───────┼───────┼───────┤
│ x_6 │  0  │  -1   │   1   │   1   │   0   │   0   │   1   │   4   │
╘═════╧═════╧═══════╧═══════╧═══════╧═══════╧═══════╧═══════╧═══════╛
Beginning Problem solving:
2th iteration:
Entering variable: x_3
Exiting  variable: x_6
Pivot: 1
╒═════╤═════╤═══════╤═══════╤═══════╤═══════╤═══════╤═══════╤═══════╕
│     │  z  │  x_1  │  x_2  │  x_3  │  x_4  │  x_5  │  x_6  │  RHS  │
╞═════╪═════╪═══════╪═══════╪═══════╪═══════╪═══════╪═══════╪═══════╡
│  z  │  1  │   3   │  -5   │   0   │   0   │   0   │  -4   │  -16  │
├─────┼─────┼───────┼───────┼───────┼───────┼───────┼───────┼───────┤
│ x_4 │  0  │   3   │  -1   │   0   │   1   │   0   │  -2   │   1   │
├─────┼─────┼───────┼───────┼───────┼───────┼───────┼───────┼───────┤
│ x_5 │  0  │   0   │   2   │   0   │   0   │   1   │   1   │   6   │
├─────┼─────┼───────┼───────┼───────┼───────┼───────┼───────┼───────┤
│ x_3 │  0  │  -1   │   1   │   1   │   0   │   0   │   1   │   4   │
╘═════╧═════╧═══════╧═══════╧═══════╧═══════╧═══════╧═══════╧═══════╛
3th iteration:
Entering variable: x_1
Exiting  variable: x_4
Pivot: 3
╒═════╤═════╤═══════╤═══════╤═══════╤═══════╤═══════╤═══════╤═══════╕
│     │  z  │  x_1  │  x_2  │  x_3  │  x_4  │  x_5  │  x_6  │  RHS  │
╞═════╪═════╪═══════╪═══════╪═══════╪═══════╪═══════╪═══════╪═══════╡
│  z  │  1  │   0   │  -4   │   0   │  -1   │   0   │  -2   │  -17  │
├─────┼─────┼───────┼───────┼───────┼───────┼───────┼───────┼───────┤
│ x_1 │  0  │   1   │ -1/3  │   0   │  1/3  │   0   │ -2/3  │  1/3  │
├─────┼─────┼───────┼───────┼───────┼───────┼───────┼───────┼───────┤
│ x_5 │  0  │   0   │   2   │   0   │   0   │   1   │   1   │   6   │
├─────┼─────┼───────┼───────┼───────┼───────┼───────┼───────┼───────┤
│ x_3 │  0  │   0   │  2/3  │   1   │  1/3  │   0   │  1/3  │ 13/3  │
╘═════╧═════╧═══════╧═══════╧═══════╧═══════╧═══════╧═══════╧═══════╛
4th iteration:
Negative Zj-Cj. Optimal Solution found. Z=-17
```
