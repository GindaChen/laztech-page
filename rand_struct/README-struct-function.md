# How to Calculate Structure Function in a decent manner?

## The crux of the problem

We want to find the **structure function $SF$** of a given function $D: R^N\rightarrow R^M$.

We define structure function
$$
SF(R) = \textit{average\{}\sum_{||\vec{q}-\vec{p}|| = R} ||D(\vec{q}) - D(\vec{p})||^2\ \  \textit{\}}
$$
where the domain function
$$
D(\vec{p}) : R^N \rightarrow R^M
$$
Typically the domain is a *subspace* of $R^N$, and $N=3$ for 3D field. and $M= 1, 3$ for either a scaler field or a vector field.



The simpliest form in 3D structure function looks like this $O(N^6)$  (never wrote an ugly code like this, huh?) 

**Naive 3D structure function**

```julia
function struct_function_3D_naive(T::tensor)
    counts = []; # How many points we have count in the certain R^2
    result = []; # The accumulated array to store the result
    
    for px in 1:nx, py in 1:ny, pz in 1:nz
        for qx in 1:nx, qy in 1:ny, qz in 1:nz
            p = Point(px, py, pz) 
            q = Point(qx, qy, qz)
            
            R = distance3D(p, q) 
            R2 = R * R;
            
            dT = norm(T[q] - T[p]) .^ 2 ;
            result[R2] += dT; # magic: we use R2 to preserve the integer index 
            counts[R2] += 1; 
        end
    end
    
    return result ./ counts # The average in each of the bucket
end
```



## The problem: speed & accuracy

In our field of study, a **cube** is typically something 3D and $N = 480, 792, 1024, 2048, ...$  

For a 3D cube, the complexity of the problem become $O(N^6)$. Wow... you think, what should we do?



## Solution 

### 1. Monte Carlo

What we have done for a simple case is to randomly choose points so taht we can obtain the suboptimal solution of the result. There are several following problems to consider ~~(parameters to take)~~:

1. How many **seed-points** are **statistically sufficient** to do the Monte Carlo?  More specifically, what is the distribution of the points (random without replacement? points should have certain distances? fulfill some constrain?)
2. For a certain seed-point $p$, how many **other points** $q$ should we choose? What is the distribution (preferably related to the distance between $p, q$ , say the larger the $R$ the more we should take to stablize the region we want -- say, when what we really want is to see the region when $R$ is larger). 

Simply speaking:

1. How big is `nRand` when choosing a cube `D` with size `nx, ny, nz` (or simply `N^3`)
2. Given the randomly chosen seed-points set `P`, for each point `p`,  how to choose points `q` ?



### 2. Parallel Programming

Parallel does not solve for everything... But for smaller cube, it is hopeful that parallelizing the problem will get us to the accurate solution in the real world (especially for those who love astronomically-accurate result). 

Then the problem becomes:

1. How to dispatch the tasks? Typically when the cube is large, it is not possible to dispatch the whole cube to each process. A crude estimation of a $792^3$ cube will take up several GB memory (when we are going to use the information we want, say $B_x, B_y, B_z$, all double precision and has $792^3$ items, we will get approximately 1G for each cube). How to elegantly distribute the problem becomes a problem.
2. How to prevent bugs. With multi-thread/multi-process, we will get a headache wiht synchronization, race condition, etc.  We could also avoid the problem by naively generate the results, stored in local environment, and then merge them later on.



## Progress

